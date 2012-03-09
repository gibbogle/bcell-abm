! Concentration fields

module fields
use global
use bdry_linked_list
implicit none
save

logical, parameter :: use_S1P = .true.
logical, parameter :: use_CXCL13 = .false.
logical, parameter :: use_CCL21 = .true.
logical, parameter :: use_OXY = .false.
real, parameter :: S1P_KDIFFUSION = 0.001
real, parameter :: S1P_KDECAY = 0.000005
real, parameter :: OXY_KDIFFUSION = 0.001
real, parameter :: OXY_KDECAY = 0.000005
real, parameter :: CCL21_KDIFFUSION = 0.001
real, parameter :: CCL21_KDECAY = 0.00001
real, parameter :: CXCL13_KDIFFUSION = 1.0
real, parameter :: CXCL13_KDECAY = 1.0
integer, parameter :: MAX_BDRY = 20000

real :: BdryS1PConc = 1.0
real :: BdryOXYConc = 1.0
real :: BdryCCL21Conc = 1.0
real :: BdryCXCL13Conc = 1.0

real, allocatable :: S1P_conc(:,:,:), S1P_grad(:,:,:,:), S1P_influx(:,:,:)
real, allocatable :: CXCL13_conc(:,:,:), CXCL13_grad(:,:,:,:), CXCL13_influx(:,:,:)
real, allocatable :: CCL21_conc(:,:,:), CCL21_grad(:,:,:,:), CCL21_influx(:,:,:)
real, allocatable :: OXY_conc(:,:,:), OXY_grad(:,:,:,:), OXY_influx(:,:,:)

contains

!----------------------------------------------------------------------------------------
! Check to see if (x,y,z) is outside the grid
!----------------------------------------------------------------------------------------
logical function outside_xyz(x,y,z)
integer :: x, y, z
outside_xyz = .true.
if (x < 1 .or. x > NX) return
if (y < 1 .or. y > NY) return
if (z < 1 .or. z > NZ) return
outside_xyz = .false.
end function

!----------------------------------------------------------------------------------------
! A boundary site is inside and has at least one Neumann neighbour outside the blob.
!----------------------------------------------------------------------------------------
logical function isbdry(x,y,z)
integer :: x, y, z
integer :: k, xx, yy, zz

if (occupancy(x,y,z)%indx(1) < 0) then
    isbdry = .false.
    return
endif
do k = 1,6
	xx = x + neumann(1,k)
	yy = y + neumann(2,k)
	zz = z + neumann(3,k)
	if (outside_xyz(xx,yy,zz)) cycle
    if (occupancy(xx,yy,zz)%indx(1) < 0) then
        isbdry = .true.
        return
    endif
enddo
isbdry = .false.
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CreateBdrylist
integer :: x, y, z
integer :: k, site(3)
type (boundary_type), pointer :: bdry

!allocate(bdrylist(MAX_BDRY))
nullify(bdrylist)
nbdry = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            if (isbdry(x,y,z)) then
                nbdry = nbdry + 1
                allocate(bdry)
                site = (/x,y,z/)
                bdry%site = site
                bdry%S1P = .false.
                bdry%CXCL13 = .false.
                bdry%CCL21 = .false.
                bdry%OXY = .false.
!                nullify(bdry%previous)
                nullify(bdry%next)
                call bdrylist_insert(bdry,bdrylist)
                call AssignBdryRole(site,bdry)
!BDRY                occupancy(x,y,z)%bdry = .true.
                occupancy(x,y,z)%bdry => bdry
            endif
        enddo
    enddo
enddo
write(*,*) 'nbdry = ',nbdry, ' NBcells: ',NBcells
!call bdrylist_print(bdrylist)
!call TestAddRemoveSites
!stop
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CreateChecklist
integer :: x, y, z, ncheck
integer :: k, site(3)
type (boundary_type), pointer :: bdry

write(*,*) 'CreateCheckList'
nullify(checklist)
ncheck = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            if (isbdry(x,y,z)) then
                ncheck = ncheck + 1
                allocate(bdry)
                site = (/x,y,z/)
                bdry%site = site
                bdry%S1P = .false.
                bdry%CXCL13 = .false.
                bdry%CCL21 = .false.
                bdry%OXY = .false.
!                nullify(bdry%previous)
                nullify(bdry%next)
                call bdrylist_insert(bdry,checklist)
                call AssignBdryRole(site,bdry)
!                occupancy(x,y,z)%bdry = .true.
            endif
        enddo
    enddo
enddo
write(*,*) 'ncheck = ',ncheck, ' NBcells: ',NBcells
end subroutine

!----------------------------------------------------------------------------------------
! Identify S1P, CCL21 and OXY influx sites.
! Identify sites where cells can enter and exit.
! Candidate B cell entry sites are at the lower surface of the ellipsoidal blob,
! in the patch delineated by the plane y = -ENTRY_ALPHA*Rb, where:
!    Rb = Ra/ELLIPSE_RATIO
! and Ra and Rb are the major and minor axis radii of the defining ellipse.
!----------------------------------------------------------------------------------------
subroutine AssignBdryRole(site,bdry)
integer :: site(3), r(3)
real :: dy2
type (boundary_type), pointer :: bdry

site = bdry%site
if (site(2) > Centre(2)) then
    bdry%OXY = .true.
    bdry%S1P = .true.
else
    bdry%CCL21 = .true.
    if (Centre(2) - site(2) < 0.3*bRadius) then
        bdry%S1P = .true.
    endif
endif
r = site - Centre
if (r(2) > -ENTRY_ALPHA*bRadius) then
    bdry%entry_ok = .false.
else
    bdry%entry_ok = .true.
endif
if (r(2) > -EXIT_ALPHA*bRadius) then
    bdry%exit_ok = .false.
else
    bdry%exit_ok = .true.
endif 
bdry%CXCL13 = .false.
end subroutine


!----------------------------------------------------------------------------------------
! A new site is to be added to the blob.  The site is chosen so as to maintain the
! ellipsoidal shape.  For each current bdry site, P, imagine a line drawn from the blob centre,
! O, through the site, and determine the point Q where this line intersects the surface of the
! ellipsoid with a = aRadius, b = bRadius.  The distance from the site to the intersection
! point is calculated, and the bdry site that maximizes this distance (the site most inside the
! ellipsoid is chosen as the starting point of the search for a site to add.
! If r = site - Centre, we can parametrise the line by t.r, where t is a scalar:
! x = t.r(1), y = t.r(2), z = t.r(3)
! and t = 1 gives the point P.
! On the ellipsoid with radii a and b:
! x^2/a^2 + (y^2 + z^2)/b^2 = 1
! which implies that t^2[r(1)^2/a^2 + (r(2)^2 + r(3)^2)/b^2] = 1
! from which we obtain t^2.
! Next the 'outside' neighbours of the chosen bdry site are examined, and the one that
! is closest to the centre is selected as the new 'inside' blob site.
! The square of the length of OP is dp^2 = r(1)^2 + r(2)^2 + r(3)^2, and the square of the length
! of OQ is dq^2 = t^2.dp^2, dq = t.dp, and the distance PQ = (t-1).dp.  We want to choose the
! bdry site for which PQ is a maximum.
! Finally the bdry sites in the vicinity are tested to see if any need to be removed from
! the bdrylist. (This is why it makes sense to use a linked list).
! Removal of bdry sites is done by FixBdrylist.
!----------------------------------------------------------------------------------------
subroutine AddSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), indx(2), k, kmin
real :: r(3), x2, y2, z2, dp, t, tmax, dpq, dmax, dmin

!write(*,*) 'AddSite'
!call bdrylist_print(bdrylist)
dmax = -1.0e10
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    t = sqrt(1/(x2/aRadius**2 + (y2+z2)/bRadius**2))
    dpq = (t-1)*dp
    if (dpq > dmax) then
        psite = site
        dmax = dpq
        tmax = t
    endif
    bdry => bdry%next
enddo
!write(*,'(a,3i4,2f8.4)') 'psite, tmax, dmax:     ',psite,tmax,dmax
if (.not.isbdry(psite(1),psite(2),psite(3)))then
    write(*,*) 'psite is not a bdry site'
    stop
endif

! Now we need to look at all 26 neighbours of psite, to determine which 'outside' neighbour
! is closest to O.
dmin = 1.0e10
kmin = 0
do k = 1,27
	if (k == 14) cycle
	site = psite + jumpvec(:,k)
	indx = occupancy(site(1),site(2),site(3))%indx
!	write(*,'(6i8)') k,site,indx
    if (indx(1) >= 0) cycle
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    if (dp < dmin) then
        dmin = dp
        kmin = k
    endif
enddo    
if (kmin == 0) then
    write(*,*) 'Error: no outside neighbours of bdry site'
    ok = .false.
    return
endif

! This is the site to convert to 'inside'
site = psite + jumpvec(:,kmin)
occupancy(site(1),site(2),site(3))%indx = 0
Nsites = Nsites + 1
aRadius = (ELLIPSE_RATIO**2*Nsites*3/(4*PI))**0.33333
bRadius = aRadius/ELLIPSE_RATIO
if (isbdry(site(1),site(2),site(3))) then   ! add it to the bdrylist
    nbdry = nbdry + 1
    allocate(bdry)
    bdry%site = site
    bdry%S1P = .false.
    bdry%CXCL13 = .false.
    bdry%CCL21 = .false.
    bdry%OXY = .false.
!    nullify(bdry%previous)
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    call AssignBdryRole(site,bdry)
!BDRY    occupancy(site(1),site(2),site(3))%bdry = .true.
    occupancy(site(1),site(2),site(3))%bdry => bdry
    call SetBdryConcs(site)
else
    write(*,*) 'Added site is not bdry: ',site
    stop
endif
ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
! Check all sites in bdrylist to ensure that they are still on the bdry, and remove any 
! that are not from bdrylist.
!----------------------------------------------------------------------------------------
subroutine FixBdrylist
integer :: site(3), sitelist(1000,3), k, n
type (boundary_type), pointer :: bdry

!write(*,*) 'FixBdrylist'
n = 0
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    if (.not.isbdry(site(1),site(2),site(3))) then
!        write(*,*) 'Not a bdry site: ',site
        n = n+1
        sitelist(n,:) = site
    endif
    bdry => bdry%next
enddo
do k = 1,n
    site = sitelist(k,:)
    call bdrylist_delete(site,bdrylist)
!BDRY    occupancy(site(1),site(2),site(3))%bdry = .false.
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    nbdry = nbdry - 1
!    write(*,*) 'Removed site from bdrylist: ',site
!    if (site(1)==28 .and. site(2)==48 .and. site(3)==51) then
!        write(*,*) '-----------------------------------------------------------'
!    endif
!    write(*,*) 'Head site: ',bdrylist%site
    if (bdrylist_present(site, bdrylist)) then
        write(*,*) 'Error: still in bdrylist: ',site
        stop
    endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! The logic is similar to AddSite.  We look for the bdry site P which is most outside
! the ellipsoid, i.e. to minimise OQ-OP, where Q is the point at
! which the line through OP intersects the ellipsoid surface.  We then choose the 'inside'
! neighbour of P that is furthest from O.
!----------------------------------------------------------------------------------------
subroutine RemoveSite(ok)
logical :: ok
type (boundary_type), pointer :: bdry
integer :: site(3), psite(3), k, kmax
real :: r(3), x2, y2, z2, dp, t, dpq, dmax, dmin

dmin = 1.0e10
psite = 0
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    t = sqrt(1/(x2/aRadius**2 + (y2+z2)/bRadius**2))
    dpq = (t-1)*dp
    if (dpq < dmin) then
        psite = site
        dmin = dpq
    endif
    bdry => bdry%next
enddo
if (psite(1) == 0) then
    write(*,*) 'No more bdry sites: ',nbdry
    ok = .false.
    return
endif
!write(*,*) 'psite, dmax: ',psite,dmin
if (.not.isbdry(psite(1),psite(2),psite(3))) then
    write(*,*) 'Not a bdry site: ',psite
    ok = .false.
    return
endif

! Now we need to look at psite and its 26 neighbours, to determine which 'inside' neighbour
! is most distant from O.
dmax = -1.0e10
kmax = 0
do k = 1,27
	site = psite + jumpvec(:,k)
    if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle   ! outside
    r = site - Centre
    x2 = r(1)**2
    y2 = r(2)**2
    z2 = r(3)**2
    dp = sqrt(x2 + y2 + z2)
    if (dp > dmax) then
        dmax = dp
        kmax = k
    endif
enddo    
if (kmax == 0) then
    write(*,*) 'Error: no inside neighbours of bdry site'
    ok = .false.
    return
endif
! This is the site to convert to 'outside'
site = psite + jumpvec(:,kmax)
! Note: any occupying cells need to be moved!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call ClearSite(site)
!call GiveChemokines(site)
occupancy(site(1),site(2),site(3))%indx = OUTSIDE_TAG
Nsites = Nsites - 1
aRadius = (ELLIPSE_RATIO**2*Nsites*3/(4*PI))**0.33333
bRadius = aRadius/ELLIPSE_RATIO
!BDRY if (occupancy(site(1),site(2),site(3))%bdry) then   ! remove it from the bdrylist
if (associated(occupancy(site(1),site(2),site(3))%bdry)) then   ! remove it from the bdrylist
    call bdrylist_delete(site, bdrylist)
!BDRY    occupancy(site(1),site(2),site(3))%bdry = .false.
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    nbdry = nbdry - 1
! Need to check for a new bdry site to replace the one removed
! Look at all the neighbours
    psite = site
    do k = 1,27
	    site = psite + jumpvec(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle   ! outside
!BDRY        if (occupancy(site(1),site(2),site(3))%bdry) cycle   ! bdry
        if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle   ! bdry
        if (isbdry(site(1),site(2),site(3))) then
            nbdry = nbdry + 1
            allocate(bdry)
            bdry%site = site
            bdry%S1P = .false.
            bdry%CXCL13 = .false.
            bdry%CCL21 = .false.
            bdry%OXY = .false.
            nullify(bdry%next)
            call bdrylist_insert(bdry,bdrylist)
            call AssignBdryRole(site,bdry)
!BDRY            occupancy(site(1),site(2),site(3))%bdry = .true.
            occupancy(site(1),site(2),site(3))%bdry => bdry
            call SetBdryConcs(site)
        endif
    enddo
endif
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddSites(n,ok)
integer :: n
logical :: ok
integer :: k

do k = 1,n
    call AddSite(ok)
    if (.not.ok) return
    call FixBdrylist
enddo
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine RemoveSites(n,ok)
integer :: n
logical :: ok
integer :: k

do k = 1,n
    call RemoveSite(ok)
    if (.not.ok) return
    call FixBdrylist
enddo
ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine TestAddRemoveSites
integer :: i, n, k, nb
logical :: ok
type (boundary_type), pointer :: bdry

write(*,*) 'TestAddRemoveSites'
n = 10000
do i = 1,10
    write(*,*) 'Adding sites'
    call AddSites(n,ok)
    write(*,*) 'nbdry = ',nbdry, ' NBcells: ',NBcells
    write(*,*) 'Removing sites'
    call RemoveSites(n,ok)
    write(*,*) 'nbdry = ',nbdry, ' NBcells: ',NBcells
enddo
call CreateCheckList
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetBdryConcs(site)
integer :: site(3)

if (use_S1P) then
    S1P_conc(site(1),site(2),site(3)) = BdryS1PConc
endif
if (use_OXY) then
    OXY_conc(site(1),site(2),site(3)) = BdryOXYConc
endif
if (use_CCL21) then
    CCL21_conc(site(1),site(2),site(3)) = BdryCCL21Conc
endif
if (use_CXCL13) then
    CXCL13_conc(site(1),site(2),site(3)) = BdryCXCL13Conc
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Move any cells from site to allow the site to be released (made 'outside')
!-----------------------------------------------------------------------------------------
subroutine ClearSite(csite)
integer :: csite(3)
integer :: i, k, cindx(2), kcell, site(3), indx(2), r
logical :: done

cindx = occupancy(csite(1),csite(2),csite(3))%indx
do i = 2,1,-1
    kcell = cindx(i)
    if (kcell == 0) cycle
    r = 0
    done = .false.
    do while (.not.done)
        r = r + 1
        do k = 1,27
            if (k == 14) cycle
	        site = csite + r*jumpvec(:,k)
	        if (outside_xyz(site(1),site(2),site(3))) cycle
	        indx = occupancy(site(1),site(2),site(3))%indx
            if (indx(1) < 0) cycle  ! outside
            if (indx(1) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(1) = kcell
                cellist(kcell)%site = site
                done = .true.
                exit
            elseif (indx(2) == 0) then  ! use this slot
                occupancy(site(1),site(2),site(3))%indx(2) = kcell
                cellist(kcell)%site = site
                occupancy(csite(1),csite(2),csite(3))%indx(i) = 0
                done = .true.
                exit
            endif
        enddo
    enddo
enddo           
end subroutine

!-----------------------------------------------------------------------------------------
! Using bdrylist, randomly select bdry sites, check that the site is an entrysite and
! suitable for cell ingress.  In fact we look at all neighbouring 'inside' sites.
!-----------------------------------------------------------------------------------------
subroutine getEntrySite(site,ok)
integer :: site(3)
logical :: ok
integer :: ibdry, k, it, bsite(3), indx(2), nt = 100, kpar=0
type (boundary_type), pointer :: bdry

it = 0
do
    it = it+1
    if (it == nt) then
        ok = .false.
        return
    endif
    ibdry = random_int(1,nbdry,kpar)
    k = 0
    bdry => bdrylist
    do while ( associated ( bdry )) 
        k = k+1
        if (k == ibdry) exit
        bdry => bdry%next
    enddo
    if (.not.bdry%entry_ok) cycle
    do k = 1,27
	    site = bdry%site + jumpvec(:,k)
	    indx = occupancy(site(1),site(2),site(3))%indx
        if (indx(1) < 0) cycle                      ! outside
        if (indx(1) == 0 .or. indx(2) == 0) then    ! free slot
            ok = .true.
            return
        endif
    enddo     
enddo
end subroutine


!----------------------------------------------------------------------------------------
! When a site is removed from the blob (freed), the mass of each chemokine in that gridcell
! must be shared out among the neighbour sites.
!----------------------------------------------------------------------------------------
subroutine GiveChemokines(site)
integer :: site(3)

end subroutine

!----------------------------------------------------------------------------------------
! We can solve for a steady-state concentration field, compute the gradient field from it,
! and use this for all chemotaxis calculations until the ellipsoid size changes.  At this
! point we need to resolve for the concentration field.
! There are two possible approaches.  Either to solve for the steady-state field again,
! using the previous solution as the starting condition, or to allow the field to evolve.
! In the latter case the concentration in some number of sites in the neighbourhood of the
! added/removed site is averaged in a way that ensures mass conservation.
! This all works, but there is now serious doubt as to whether there is such a thing as
! S1P chemotaxis.  According to Irina Grigorova, S1P1 on a T cell enables it to cross
! the endothelial boundary into a sinus (high-S1P), but the cell must reach the sinus
! by random motion.
!----------------------------------------------------------------------------------------
subroutine InitFields
integer :: i

if (use_CXCL13) then
	call init_CXCL13
endif
if (use_S1P) then
	call init_S1P
endif
if (use_OXY) then
	call init_OXY
endif
if (use_CCL21) then
	call init_CCL21
endif
call ShowConcs

!do i = 1,100
!    call evolveS1P(6,0.0)
!enddo
!call ShowConcs
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowConcs
integer :: x, y, z

x = NX/2
z = NZ/2
do y = 1,NY
    if (occupancy(x,y,z)%indx(1) < 0) cycle
    write(*,'(i4,4e12.4)') y,S1P_conc(x,y,z),CCL21_conc(x,y,z) !,OXY_conc(x,y,z)
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Set up S1P_influx
! The influx array is used both to flag sites that are within the blob and also
! to convey chemokine influx rates at the boundary.
!   influx = -1  outside blob
!             0  inside blob, no chemokine influx
!            >0  inside blob, chemokine influx (i.e. a bdry gridcell)
!----------------------------------------------------------------------------------------
!subroutine SetupS1P_influx
!integer :: x, y, z, indx(2), site(3), nc
!type (boundary_type), pointer :: bdry

!S1P_influx = -1
!do x = 1,NX
!    do y = 1,NY
!        do z = 1,NZ
!            indx = occupancy(x,y,z)%indx
!            if (indx(1) < 0) cycle
!            S1P_influx(x,y,z) = 0
!        enddo
!    enddo
!enddo
!    bdry => bdrylist
!    do while ( associated ( bdry )) 
!        site = bdry%site
!        nc = 0
!        if (bdry%S1P) then
!            do k = 1,6
!                x = site(1) + neumann(1,k)
!                y = site(2) + neumann(2,k)
!                z = site(3) + neumann(3,k)
!                if (outside_xyz(x,y,z)) cycle
!                if (occupancy(x,y,z)%indx(1) < 0) then
!                    nc = nc + 1
!                endif
!            enddo
!            S1P_influx(site(1),site(2),site(3)) = nc*S1P_rate
!        bdry => bdry%next
!    enddo

!end subroutine

!----------------------------------------------------------------------------------------
! The flux of S1P into the gridcell (x,y,z) is proportional to influx(x,y,z)
! The first task is to determine the equilibrium S1P concentration field.
! One way is to solve the diffusion equation over a time long enough to reach steady-state.
! S1P both diffuses and decays.
! Treat the grid spacing as the unit of distance, therefore gridcell volume = 1.
! Note that with gridcell volume = 1 the mass of S1P = concentration.
! Note that influx(x,y,z) is used to flag the region of diffusion of S1P,
! by setting influx(x,y,z) = -1 for gridcells that are not within the follicle or not
! available for chemokine diffusion for some other reason (e.g. FDC).
! influx(x,y,z) >= 0 for gridcells within the region.
!
! A problem:
! The set of gridcells (sites) that experience chemokine influx varies as the follicle
! expands and contracts (for boundary influx) and as chemokine-secreting cells move.
! For now we consider only boundary influx.
! We need to maintain a list of boundary sites receiving chemokine influx, for each
! constituent.  The easiest method is to maintain a list of all boundary sites,
! with associated info indicating chemokine influx and other data, if required.
!----------------------------------------------------------------------------------------
subroutine init_S1P
integer :: x, y, z, xx, yy, zz, k, nb
real :: g(3), gamp, gmax
logical :: steady = .true.

write(logmsg,*) 'Initializing S1P'
call logger(logmsg)
allocate(S1P_conc(NX,NY,NZ))
allocate(S1P_grad(3,NX,NY,NZ))
allocate(S1P_influx(NX,NY,NZ))

!call SetupS1P_influx

call SolveSteadystate(S1P,S1P_KDIFFUSION,S1P_KDECAY,S1P_conc)
! Now compute the gradient field.
call gradient(S1P_conc,S1P_grad)
!write(*,*) 'S1P gradient: ',S1P_grad(:,25,25,25)
gmax = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			g = S1P_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,*) 'Max S1P gradient: ',gmax
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine init_OXY
integer :: x, y, z, xx, yy, zz, k, nb
real :: g(3), gamp, gmax
logical :: steady = .true.

write(logmsg,*) 'Initializing OXY'
call logger(logmsg)
allocate(OXY_conc(NX,NY,NZ))
allocate(OXY_grad(3,NX,NY,NZ))
allocate(OXY_influx(NX,NY,NZ))

OXY_influx = -1

! Set up OXY_influx

call SolveSteadystate(OXY,OXY_KDIFFUSION,OXY_KDECAY,OXY_conc)
! Now compute the gradient field.
call gradient(OXY_conc,OXY_grad)
gmax = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			g = OXY_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,*) 'Max OXY gradient: ',gmax
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine init_CCL21
integer :: x, y, z, xx, yy, zz, k, nb
real :: g(3), gamp, gmax
logical :: steady = .true.

write(logmsg,*) 'Initializing CCL21'
call logger(logmsg)
allocate(CCL21_conc(NX,NY,NZ))
allocate(CCL21_grad(3,NX,NY,NZ))
allocate(CCL21_influx(NX,NY,NZ))

CCL21_influx = -1

! Set up CCL21_influx

call SolveSteadystate(CCL21,CCL21_KDIFFUSION,CCL21_KDECAY,CCL21_conc)
! Now compute the gradient field.
call gradient(CCL21_conc,CCL21_grad)
gmax = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			g = CCL21_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,*) 'Max CCL21 gradient: ',gmax
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine init_CXCL13
integer :: isignal, site(3), x, y, z, it, nt, isig, iblast
real :: dt, sum
real :: g(3), gamp, gmax
!real, parameter :: Kdecay = 0.00001, Kdiffusion = 0.001

write(logmsg,*) 'Initializing CXCL13: KDECAY: ', CXCL13_KDECAY
!call formulate(CXCL12_KDIFFUSION,CXCL12_KDECAY)

call logger(logmsg)
allocate(CXCL13_conc(NX,NY,NZ))
allocate(CXCL13_grad(3,NX,NY,NZ))
allocate(CXCL13_influx(NX,NY,NZ))

CXCL13_grad = 0
CXCL13_influx = -1

! Set up CXCL13_influx

call SolveSteadystate(CXCL13,CXCL13_KDIFFUSION,CXCL13_KDECAY,CXCL13_conc)
call gradient(CXCL13_conc,CXCL13_grad)
! Testing convergence
!nt = 10
!dt = 100*DELTA_T
!do it = 1,nt
!	call evolve(CXCL13_influx,CXCL13_KDIFFUSION,CXCL13_KDECAY,CXCL13_conc,dt)
!	sum = 0
!	do isig = 1,nsignal
!		site = signal(isig)%site
!		sum = sum + CXCL13_conc(site(1),NBY+1,site(3))
!	enddo
!	write(logmsg,*) 'Mean patch CXCL13: ',sum/nsignal
!	call logger(logmsg)
!enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine EvolveS1P(nt,totsig)
integer :: nt
real :: totsig
integer :: x, y, z, site(3), isig, i, iblast
real :: dt, sum, sig

dt = nt*DELTA_T
call evolve(S1P,S1P_KDIFFUSION,S1P_KDECAY,S1P_conc,dt)
call gradient(S1P_conc,S1P_grad)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine EvolveOXY(nt,totsig)
integer :: nt
real :: totsig
integer :: x, y, z, site(3), isig, i, iblast
real :: dt, sum, sig

dt = nt*DELTA_T
call evolve(OXY,OXY_KDIFFUSION,OXY_KDECAY,OXY_conc,dt)
call gradient(OXY_conc,OXY_grad)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine EvolveCCL21(nt,totsig)
integer :: nt
real :: totsig
integer :: x, y, z, site(3), isig, i, iblast
real :: dt, sum, sig

dt = nt*DELTA_T
call evolve(CCL21,CCL21_KDIFFUSION,CCL21_KDECAY,CCL21_conc,dt)
call gradient(CCL21_conc,CCL21_grad)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine EvolveCXCL13(nt,totsig)
integer :: nt
real :: totsig
integer :: x, y, z, site(3), isig, i, iblast
real :: dt, sum, sig

dt = nt*DELTA_T
call evolve(CXCL13,CXCL13_KDIFFUSION,CXCL13_KDECAY,CXCL13_conc,dt)
call gradient(CXCL13_conc,CXCL13_grad)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine SolveSteadystate(ichemo,Kdiffusion,Kdecay,C)
integer :: ichemo
real :: C(:,:,:)
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, maxchange_par, total_par
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-6		! max change in C at any site as fraction of average C
integer :: nc, nc_par, k, it, n, kpar
integer :: xlim(2,16), dx, x1, x2, xfr, xto	! max 16 threads
real, allocatable :: C_par(:,:,:)
integer :: nt = 10000
real, allocatable :: Ctemp(:,:,:)

dx = NX/Mnodes + 0.5
do k = 1,Mnodes
	xlim(1,k) = (k-1)*dx + 1
	xlim(2,k) = k*dx
enddo
xlim(2,Mnodes) = NX

C = 0
do it = 1,nt
    if (mod(it,10) == 0) then
        write(*,*) 'iteration: ',it
    endif
	maxchange = 0
	total = 0
	nc = 0
	!$omp parallel do private(x1,x2,n,C_par,maxchange_par,total_par,nc_par)
	do kpar = 0,Mnodes-1
		x1 = xlim(1,kpar+1)
		x2 = xlim(2,kpar+1)
		n = x2 - x1 + 1
		allocate(C_par(n,NY,NZ))
		call par_steadystate(ichemo,C,Kdiffusion,Kdecay,C_par,x1,x2,maxchange_par,total_par,nc_par,kpar)
		C(x1:x2,:,:) = C_par(1:n,:,:)
		deallocate(C_par)
		nc = nc + nc_par
		total = total + total_par
		maxchange = max(maxchange,maxchange_par)
	enddo
	if (maxchange < tol*total/nc) then
		write(logmsg,*) 'Convergence reached: it: ',it
		call logger(logmsg)
		exit
	endif
enddo
end subroutine	

!----------------------------------------------------------------------------------------
! Note that the array Ctemp() is defined with z:z1-z2, i.e. with
! a range of n = z2 - z1 + 1 values.  The read accesses of C() and influx() are shared 
! in this version.  Does this have a big penalty?
!----------------------------------------------------------------------------------------
subroutine par_steadystate1(C,Kdiffusion,Kdecay,Ctemp,influx,z1,z2,maxchange,total,nc,kpar)
real :: C(:,:,:), Ctemp(:,:,:), influx(:,:,:)
integer :: z1, z2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, zpar

dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3

maxchange = 0
total = 0
nc = 0
do zpar = 1,z2-z1+1
	z = zpar + z1-1
	do y = 1,NY
		do x = 1,NX
			if (influx(x,y,z) < 0) cycle
			if (influx(x,y,z) > 0) then
				dC = influx(x,y,z)*dx2diff
			else
				dC = 0
			endif
			sum = 0
			nb = 0
			do k = 1,6
				xx = x + neumann(1,k)
				yy = y + neumann(2,k)
				zz = z + neumann(3,k)
				if (outside_xyz(xx,yy,zz)) cycle
				if (influx(xx,yy,zz) < 0) cycle
				nb = nb + 1
				sum = sum + C(xx,yy,zz)
			enddo
!			Ctemp(x,y,z) = alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) + (1-alpha)*C(x,y,z)
			Ctemp(x,y,zpar) = alpha*(DELTA_X*Kdiffusion*sum + influx(x,y,z))/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
			dC = abs(Ctemp(x,y,zpar) - C(x,y,z))
			maxchange = max(dC,maxchange)
			nc = nc + 1
			total = total + Ctemp(x,y,zpar)
		enddo
	enddo
enddo

end subroutine

!----------------------------------------------------------------------------------------
! A different approach from that used in bone-abm.
! The bdry gridcells that are adjacent to a chemokine-rich region (i.e. there are 
! chemokine-secreting cells near the boundary) are given a fixed concentration.
! All boundaries are no-flux.
! This is a simplified approach.  We need only keep track of the boundary.
! Note: It would speed things up significantly if instead of simply flagging a
! site as %bdry = .true., we make %bdry a pointer into the bdrylist.  How hard is this?
! Seems to be OK.
! The concentrations in a specified-conc bdry site are input parameters.
!----------------------------------------------------------------------------------------
subroutine par_steadystate_z(ichemo,C,Kdiffusion,Kdecay,Ctemp,z1,z2,maxchange,total,nc,kpar)
integer :: ichemo
real :: C(:,:,:), Ctemp(:,:,:)
integer :: z1, z2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, zpar, indx(2)
logical :: bdry_conc

!write(*,*) 'par_steadystate: ',ichemo,z1,z2
dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3

maxchange = 0
total = 0
nc = 0
do zpar = 1,z2-z1+1
	z = zpar + z1-1
	do y = 1,NY
		do x = 1,NX
		    indx = occupancy(x,y,z)%indx
!			if (influx(x,y,z) < 0) cycle
            if (indx(1) < 0) cycle      ! outside
!			if (influx(x,y,z) > 0) then
!				dC = influx(x,y,z)*dx2diff
!			else
!				dC = 0
!			endif
            bdry_conc = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for S1P bdry site - no change to the concentration at such a site
                if (ichemo == S1P .and. occupancy(x,y,z)%bdry%S1P) then
                    C(x,y,zpar) = BdryS1PConc
                    Ctemp(x,y,zpar) = C(x,y,zpar)
                    bdry_conc = .true.
                endif
                if (ichemo == OXY .and. occupancy(x,y,z)%bdry%OXY) then
                    C(x,y,zpar) = BdryOXYConc
                    Ctemp(x,y,zpar) = C(x,y,zpar)
                    bdry_conc = .true.
                endif
                if (ichemo == CCL21 .and. occupancy(x,y,z)%bdry%CCL21) then
                    C(x,y,zpar) = BdryCCL21Conc
                    Ctemp(x,y,zpar) = C(x,y,zpar)
                    bdry_conc = .true.
                endif
                if (ichemo == CXCL13 .and. occupancy(x,y,z)%bdry%CXCL13) then
                    C(x,y,zpar) = BdryCXCL13Conc
                    Ctemp(x,y,zpar) = C(x,y,zpar)
                    bdry_conc = .true.
                endif
            endif
            if (.not.bdry_conc) then
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
!    				if (influx(xx,yy,zz) < 0) cycle
    				if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
!    			Ctemp(x,y,z) = alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) + (1-alpha)*C(x,y,z)
!			    Ctemp(x,y,zpar) = alpha*(DELTA_X*Kdiffusion*sum + influx(x,y,z))/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
			    Ctemp(x,y,zpar) = alpha*(DELTA_X*Kdiffusion*sum)/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
!    if (Ctemp(x,y,zpar) > 0) then
!        write(*,'(3i4,f8.4)') x,y,zpar,Ctemp(x,y,zpar)
!    endif
			    dC = abs(Ctemp(x,y,zpar) - C(x,y,z))
			    maxchange = max(dC,maxchange)
			endif
			nc = nc + 1
			total = total + Ctemp(x,y,zpar)
		enddo
	enddo
enddo
!write(*,*) 'maxchange, average: ',maxchange,total/nc

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine par_steadystate(ichemo,C,Kdiffusion,Kdecay,Ctemp,x1,x2,maxchange,total,nc,kpar)
integer :: ichemo
real :: C(:,:,:), Ctemp(:,:,:)
integer :: x1, x2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, xpar, indx(2)
logical :: bdry_conc

!write(*,*) 'par_steadystate: ',ichemo,z1,z2
dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3

maxchange = 0
total = 0
nc = 0
do xpar = 1,x2-x1+1
	x = xpar + x1-1
	do y = 1,NY
		do z = 1,NZ
		    indx = occupancy(x,y,z)%indx
            if (indx(1) < 0) cycle      ! outside
            bdry_conc = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for chemo bdry site - no change to the concentration at such a site
                if (ichemo == S1P .and. occupancy(x,y,z)%bdry%S1P) then
                    C(x,y,z) = BdryS1PConc
                    Ctemp(xpar,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == OXY .and. occupancy(x,y,z)%bdry%OXY) then
                    C(x,y,z) = BdryOXYConc
                    Ctemp(xpar,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == CCL21 .and. occupancy(x,y,z)%bdry%CCL21) then
                    C(x,y,z) = BdryCCL21Conc
                    Ctemp(xpar,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == CXCL13 .and. occupancy(x,y,z)%bdry%CXCL13) then
                    C(x,y,z) = BdryCXCL13Conc
                    Ctemp(xpar,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
            endif
            if (.not.bdry_conc) then
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
    				if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
			    Ctemp(xpar,y,z) = alpha*(DELTA_X*Kdiffusion*sum)/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
			    dC = abs(Ctemp(xpar,y,z) - C(x,y,z))
			    maxchange = max(dC,maxchange)
			endif
			nc = nc + 1
			total = total + Ctemp(xpar,y,z)
		enddo
	enddo
enddo
!write(*,*) 'maxchange, average: ',maxchange,total/nc

end subroutine
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine gradient(C,grad)
!real :: influx(:,:,:)
real :: C(:,:,:), grad(:,:,:,:)
integer :: x, y, z, xx, yy, zz, x1, x2, y1, y2, z1, z2, i, k
real :: g(3)
logical :: missed
real, parameter :: MISSING_VAL = 1.0e10

grad = 0
do z = 1,NZ
	do y = 1,NY
		do x = 1,NX
!			if (influx(x,y,z) < 0) cycle
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			x1 = x - 1
			x2 = x + 1
			if (x1 < 1 .or. x2 > NX) then
				g(1) = 0
!			elseif (influx(x1,y,z) >= 0 .and. influx(x2,y,z) >= 0) then
			elseif (occupancy(x1,y,z)%indx(1) >= 0 .and. occupancy(x2,y,z)%indx(1) >= 0) then
				g(1) = (C(x2,y,z) - C(x1,y,z))/(2*DELTA_X)
			else
				g(1) = MISSING_VAL
			endif
			y1 = y - 1
			y2 = y + 1
			if (y1 < 1 .or. y2 > NY) then
				g(2) = 0
!			elseif (influx(x,y1,z) >= 0 .and. influx(x,y2,z) >= 0) then
			elseif (occupancy(x,y1,z)%indx(1) >= 0 .and. occupancy(x,y2,z)%indx(1) >= 0) then
				g(2) = (C(x,y2,z) - C(x,y1,z))/(2*DELTA_X)
			else
				g(2) = MISSING_VAL
			endif
			z1 = z - 1
			z2 = z + 1
			if (z1 < 1 .or. z2 > NZ) then
				g(3) = 0
!			elseif (influx(x,y,z1) >= 0 .and. influx(x,y,z2) >= 0) then
			elseif (occupancy(x,y,z1)%indx(1) >= 0 .and. occupancy(x,y,z2)%indx(1) >= 0) then
				g(3) = (C(x,y,z2) - C(x,y,z1))/(2*DELTA_X)
			else
				g(3) = MISSING_VAL
			endif
			grad(:,x,y,z) = g
		enddo
	enddo
enddo
do z = 1,NZ
	do y = 1,NY
		do x = 1,NX
!			if (influx(x,y,z) < 0) cycle
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			do i = 1,3
				if (grad(i,x,y,z) == MISSING_VAL) then
					missed = .true.
					grad(i,x,y,z) = 0
					do k = 1,6
						xx = x + neumann(1,k)
						yy = y + neumann(2,k)
						zz = z + neumann(3,k)
						if (outside_xyz(xx,yy,zz)) cycle
!						if (influx(xx,yy,zz) < 0) cycle
						if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
						if (grad(i,xx,yy,zz) /= MISSING_VAL) then
							grad(i,x,y,z) = grad(i,xx,yy,zz)
							missed = .false.
							exit
						endif
					enddo
					if (missed) then
!						write(*,*) 'Missing gradient at: ',x,y,z,grad(:,x,y,z)
					endif
				endif
			enddo
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine evolve(ichemo,Kdiffusion,Kdecay,C,dt)
integer :: ichemo
real :: C(:,:,:)
real :: Kdiffusion, Kdecay, dt
real :: dx2diff, total, maxchange, C0, dC, sum, dV, dMdt
real, parameter :: alpha = 0.99
integer :: x, y, z, xx, yy, zz, nb, nc, k, it
real, allocatable :: Ctemp(:,:,:)
logical :: bdry_conc

dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
do z = 1,NZ
	do y = 1,NY
		do x = 1,NX
!			if (influx(x,y,z) < 0) cycle
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			C0 = C(x,y,z)
            bdry_conc = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for S1P bdry site - no change to the concentration at such a site
                if (ichemo == S1P .and. occupancy(x,y,z)%bdry%S1P) then
                    C(x,y,z) = BdryS1PConc
                    Ctemp(x,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == OXY .and. occupancy(x,y,z)%bdry%OXY) then
                    C(x,y,z) = BdryOXYConc
                    Ctemp(x,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == CCL21 .and. occupancy(x,y,z)%bdry%CCL21) then
                    C(x,y,z) = BdryCCL21Conc
                    Ctemp(x,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
                if (ichemo == CXCL13 .and. occupancy(x,y,z)%bdry%CXCL13) then
                    C(x,y,z) = BdryCXCL13Conc
                    Ctemp(x,y,z) = C(x,y,z)
                    bdry_conc = .true.
                endif
            endif
            if (.not.bdry_conc) then
			    sum = 0
			    nb = 0
			    do k = 1,6
				    xx = x + neumann(1,k)
				    yy = y + neumann(2,k)
				    zz = z + neumann(3,k)
				    if (outside_xyz(xx,yy,zz)) cycle
!    				if (influx(xx,yy,zz) < 0) cycle
				    if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
				    nb = nb + 1
				    sum = sum + C(xx,yy,zz)
			    enddo
			    dMdt = Kdiffusion*DELTA_X*(sum - nb*C0) - Kdecay*C0*dV ! + influx(x,y,z)
			    Ctemp(x,y,z) = (C0*dV + dMdt*dt)/dV
			endif
		enddo
	enddo
enddo
C = Ctemp
deallocate(Ctemp)

end subroutine


end module