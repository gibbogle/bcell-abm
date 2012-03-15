! Chemokine concentration fields

module fields
use global
use bdry_linked_list
implicit none
save

type chemo_type
	character(8) :: name
	logical :: used
	real :: bdry_conc
	real :: diff_coeff
	real :: halflife
	real :: decay_rate
	real :: strength
	real, allocatable :: conc(:,:,:)
	real, allocatable :: grad(:,:,:,:)
end type

integer, parameter :: MAX_BDRY = 20000

type(chemo_type) :: chemo(MAX_CHEMO)
logical :: use_S1P = .true.
logical :: use_CXCL13 = .false.
logical :: use_CCL21 = .true.
logical :: use_OXY = .false.

! The following are no longer used
real, parameter :: S1P_KDIFFUSION = 0.001
real, parameter :: S1P_KDECAY = 0.000005
real, parameter :: OXY_KDIFFUSION = 0.001
real, parameter :: OXY_KDECAY = 0.000005
real, parameter :: CCL21_KDIFFUSION = 0.001
real, parameter :: CCL21_KDECAY = 0.00001
real, parameter :: CXCL13_KDIFFUSION = 1.0
real, parameter :: CXCL13_KDECAY = 1.0

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
                bdry%chemo_influx = .false.
!                bdry%S1P = .false.
!                bdry%CXCL13 = .false.
!                bdry%CCL21 = .false.
!                bdry%OXY = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,bdrylist)
                call AssignBdryRole(site,bdry)
                occupancy(x,y,z)%bdry => bdry
            endif
        enddo
    enddo
enddo
!call bdrylist_print(bdrylist)
!call TestAddRemoveSites
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
                bdry%chemo_influx = .false.
!                bdry%S1P = .false.
!                bdry%CXCL13 = .false.
!                bdry%CCL21 = .false.
!                bdry%OXY = .false.
                nullify(bdry%next)
                call bdrylist_insert(bdry,checklist)
                call AssignBdryRole(site,bdry)
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
integer :: site(3), dy
type (boundary_type), pointer :: bdry
real, parameter :: S1Pfraction = 0.25

bdry%chemo_influx = .false.
!bdry%S1P = .false.
!bdry%CCL21 = .false.
!bdry%OXY = .false.
!bdry%CXCL13 = .false.
site = bdry%site
dy = site(2) - Centre(2)
if (dy > 0) then
	bdry%chemo_influx(S1P) = .true.
	bdry%chemo_influx(OXY) = .true.
!    bdry%S1P = .true.
!    bdry%OXY = .true.
else
    if (dy > -S1Pfraction*bRadius) then
		bdry%chemo_influx(S1P) = .true.
!        bdry%S1P = .true.
    else
		bdry%chemo_influx(CCL21) = .true.
!	    bdry%CCL21 = .true.
    endif
endif
if (dy > -ENTRY_ALPHA*bRadius) then
    bdry%entry_ok = .false.
else
    bdry%entry_ok = .true.
endif
if (dy > -EXIT_ALPHA*bRadius) then
    bdry%exit_ok = .false.
else
    bdry%exit_ok = .true.
endif 
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
integer :: site(3), psite(3), indx(2), k, kmin, x, y, z
real :: r(3), x2, y2, z2, dp, t, tmax, dpq, dmax, dmin
real :: cmin, cmax

!write(*,*) 'AddSite'
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
if (.not.isbdry(psite(1),psite(2),psite(3)))then
    write(logmsg,*) 'psite is not a bdry site'
	call logger(logmsg)
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
    write(logmsg,*) 'Error: no outside neighbours of bdry site'
	call logger(logmsg)
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
    bdry%chemo_influx = .false.
!    bdry%S1P = .false.
!    bdry%CXCL13 = .false.
!    bdry%CCL21 = .false.
!    bdry%OXY = .false.
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    call AssignBdryRole(site,bdry)
    occupancy(site(1),site(2),site(3))%bdry => bdry
    call SetBdryConcs(site)
else
    write(logmsg,*) 'Added site is not bdry: ',site
	call logger(logmsg)
    stop
endif
ok = .true.
return

!cmin = 1.0e10
!cmax = 0
!do x = 1,NX
!	do y = 1,NY
!		do z = 1,NZ
!			if (occupancy(x,y,z)%indx(1) >= 0) then
!				cmin = min(cmin,S1P_conc(x,y,z))
!				cmax = max(cmin,S1P_conc(x,y,z))
!			endif
!		enddo
!	enddo
!enddo
!write(*,*) 'S1P_conc: min, max: ',cmin,cmax
!stop
end subroutine

!----------------------------------------------------------------------------------------
! Check all sites in bdrylist to ensure that they are still on the bdry, and remove any 
! that are not from bdrylist.
! When a site is removed from the list its chemokine concs and gradients are replaced by
! neighbourhood averages.
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
        n = n+1
        sitelist(n,:) = site
    endif
    bdry => bdry%next
enddo
do k = 1,n
    site = sitelist(k,:)
    call bdrylist_delete(site,bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    nbdry = nbdry - 1
    if (bdrylist_present(site, bdrylist)) then
        write(logmsg,*) 'Error: FixBdrylist: still in bdrylist: ',site
		call logger(logmsg)
        stop
    endif
    call SetConcs(site)
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
    write(logmsg,*) 'No more bdry sites: ',nbdry
	call logger(logmsg)
    ok = .false.
    return
endif
if (.not.isbdry(psite(1),psite(2),psite(3))) then
    write(logmsg,*) 'Not a bdry site: ',psite
	call logger(logmsg)
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
    write(logmsg,*) 'Error: no inside neighbours of bdry site'
	call logger(logmsg)
    ok = .false.
    return
endif
! This is the site to convert to 'outside'
site = psite + jumpvec(:,kmax)
call ClearSite(site)
occupancy(site(1),site(2),site(3))%indx = OUTSIDE_TAG
Nsites = Nsites - 1
aRadius = (ELLIPSE_RATIO**2*Nsites*3/(4*PI))**0.33333
bRadius = aRadius/ELLIPSE_RATIO
if (associated(occupancy(site(1),site(2),site(3))%bdry)) then   ! remove it from the bdrylist
    call bdrylist_delete(site, bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
    nbdry = nbdry - 1
! Need to check for a new bdry site to replace the one removed
! Look at all the neighbours
    psite = site
    do k = 1,27
	    site = psite + jumpvec(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle   ! outside
        if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle   ! bdry
        if (isbdry(site(1),site(2),site(3))) then
            nbdry = nbdry + 1
            allocate(bdry)
            bdry%site = site
            bdry%chemo_influx = .false.
!            bdry%S1P = .false.
!            bdry%CXCL13 = .false.
!            bdry%CCL21 = .false.
!            bdry%OXY = .false.
            nullify(bdry%next)
            call bdrylist_insert(bdry,bdrylist)
            call AssignBdryRole(site,bdry)
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
! Check consistency between bdrylist and occupancy%bdry
!-----------------------------------------------------------------------------------------
subroutine CheckBdryList
type (boundary_type), pointer :: bdry => null(), ocbdry => null()
integer :: x, y, z, site(3), dy, nb1, nb2, nbx

nb1 = 0
nbx = 0
bdry => bdrylist
do while ( associated ( bdry )) 
	nb1 = nb1 + 1
	if (bdry%exit_ok) nbx = nbx + 1
    site = bdry%site
    ocbdry => occupancy(site(1),site(2),site(3))%bdry
    if (.not.associated(ocbdry,bdry)) then
		write(logmsg,*) 'Error: CheckBdryList: inconsistent bdry pointers at site: ',site
		call logger(logmsg)
		stop
	endif
    bdry => bdry%next
enddo
nb2 = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			bdry => occupancy(x,y,z)%bdry
			if (associated(bdry)) then
				nb2 = nb2 + 1
				dy = y - Centre(2)
!				if (dy > -EXIT_ALPHA*bRadius) then
!					if (bdry%exit_ok) then
!						write(*,*) 'Error: bdry%exit_ok wrongly true: ',x,y,z, dy, -EXIT_ALPHA*bRadius
!					endif
!				else
!					if (.not.bdry%exit_ok) then
!						write(*,*) 'Error: bdry%exit_ok wrongly false: ',x,y,z, dy, -EXIT_ALPHA*bRadius
!					endif
!				endif
			else
				if (isbdry(x,y,z)) then
					write(logmsg,*) 'Error: boundary site not in list: ',x,y,z
					call logger(logmsg)
					stop
				endif
			endif
		enddo
	enddo				
enddo
if (nb1 /= nb2) then
	write(logmsg,*) 'Error: inconsistent boundary site counts: ',nb1,nb2
	call logger(logmsg)
	stop
endif
write(logmsg,*) 'bdrylist and occupancy%bdry are consistent: ',nb1,nbx
call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! The added site needs approximate values for chemokine concentrations and gradients,
! as an interim measure until the next steady-state computation.  Actually only the
! gradient is used (for chemotaxis) but the concentration serves as the starting value
! when the steady-state solution is computed.
!-----------------------------------------------------------------------------------------
subroutine SetBdryConcs(site)
integer :: site(3)
integer :: i
type (boundary_type), pointer :: bdry
real :: cbnd

bdry => occupancy(site(1),site(2),site(3))%bdry
if (.not.associated(bdry)) then
	write(logmsg,*) 'Error: SetBdryConcs: not a bdry site'
	call logger(logmsg)
	stop
endif

do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		if (bdry%chemo_influx(i)) then
			cbnd = chemo(i)%bdry_conc
		else
			cbnd = -1
		endif 
		call AverageBdryConc(bdry,chemo(i)%conc,chemo(i)%grad,site,cbnd)
	endif
enddo

!if (use_S1P) then
!	if (bdry%S1P) then
!		cbnd = BdryS1PConc
!	else
!		cbnd = -1
!	endif 
!	call AverageBdryConc(bdry,S1P_conc,S1P_grad,site,cbnd)
!endif
!if (use_OXY) then
!	if (bdry%OXY) then
!		cbnd = BdryOXYConc
!	else
!		cbnd = -1
!	endif 
!	call AverageBdryConc(bdry,OXY_conc,OXY_grad,site,cbnd)
!endif
!if (use_CCL21) then
!	if (bdry%CCL21) then
!		cbnd = BdryCCL21Conc
!	else
!		cbnd = -1
!	endif 
!	call AverageBdryConc(bdry,CCL21_conc,CCL21_grad,site,cbnd)
!endif
!if (use_CXCL13) then
!	if (bdry%CXCL13) then
!		cbnd = BdryCXCL13Conc
!	else
!		cbnd = -1
!	endif 
!	call AverageBdryConc(bdry,CXCL13_conc,CXCL13_grad,site,cbnd)
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! This subroutine is for the case of a site that was boundary and is now in the blob
! interior.  The chemokine gradients are adjusted to reduce the inaccuracy in chemotactic 
! effects until the new steady-state is computed, while the aim in adjusting the concs
! is to speed up convergence in the steady-state computation.  In fact it doesn't have
! much effect on the convergence.  There might be a better way to do the adjustment than
! the simple neighbourhood averaging that is used.
!-----------------------------------------------------------------------------------------
subroutine SetConcs(site)
integer :: site(3)
integer :: i

do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		call AverageConc(chemo(i)%conc,chemo(i)%grad,site)
	endif
enddo
!if (use_S1P) then
!	call AverageConc(S1P_conc,S1P_grad,site)
!endif
!if (use_OXY) then
!	call AverageConc(OXY_conc,OXY_grad,site)
!endif
!if (use_CCL21) then
!	call AverageConc(CCL21_conc,CCL21_grad,site)
!endif
!if (use_CXCL13) then
!	call AverageConc(CXCL13_conc,CXCL13_grad,site)
!endif
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
subroutine GetEntrySite(site,ok)
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
! If the site is a source boundary for the chemokine, cbnd = boundary concentration,
! otherwise cbnd = -1, which is a flag to compute the estimate of concentration by
! averaging the over neighbour sites.
! For the gradient, averaging is always used.
!----------------------------------------------------------------------------------------
subroutine AverageBdryConc(bdry,C,G,site,cbnd)
integer :: site(3)
real :: C(:,:,:), G(:,:,:,:)
real :: cbnd
integer :: x, y, z, k, nave
real :: cave, gave(3)
type (boundary_type), pointer :: bdry

cave = 0
gave = 0
nave = 0
do k = 1,27
	if (k == 14) cycle
    x = site(1) + jumpvec(1,k)
    y = site(2) + jumpvec(2,k)
    z = site(3) + jumpvec(3,k)
    if (outside_xyz(x,y,z)) cycle
    if (occupancy(x,y,z)%indx(1) < 0) cycle
    nave = nave + 1
    cave = cave + C(x,y,z)
    gave = gave + G(:,x,y,z)
enddo
if (cbnd >= 0) then
	C(site(1),site(2),site(3)) = cbnd
else
	if (nave > 0) then
		C(site(1),site(2),site(3)) = cave/nave
	else
		C(site(1),site(2),site(3)) = 0
	endif
endif
if (nave > 0) then
	G(:,site(1),site(2),site(3)) = gave/nave
else
	G(:,site(1),site(2),site(3)) = 0
endif
end subroutine

!----------------------------------------------------------------------------------------
! The concentration and gradient are estimated by averaging the over neighbour sites.
!----------------------------------------------------------------------------------------
subroutine AverageConc(C,G,site)
integer :: site(3)
real :: C(:,:,:), G(:,:,:,:)
integer :: x, y, z, k, nave
real :: cave, gave(3)

cave = 0
gave = 0
nave = 0
do k = 1,27
	if (k == 14) cycle
    x = site(1) + jumpvec(1,k)
    y = site(2) + jumpvec(2,k)
    z = site(3) + jumpvec(3,k)
    if (outside_xyz(x,y,z)) cycle
    if (occupancy(x,y,z)%indx(1) < 0) cycle
    nave = nave + 1
    cave = cave + C(x,y,z)
    gave = gave + G(:,x,y,z)
enddo
if (nave > 0) then
	C(site(1),site(2),site(3)) = cave/nave
	G(:,site(1),site(2),site(3)) = gave/nave
else
	C(site(1),site(2),site(3)) = 0
	G(:,site(1),site(2),site(3)) = 0
endif
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
integer :: i, x, y, z
real :: g(3), gamp, gmax
logical, save :: first = .true.

do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		if (first) then
			allocate(chemo(i)%conc(NX,NY,NZ))
			allocate(chemo(i)%grad(3,NX,NY,NZ))
			chemo(i)%conc = 0
		endif
		write(logmsg,*) 'Solving steady-state: ',chemo(i)%name
		call logger(logmsg)
		write(*,'(a,i2,2e12.3)') 'diff_coeff, decay_rate: ', i,chemo(i)%diff_coeff,chemo(i)%decay_rate
		call SolveSteadystate(i,chemo(i)%diff_coeff,chemo(i)%decay_rate,chemo(i)%conc)
		! Now compute the gradient field.
		call gradient(chemo(i)%conc,chemo(i)%grad)
		gmax = 0
		do x = 1,NX
			do y = 1,NY
				do z = 1,NZ
					g = chemo(i)%grad(:,x,y,z)
					gamp = sqrt(dot_product(g,g))
					gmax = max(gamp,gmax)
				enddo
			enddo
		enddo
		write(logmsg,*) 'Max gradient: ',gmax
		call logger(logmsg)
	endif
enddo
if (first) then
	call ShowConcs
endif
first = .false.
return

!if (use_S1P) then
!	if (first) then
!		allocate(S1P_conc(NX,NY,NZ))
!		allocate(S1P_grad(3,NX,NY,NZ))
!		S1P_conc = 0
!	endif
!	call init_S1P
!endif
!if (use_OXY) then
!	if (first) then
!		allocate(OXY_conc(NX,NY,NZ))
!		allocate(OXY_grad(3,NX,NY,NZ))
!		OXY_conc = 0
!	endif
!	call init_OXY
!endif
!if (use_CCL21) then
!	if (first) then
!		allocate(CCL21_conc(NX,NY,NZ))
!		allocate(CCL21_grad(3,NX,NY,NZ))
!		CCL21_conc = 0
!	endif
!	call init_CCL21
!endif
!if (use_CXCL13) then
!	if (first) then
!		allocate(CXCL13_conc(NX,NY,NZ))
!		allocate(CXCL13_grad(3,NX,NY,NZ))
!		CXCL13_conc = 0
!	endif
!	call init_CXCL13
!endif
if (first) then
	call ShowConcs
endif
first = .false.

!do i = 1,100
!    call evolveS1P(6,0.0)
!enddo
!call ShowConcs
end subroutine

!----------------------------------------------------------------------------------------
! Recompute steady-state concentration fields if there has been a significant change in
! the B cell population.
!----------------------------------------------------------------------------------------
subroutine UpdateFields
integer, save :: NBlast = 0
integer :: x, y, z, site(3)
real :: delNB
real :: cmin, cmax, dc, dcmax
real, allocatable :: S1P_old(:,:,:)
logical :: doit

if (NBlast == 0) then
	NBlast = NBcells0
endif
doit = .false.
if (inflammation_level == 0) then
	if (mod(istep,4*60*6) == 0) then	! every 6 hours in no inflammation case
		doit = .true.
	endif
else
	delNB = abs(NBcells - NBlast)
	if (delNB/NBcells > 0.05) then
		doit = .true.
	endif
endif
if (.not.doit) return
!cmin = 1.0e10
!cmax = 0
!do x = 1,NX
!	do y = 1,NY
!		do z = 1,NZ
!			if (occupancy(x,y,z)%indx(1) >= 0) then
!				cmin = min(cmin,S1P_conc(x,y,z))
!				cmax = max(cmin,S1P_conc(x,y,z))
!			endif
!		enddo
!	enddo
!enddo
!write(*,*) 'S1P_conc: min, max: ',cmin,cmax
!allocate(S1P_old(NX,NY,NZ))
!S1P_old = S1P_conc
call InitFields
NBlast = NBcells
!dcmax = 0
!do x = 1,NX
!	do y = 1,NY
!		do z = 1,NZ
!			if (occupancy(x,y,z)%indx(1) >= 0) then
!				dc = (abs(S1P_old(x,y,z) - S1P_conc(x,y,z)))/S1P_old(x,y,z)
!				if (dc > dcmax) then
!					dcmax = dc
!					site = (/x,y,z/)
!				endif
!			endif
!		enddo
!	enddo
!enddo
!write(*,*) 'Max fractional change in S1P_conc: ',dcmax,S1P_old(site(1),site(2),site(3)),S1P_conc(site(1),site(2),site(3))
!write(*,*) 'site: ',site
!x = site(1)
!y = site(2)
!z = site(3)
!write(*,*) 'old conc:'
!write(*,'(a,4e12.3)') 'above: ',S1P_old(x,y+1,z)
!write(*,'(a,4e12.3)') 'level: ',S1P_old(x-1,y,z),S1P_old(x+1,y,z),S1P_old(x,y,z-1),S1P_old(x,y,z+1)
!write(*,'(a,4e12.3)') 'below: ',S1P_old(x,y-1,z)
!write(*,*) 'new conc:'
!write(*,'(a,4e12.3)') 'above: ',S1P_conc(x,y+1,z)
!write(*,'(a,4e12.3)') 'level: ',S1P_conc(x-1,y,z),S1P_conc(x+1,y,z),S1P_conc(x,y,z-1),S1P_conc(x,y,z+1)
!write(*,'(a,4e12.3)') 'below: ',S1P_conc(x,y-1,z)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowConcs
integer :: x, y, z, i

x = NX/2
z = NZ/2
do i = 1,MAX_CHEMO
	if (chemo(i)%used) then
		write(nfout,'(a,a)') chemo(i)%name,'  conc        gradient' 
		do y = 1,NY
			if (occupancy(x,y,z)%indx(1) < 0) cycle
		    write(nfout,'(i4,4e12.4)') y,chemo(i)%conc(x,y,z),chemo(i)%grad(:,x,y,z)
		enddo
	endif
enddo
!do y = 1,NY
!    if (occupancy(x,y,z)%indx(1) < 0) cycle
!    write(nfout,'(i4,4e12.4)') y,S1P_conc(x,y,z),CCL21_conc(x,y,z) !,OXY_conc(x,y,z)
!enddo
!write(nfout,*) 'Grad   S1P                                 CCL21'
!do y = 1,NY
!    if (occupancy(x,y,z)%indx(1) < 0) cycle
!    write(nfout,'(i4,6e12.4)') y,S1P_grad(:,x,y,z),CCL21_grad(:,x,y,z)
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
! On entry C contains the starting concentration field, which may be 0.
!----------------------------------------------------------------------------------------
subroutine SolveSteadystate(ichemo,Kdiffusion,Kdecay,C)
integer :: ichemo
real :: C(:,:,:)
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, maxchange_par, total_par
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-5		! max change in C at any site as fraction of average C
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

do it = 1,nt
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
! A different approach from that used in bone-abm.
! The bdry gridcells that are adjacent to a chemokine-rich region (i.e. there are 
! chemokine-secreting cells near the boundary) are given a fixed concentration.
! All boundaries are no-flux.
! This is a simplified approach.  We need only keep track of the boundary.
! The concentrations in a specified chemokine influx bdry site are input parameters.
! Note:  It is only the ratio Kdecay/Kdiffusion that matters.
!----------------------------------------------------------------------------------------
subroutine par_steadystate(ichemo,C,Kdiffusion,Kdecay,Ctemp,x1,x2,maxchange,total,nc,kpar)
integer :: ichemo
real :: C(:,:,:), Ctemp(:,:,:)
integer :: x1, x2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, xpar, indx(2), i
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
                do i = 1,MAX_CHEMO
	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(xpar,y,z) = C(x,y,z)
				        bdry_conc = .true.
					endif
				enddo
                
!                if (ichemo == S1P .and. occupancy(x,y,z)%bdry%S1P) then
!                    C(x,y,z) = BdryS1PConc
!                    Ctemp(xpar,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == OXY .and. occupancy(x,y,z)%bdry%OXY) then
!                    C(x,y,z) = BdryOXYConc
!                    Ctemp(xpar,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == CCL21 .and. occupancy(x,y,z)%bdry%CCL21) then
!                    C(x,y,z) = BdryCCL21Conc
!                    Ctemp(xpar,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == CXCL13 .and. occupancy(x,y,z)%bdry%CXCL13) then
!                    C(x,y,z) = BdryCXCL13Conc
!                    Ctemp(xpar,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
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
real :: C(:,:,:), grad(:,:,:,:)
integer :: x, y, z, xx, yy, zz, x1, x2, y1, y2, z1, z2, i, k
real :: g(3)
logical :: missed
real, parameter :: MISSING_VAL = 1.0e10

grad = 0
do z = 1,NZ
	do y = 1,NY
		do x = 1,NX
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			x1 = x - 1
			x2 = x + 1
			if (x1 < 1 .or. x2 > NX) then
				g(1) = 0
			elseif (occupancy(x1,y,z)%indx(1) >= 0 .and. occupancy(x2,y,z)%indx(1) >= 0) then
				g(1) = (C(x2,y,z) - C(x1,y,z))/(2*DELTA_X)
			else
				g(1) = MISSING_VAL
			endif
			y1 = y - 1
			y2 = y + 1
			if (y1 < 1 .or. y2 > NY) then
				g(2) = 0
			elseif (occupancy(x,y1,z)%indx(1) >= 0 .and. occupancy(x,y2,z)%indx(1) >= 0) then
				g(2) = (C(x,y2,z) - C(x,y1,z))/(2*DELTA_X)
			else
				g(2) = MISSING_VAL
			endif
			z1 = z - 1
			z2 = z + 1
			if (z1 < 1 .or. z2 > NZ) then
				g(3) = 0
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
integer :: x, y, z, xx, yy, zz, nb, nc, k, it, i
real, allocatable :: Ctemp(:,:,:)
logical :: bdry_conc

dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
do z = 1,NZ
	do y = 1,NY
		do x = 1,NX
			if (occupancy(x,y,z)%indx(1) < 0) cycle
			C0 = C(x,y,z)
            bdry_conc = .false.
            if (associated(occupancy(x,y,z)%bdry)) then
                ! Check for S1P bdry site - no change to the concentration at such a site
                do i = 1,MAX_CHEMO
	                if (ichemo == i .and. occupancy(x,y,z)%bdry%chemo_influx(i)) then
		                C(x,y,z) = chemo(i)%bdry_conc
			            Ctemp(x,y,z) = C(x,y,z)
				        bdry_conc = .true.
					endif
				enddo

!                if (ichemo == S1P .and. occupancy(x,y,z)%bdry%S1P) then
!                    C(x,y,z) = BdryS1PConc
!                    Ctemp(x,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == OXY .and. occupancy(x,y,z)%bdry%OXY) then
!                    C(x,y,z) = BdryOXYConc
!                    Ctemp(x,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == CCL21 .and. occupancy(x,y,z)%bdry%CCL21) then
!                    C(x,y,z) = BdryCCL21Conc
!                    Ctemp(x,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
!                if (ichemo == CXCL13 .and. occupancy(x,y,z)%bdry%CXCL13) then
!                    C(x,y,z) = BdryCXCL13Conc
!                    Ctemp(x,y,z) = C(x,y,z)
!                    bdry_conc = .true.
!                endif
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