! To handle FDCs

module FDC

use global
implicit none
save

contains

!--------------------------------------------------------------------------------
! The simulation currently starts with FDCs already clustered in a region of the
! follicle that is below the centre, i.e. closer to the T zone interface.
! The FDCs are assumed to be close together, but not so close that the lymphocytes
! cannot move among them.
! In the absence of a better idea, the FDCs are given the same shape as the DCs,
! i.e. they occupy 7 sites in the 6-point 3D von Neumann configuration.
! Start with an FDC at the centre of the desired GC region, then add others at
! nearby locations.
!--------------------------------------------------------------------------------
subroutine placeFDCs(ok)
logical :: ok
integer :: NFDCrequired, site(3), nassigned, ifdc, ilim, dx, dy, dz, d2, x, y, z, k, kpar=0
integer, parameter :: MAX_LIM = 30
logical :: success, done, checked(-MAX_LIM:MAX_LIM,-MAX_LIM:MAX_LIM,-MAX_LIM:MAX_LIM)
real :: R, separation_limit, GC_centre(3)
real :: FDC_SEPARATION
type(vector3_type) :: eradius
integer, parameter :: method = 3	! 1, 2 or 3

if (.not.use_FDCs) then
	NFDC = 0
	ok = .true.
	return
endif
k = 0
DCoffset(:,1) = 0
DCoffset(:,2) = (/-1,0,0/)
DCoffset(:,3) = (/1,0,0/)
DCoffset(:,4) = (/0,-1,0/)
DCoffset(:,5) = (/0,1,0/)
DCoffset(:,6) = (/0,0,-1/)
DCoffset(:,7) = (/0,0,1/)
NFDCrequired = BASE_NFDC
if (method <= 2) then
    ! Place the first FDC at a point midway between the blob centre and the T zone bdry
    GC_centre = Centre + (/0.,-Radius%y/2, 0./)
    FDC_SEPARATION = 1.8
elseif (method == 3) then
    GC_centre = Centre + (/0.,Radius%y/3, 0./)
    eradius%y = (2./3.)*Radius%y
    eradius%x = 1.5*ELLIPSE_RATIO*eradius%y
    eradius%z = eradius%x
    FDC_SEPARATION = 2.8
endif
write(logmsg,*) 'GC_centre: ',GC_centre
call logger(logmsg)
FDClist(1)%ID = 1
FDClist(1)%site = GC_centre
FDClist(1)%alive = .true.
call AssignFDCsites(1,nassigned,ok)
if (.not.ok) then
	return
endif
FDClist(1)%nsites = nassigned
NFDC = 1
nlist = nlist - nassigned
write(logmsg,*) 'FDC site: ',1,FDClist(1)%site,nassigned
call logger(logmsg)
! Now place the rest of the FDCs.
! The criteria of FDC location are:
!  for methods 1 and 2:
!    near the proposed GC centre
!    not nearer than FDC_SEPARATION (number of sites) from any other FDC
!   (method 1 generates a rather cubic FDC region)
!  for method 3:
!    within an ellipsoid with specified:
!    centre
!    major radius (xy plane)
!    minor radius (z axis)
!
if (method == 1) then
	checked = .false.
	ilim = 1
	outer_loop: do
		if (ilim == 10) then
			call logger('Error: PlaceFDCs: failed to place all FDCs')
			ok = .false.
			return
		endif
		ilim = ilim + 1
		do dz = -ilim,ilim
			do dy = -ilim, ilim
				do dx = -ilim,ilim
					if (checked(dx,dy,dz)) cycle
					checked(dx,dy,dz) = .true.
					site = GC_centre + (/dx,dy,dz/)
					R = par_uni(kpar)
					separation_limit = FDC_SEPARATION*(0.8 + 0.3*R)
					if (FDCSiteAllowed(site,separation_limit)) then
						NFDC = NFDC + 1
						FDClist(NFDC)%ID = NFDC
						FDClist(NFDC)%site = site
						FDClist(NFDC)%alive = .true.
						call AssignFDCsites(NFDC,nassigned,ok)
						if (.not.ok) then
							return
						endif
						FDClist(NFDC)%nsites = nassigned
						nlist = nlist - nassigned
						if (NFDC == NFDCrequired) exit outer_loop
					endif
				enddo
			enddo
		enddo
	enddo outer_loop
elseif (method == 2) then
	done = .false.
	checked = .false.
	checked(0,0,0) = .true.
	ilim = 2
	do
		ilim = ilim + 1
		if (ilim > MAX_LIM) then
		    write(logmsg,*) 'Error: PlaceFDCs: ilim exceeds MAX_LIM: reduce FDC_SEPARATION'
		    call logger(logmsg)
		    ok = .false.
		    return
		endif
		do k = 1,1000
			dx = random_int(-ilim,ilim,kpar)
			dy = random_int(-ilim,ilim,kpar)
			dz = random_int(-ilim,ilim,kpar)
			if (checked(dx,dy,dz)) cycle
			d2 = dx*dx + dy*dy + dz*dz
			if (d2 > ilim*ilim) cycle
			site = GC_centre + (/dx,dy,dz/)
			if (.not.InsideEllipsoid(site,Centre,Radius)) cycle
			R = par_uni(kpar)
			separation_limit = FDC_SEPARATION*(0.8 + 0.3*R)
			if (FDCSiteAllowed(site,separation_limit)) then
				NFDC = NFDC + 1
				FDClist(NFDC)%ID = NFDC
				FDClist(NFDC)%site = site
				FDClist(NFDC)%alive = .true.
				call AssignFDCsites(NFDC,nassigned,ok)
				if (.not.ok) then
					return
				endif
				FDClist(NFDC)%nsites = nassigned
				nlist = nlist - nassigned
				checked(dx,dy,dz) = .true.
				write(logmsg,'(a,10i4)') 'FDC site: ilim: ',NFDC,ilim,dx,dy,dz,NFDC,FDClist(NFDC)%site,nassigned
				call logger(logmsg)
				if (NFDC == NFDCrequired) then
					done = .true.
					exit
				endif
			endif				
		enddo
		if (done) exit
	enddo
elseif (method == 3) then
	done = .false.
	checked = .false.
	checked(0,0,0) = .true.
	ilim = 2
	do
		ilim = ilim + 1
		if (ilim > MAX_LIM) then
		    write(logmsg,*) 'Error: PlaceFDCs: ilim exceeds MAX_LIM: reduce FDC_SEPARATION'
		    call logger(logmsg)
		    ok = .false.
		    return
		endif
		do k = 1,1000
			dx = random_int(-ilim,ilim,kpar)
			dy = random_int(-ilim,ilim,kpar)
			dz = random_int(-ilim,ilim,kpar)
			if (checked(dx,dy,dz)) cycle
			d2 = dx*dx + dy*dy + dz*dz
			if (d2 > ilim*ilim) cycle
			site = GC_centre + (/dx,dy,dz/)
			if (.not.InsideEllipsoid(site,Centre,Radius)) cycle
			if (.not.InsideEllipsoid(site,GC_centre,eradius)) cycle
			R = par_uni(kpar)
			separation_limit = FDC_SEPARATION*(0.8 + 0.3*R)
			if (FDCSiteAllowed(site,separation_limit)) then
				NFDC = NFDC + 1
				FDClist(NFDC)%ID = NFDC
				FDClist(NFDC)%site = site
				FDClist(NFDC)%alive = .true.
				call AssignFDCsites(NFDC,nassigned,ok)
				if (.not.ok) then
					return
				endif
				FDClist(NFDC)%nsites = nassigned
				nlist = nlist - nassigned
				checked(dx,dy,dz) = .true.
				write(logmsg,'(a,9i4)') 'FDC site: ilim: ',ilim,dx,dy,dz,NFDC,FDClist(NFDC)%site,nassigned
				call logger(logmsg)
				if (NFDC == NFDCrequired) then
					done = .true.
					exit
				endif
			endif				
		enddo
		if (done) exit
	enddo
endif

call AssignFDCBdrySites

write(logmsg,'(a,2i6)') 'Number of FDCs, B cells: ',NFDC,nlist
call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Check that a proposed FDC site is far enough from all existing FDCs
!--------------------------------------------------------------------------------
logical function FDCSiteAllowed(site0,rlimit)
integer :: site0(3)
real :: rlimit
integer :: ifdc, site(3), k
real :: r(3), d

do ifdc = 1,NFDC
	r = FDClist(ifdc)%site - site0
	d = norm(r)
	if (d < rlimit) then
		FDCSiteAllowed = .false.
		return
	endif
enddo
do k = 1,NDCsites
	site = site0 + DCoffset(:,k)
	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) then
		FDCSiteAllowed = .false.
		return
	endif
enddo
FDCSiteAllowed = .true.
end function

!--------------------------------------------------------------------------------
! Sites adjacent to an FDC are flagged - these become boundary sites with
! specified chemokine concentration in SolveSteadystate_A(),
! or sites that receive FDC chemokine secretion at specified rates in
! SolveSteadystate_B().  
! occupancy(:,:,:)%FDC_nbdry records the number of adjacent FDCs
!--------------------------------------------------------------------------------
subroutine AssignFDCBdrySites
integer :: i, j, ifdc, site(3), fsite(3), bsite(3), dx, dy, dz
logical :: bdry(-2:2,-2:2,-2:2)

occupancy(:,:,:)%FDC_nbdry = 0

do ifdc = 1,NFDC
    bdry = .false.
    do i = 2,7
	    fsite = DCoffset(:,i)
	    do j = 1,27
		    if (j == 14) cycle
		    bsite = fsite + jumpvec(:,j)
		    bdry(bsite(1),bsite(2),bsite(3)) = .true.
	    enddo
    enddo
    do dx = -2,2
	    do dy = -2,2
		    do dz = -2,2
			    if (bdry(dx,dy,dz)) then
				    site = FDClist(ifdc)%site + (/dx,dy,dz/)
				    occupancy(site(1),site(2),site(3))%FDC_nbdry = occupancy(site(1),site(2),site(3))%FDC_nbdry + 1
			    endif
		    enddo
	    enddo
    enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! All associations of sites with DC are recomputed.
! That is for every site (x,y,z), create the list of DC that are near this site:
! occupancy(x,y,z)%DC(1:3), where occupancy(x,y,z)%DC(0) = number of nearby DC (if >= 0)
! To avoid having to explicitly select the closest DC (when there are more than DCDIM-1 near
! a site), the order of scanning the DClist is randomized.
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine reassign_DC(kpar,ok)
integer :: kpar
logical :: ok
integer, allocatable :: perm(:)
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax
integer :: idc, kdc, k, x, y, z, y2, z2, d2, site(3), nassigned
logical :: added


occupancy(:,:,:)%DC(0) = 0
allocate(perm(MAX_DC))
do k = 1, NDC
    perm(k) = k
enddo

call permute(perm,NDC,kpar)
NDCalive = 0
do k = 1,NDC
    idc = perm(k)
    if (.not.DClist(idc)%alive) cycle
	if (DClist(idc)%nsites < NDCsites) then
		call AssignFDCsites(idc,nassigned,ok)
		if (.not.ok) return
		DClist(idc)%nsites = nassigned
	endif
    NDCalive = NDCalive + 1
    site = DClist(idc)%site
    xdc = site(1)
    ydc = site(2)
    zdc = site(3)
    xmin = xdc - DCRadius
    xmax = xdc + DCRadius
    ymin = ydc - DCRadius
    ymax = ydc + DCRadius
    zmin = zdc - DCRadius
    zmax = zdc + DCRadius
    xmin = max(1,xmin)
    xmax = min(NX,xmax)
    ymin = max(1,ymin)
    ymax = min(NY,ymax)
    zmin = max(1,zmin)
    zmax = min(NZ,zmax)
    do z = zmin,zmax
        z2 = (z-zdc)*(z-zdc)
        do y = ymin,ymax
            y2 = (y-ydc)*(y-ydc)
            do x = xmin,xmax
                d2 = (x-xdc)*(x-xdc) + y2 + z2
                added = .false.
	            if (d2 <= DCRadius*DCRadius) then
	                kdc = occupancy(x,y,z)%DC(0)
	                if (kdc < 0) cycle			! this is a DC site
		            if (kdc < DCDIM-1) then		! not all possible DC assigned
						kdc = kdc + 1
						occupancy(x,y,z)%DC(kdc) = idc
						occupancy(x,y,z)%DC(0) = kdc
						added = .true.
					endif
				endif
            enddo
        enddo
    enddo
enddo
deallocate(perm)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Test moveDC() with a randomly selected DC, random step direction.
!-----------------------------------------------------------------------------------------
subroutine test_moveDC
integer :: idc, dir
integer :: kpar=0

write(*,*) 'test_moveDC'
if (NDCsites /= 7) then
    write(*,*) 'Valid only for NDCsites = 7'
    stop
endif
do
    idc = random_int(1,NDC,kpar)
    if (DClist(idc)%alive) exit
enddo
dir = random_int(1,6,kpar)
call moveDC(idc,dir)
end subroutine

!-----------------------------------------------------------------------------------------
! A DC is allowed to move only if it has grown to its full extent,
! i.e. if DClist()%nsites = NDCsites
! For now allow only moves in the 6 principal directions (Neumann).
! Currently valid only for NDCsites = 7.
! Cells in the 5 sites in the path of the DC step are moved.
! For 4 sites, the shift is by one site, for the one site in line with the DC centre the
! shift is 3 sites.
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine moveDC(idc,dir)
integer :: idc, dir
integer :: k, i, kcell, indx(2), site0(3), site1(3), site2(3), step(3)

if (DClist(idc)%nsites /= NDCsites) return
step = neumann(:,dir)
site0 = DClist(idc)%site
do k = 2,NDCsites
    if (all(DCoffset(:,k) == step)) then       ! this is in the direction of step
        ! move site contents by 3 sites in opposite direction to step
        site1 = site0 - step       ! old DC site
        site2 = site0 + 2*step     ! new DC site
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
        occupancy(site0(1),site0(2),site0(3))%indx = -idc
        site2 = site0  + step
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    elseif (all(DCoffset(:,k) == -step)) then   ! this is the direction opposite to step
        ! do nothing
    else                                        ! one of the other 4 directions
        ! move site contents by 1 site in opposite direction to step
        site1 = site0 + DCoffset(:,k)       ! old DC site, new T cell site
        site2 = site1 + step                ! new DC site, old T cell site
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    endif
enddo
DClist(idc)%site = site0 + step
end subroutine

!-----------------------------------------------------------------------------------------
! NOT USED
!-----------------------------------------------------------------------------------------
subroutine check_DCproximity
integer :: x,y,z,k,cnt(0:DCDIM-1)
integer :: fdc(DCDIM-1)

cnt = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            k = occupancy(x,y,z)%DC(0)
            if (k >= 0) then
                cnt(k) = cnt(k) + 1
                fdc = occupancy(x,y,z)%DC(1:DCDIM-1)
                if (k >= 2) then
                    call randomizeDC(fdc,k)
                    occupancy(x,y,z)%DC(1:DCDIM-1) = fdc
                endif
            endif
        enddo
    enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine randomizeDC(dc,k)
integer :: k, dc(:)
integer :: i,p(DCDIM-1),tmp(DCDIM-1)
integer :: kpar=0

do i = 1,k
    p(i) = i
enddo
call permute(p,k,kpar)
tmp(1:k) = dc(1:k)
do i = 1,k
    dc(i) = tmp(p(i))
enddo
!write(*,*) tmp(1:k),dc(1:k)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AssignFDCsites(ifdc,nassigned,ok)
integer :: ifdc, nassigned
logical :: ok
integer :: k, site(3), site1(3)

ok = .false.
site = FDClist(ifdc)%site
occupancy(site(1),site(2),site(3))%indx = -ifdc     ! This site holds an FDC (centre)
nassigned = 1
do k = 2,NDCsites
    site1 = site + DCoffset(:,k)
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == -ifdc) then
		nassigned = nassigned + 1
		cycle
	endif
    ! Reduce available blob sites only if the site is within the blob.
    ! A DC site should never fall outside, I think.  Check this.
    if (.not.inside_xyz(site1)) then
		call logger("Error: AssignFDCsites: site outside grid")
		return
	endif
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) then
		call logger('AssignFDCsites: FDC hits boundary')
		return
	endif
	if (occupancy(site1(1),site1(2),site1(3))%indx(1) /= 0) cycle
	if (occupancy(site1(1),site1(2),site1(3))%indx(2) /= 0) cycle
	nassigned = nassigned + 1
    occupancy(site1(1),site1(2),site1(3))%indx = -ifdc     ! This site holds a FDC (soma)
enddo
ok = .true.
end subroutine

end module
