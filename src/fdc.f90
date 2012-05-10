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
integer :: NFDCrequired, site(3), GC_centre(3), nassigned, ifdc, ilim, dx, dy, dz, x, y, z, k
logical :: success, checked(-10:10,-10:10,-10:10)
real, parameter :: FDC_SEPARATION = 4

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
! Place the first FDC at a point midway between the blob centre and the T zone bdry
GC_centre = Centre + (/0.,-bRadius/2, 0./)
write(logmsg,*) 'GC_centre: ',GC_centre
call logger(logmsg)
FDClist(1)%ID = 1
FDClist(1)%site = GC_centre
FDClist(1)%alive = .true.
call AssignFDCsites(1,nassigned,ok)
if (.not.ok) then
	return
endif
!write(*,*) 'nassigned: ',nassigned
FDClist(1)%nsites = nassigned
NFDC = 1
nlist = nlist - nassigned
write(logmsg,*) 'FDC site: ',1,FDClist(1)%site,nassigned
call logger(logmsg)
! Now place the rest of the FDCs.
! The criteria of FDC location are:
!   near the proposed GC centre
!   not nearer than FDC_SEPARATION (number of sites) from any other FDC
checked = .false.
ilim = 1
outer_loop: do
	ilim = ilim + 1
	do dz = -ilim,ilim
		do dy = -ilim, ilim
			do dx = -ilim,ilim
				if (checked(dx,dy,dz)) cycle
				site = GC_centre + (/dx,dy,dz/)
				checked(dx,dy,dz) = .true.
				if (FDCSiteAllowed(site,FDC_SEPARATION)) then
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
					write(logmsg,'(a,9i4)') 'FDC site: ilim: ',ilim,dx,dy,dz,NFDC,FDClist(NFDC)%site,nassigned
					call logger(logmsg)
					if (NFDC == NFDCrequired) exit outer_loop
				endif
			enddo
		enddo
	enddo
enddo outer_loop

do ifdc = 1,NFDC
	call AssignFDCBdry(ifdc)
enddo

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
subroutine AssignFDCBdry(ifdc)
integer :: ifdc
integer :: i, j, site(3), fsite(3), bsite(3), dx, dy, dz
logical :: bdry(-2:2,-2:2,-2:2)

occupancy(:,:,:)%FDC_nbdry = 0
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
end subroutine


!--------------------------------------------------------------------------------
! Place n DCs
! Note that we may have NDCsites > 1 !!!
! The procedure for placing DCs is as follows:
! The central DC is placed on a site that satisfies these conditions:
! (1) The site is not OUTSIDE_TAG or a DC site
! (2) There are not two T cells at this site
! (3) The site is not too near another DC, i.e. it is more than DC_DCprox*DCRadius
!     from any DC.
! (4) The site is not too near an exit, i.e. it is more than exit_DCprox from any exit.
! (5) The site is not too near the blob boundary, i.e. it is more than
!     bdry_FDCprox*FDCRadius from the sphere with radius = Radius
! (6) Either the site is free, or it is occupied by a T cell that can be moved
!     to a neighbouring site.
! A site meeting these requirements is selected for the DC, then as many as
! possible of the neighbouring sites are also allocated to the NDCsites-1 other
! sites occupied by a DC, by subroutine addDCsite().  The count of DC sites
! allocated is stored in DC%nsites.
! NOTE: This is just the code from the paracortex model for DC placement.
!--------------------------------------------------------------------------------
subroutine place_FDCs1(n,nadded)
integer :: n, nadded
integer :: i, x, y, z, site1(3), site2(3), freeslot, err
integer :: indx(2), jslot, idc, k, kcell
integer :: xmin, xmax, ymin, ymax, zmin, zmax
integer :: kpar = 0
real(DP) :: R
real :: tnow, dist, rvec(3), prox, tmins
logical :: OK
type(FDC_type) :: FDC
integer, parameter :: kmax = 100000

tnow = istep*DELTA_T
xmin = x0 - aRadius
xmax = x0 + aRadius + 1
ymin = y0 - aRadius
ymax = y0 + aRadius + 1
zmin = z0 - aRadius
zmax = z0 + aRadius + 1
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
nadded = 0
do i = 1,n
	OK = .true.
	k = 0
    do
		k = k+1
		if (k > kmax) then
			write(logmsg,*) 'Error: place_DCs: unable to find space for a DC'
			call logger(logmsg)
			OK = .false.
			stop
		endif
        R = par_uni(kpar)
        x = xmin + R*(xmax-xmin)
        R = par_uni(kpar)
        y = ymin + R*(ymax-ymin)
        R = par_uni(kpar)
        z = zmin + R*(zmax-zmin)
        site1 = (/x,y,z/)   ! global location
        indx = occupancy(x,y,z)%indx
        if (indx(1) < 0) cycle                          ! OUTSIDE_TAG or DC
        if (indx(1) /= 0 .and. indx(2) /= 0) cycle      ! two T cells at this site
        if (tooNearFDC(site1,NDC,DC_DCprox*DCRadius)) cycle
!        if (tooNearExit(site1,exit_DCprox)) cycle
!        prox = 0.5*DC_DCprox*DCRadius
        prox = bdry_FDCprox*FDCRadius
        if (.not.tooNearBdry(site1,prox)) then
            jslot = 0
            if (indx(1) /= 0) then
                jslot = 1
            elseif (indx(2) /= 0) then
                jslot = 2
            endif
            if (jslot == 0) then ! free site, use it
                exit
            else  ! one T cell here in jslot, it must be moved
                ! Can the T cell be bumped to a neighbour site?
                call get_free_slot(NX,site1,site2,freeslot)
                if (freeslot == 0) cycle    ! cannot be bumped, try again
               ! Move the cell in site1/jslot to site2/freeslot
                kcell = indx(jslot)
!                write(*,*) 'place_DCs: bump cell: ',kcell,site1,site2
                occupancy(x,y,z)%indx = 0
                occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
                cellist(kcell)%site = site2
                ! Now site1 is free to use
                exit
            endif
        endif
    enddo
    if (.not.OK) exit
    idc = 0
!    if (reuse_DC_index) then
!	    do k = 1,NDC
!		    if (.not.DClist(k)%alive) then
!		        idc = k
!		        exit
!		    endif
!		enddo
!	endif
    if (idc == 0) then    ! If there isn't a free spot in FDClist()
        NFDC = NFDC + 1
        if (NFDC > MAX_FDC) then
			write(logmsg,'(a,i6)') 'Error: place_FDCs: number of FDCs exceeds limit: ',max_FDC
			call logger(logmsg)
			ok = .false.
			return
		endif
        idc = NFDC
    endif
    FDC%ID = idc
    FDC%alive = .true.
!    FDC%capable = .true.
    FDC%site = site1
    FDC%nsites = 1
!    FDC%stimulation = 0
!    FDC%nbound = 0
!    FDC%ncogbound = 0
!    DC%density = DCdensity(kpar)
!    DC%dietime = tnow + DClifetime(kpar)
    occupancy(site1(1),site1(2),site1(3))%indx(1) = -idc
    do k = 2,NDCsites
        call addFDCsite1(idc,site1,k,err)
        if (err == 0) then
            FDC%nsites = FDC%nsites + 1
        else
!            write(*,*) 'addDCsite: idc,k,err: ',idc,k,err
        endif
    enddo
    FDClist(idc) = FDC
    nadded = nadded + 1
!    write(*,*) 'Added DC at: ',site1,' with nsites: ',DC%nsites
    ! now the DC proximity data in occupancy()%DC must be updated
    ! this is done by reassign_DC() called from balancer()
!    if (DClist(idbug)%ncogbound /= ndbug) then
!        write(*,*) 'place_DCs (c): ndbug changed: ',i,DClist(idbug)%ncogbound
!        stop
!    endif
enddo

! Now check the DC locations wrt the blob perimeter
do k = 1,NFDC
    if (FDClist(k)%alive) then
        rvec = FDClist(k)%site - (/x0,y0,z0/)
        dist = norm(rvec)
        if (dist > aRadius) then
            write(logmsg,*) 'Place_FDCs: warning: FDC distance: ',k,dist,aRadius
			call logger(logmsg)
        endif
    endif
enddo
!write(*,*) 'place_DCs (2): ',DClist(idbug)%ncogbound
end subroutine


!-----------------------------------------------------------------------------------------
! Convert one of the sites near a DC into a DC peripheral site for DC idc.
! The index k indicates which peripheral site of 2:NDCsites to convert.
! A site is a candidate for a DC site provided:
! (1) It is not outside the blob, or already a DC site
! (2) It does not hold two T cells
! (3) It is not too close the the blob boundary
! If a candidate site is free, it is used.  If it holds a single T cell, the cell
! is moved if possible to a neighbouring site.
! Note that when a T cell is bumped it retains its binding (if any) to a DC, even though
! it may have been moved to a site that is - strictly speaking - not within the SOI of
! the DC.
!-----------------------------------------------------------------------------------------
subroutine addFDCsite1(idc,site0,k,err)
integer :: idc, site0(3), k, err
integer :: indx(2), jslot, kcell, freeslot, site1(3), site2(3)
real :: prox
logical :: OK

site1 = site0 + DCoffset(:,k)
indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) < 0) then                       ! OUTSIDE_TAG or DC
    err = 1
    return
endif
if (indx(1) /= 0 .and. indx(2) /= 0) then   ! two T cells at this site
    err = 2
    return
endif
prox = bdry_FDCprox*FDCRadius
if (tooNearBdry(site1,prox)) then                ! too close to the blob boundary
    err = 4
    return
endif

OK = .false.
jslot = 0
if (indx(1) /= 0) then
    jslot = 1
elseif (indx(2) /= 0) then
    jslot = 2
endif
if (jslot == 0) then ! free site, use it
    OK = .true.
else  ! one T cell here in jslot, it must be moved
    ! Can the T cell be bumped to a neighbour site?
    call get_free_slot(NX,site1,site2,freeslot)
    if (freeslot == 0) then    ! cannot be bumped
        err = 5
        return
    endif
   ! Move the cell in site1/jslot to site2/freeslot
    kcell = indx(jslot)
    occupancy(site1(1),site1(2),site1(3))%indx = 0
    occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
    cellist(kcell)%site = site2
    ! Now site1 is free to use
    OK = .true.
endif
!if (OK) then
    err = 0
    occupancy(site1(1),site1(2),site1(3))%indx = -idc
!else
!    err = 6
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! All associations of sites with DC are recomputed.
! That is for every site (x,y,z), create the list of DC that are near this site:
! occupancy(x,y,z)%DC(1:3), where occupancy(x,y,z)%DC(0) = number of nearby DC (if >= 0)
! To avoid having to explicitly select the closest DC (when there are more than DCDIM-1 near
! a site), the order of scanning the DClist is randomized.
!-----------------------------------------------------------------------------------------
subroutine reassign_DC(kpar,ok)
integer :: kpar
logical :: ok
integer, allocatable :: perm(:)
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax
integer :: idc, kdc, k, x, y, z, y2, z2, d2, site(3), nassigned
logical :: added

!write(*,*) 'reassign_DC'

occupancy(:,:,:)%DC(0) = 0
!occupancy(:,:,:)%cDC(0) = 0
allocate(perm(MAX_DC))
do k = 1, NDC
    perm(k) = k
enddo

call permute(perm,NDC,kpar)
NDCalive = 0
do k = 1,NDC
    idc = perm(k)
    if (.not.DClist(idc)%alive) cycle
!    write(*,*) 'idc: ',idc
!	write(*,*) '(4) cell 16243: ',cellist(16243)%site,occupancy(70,63,88)%indx
	if (DClist(idc)%nsites < NDCsites) then
		call AssignFDCsites(idc,nassigned,ok)
		if (.not.ok) return
!		write(*,*) 'assigned DC sites for: ',idc,nassigned
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
                ! The following procedure is correct only if DCRadius <= chemo_radius
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
!	            if (.not.added .and. (d2 <= chemo_radius*chemo_radius)) then
!	                kdc = occupancy(x,y,z)%cDC(0)
!		            if (kdc == cDCDIM-1) cycle   ! all possible DC assigned
!		            kdc = kdc + 1
!		            occupancy(x,y,z)%cDC(kdc) = idc
!	                occupancy(x,y,z)%cDC(0) = kdc
!	            endif
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
!-----------------------------------------------------------------------------------------
subroutine moveDC(idc,dir)
integer :: idc, dir
integer :: k, i, kcell, indx(2), site0(3), site1(3), site2(3), step(3)

if (DClist(idc)%nsites /= NDCsites) return
step = neumann(:,dir)
site0 = DClist(idc)%site
!write(*,*) 'moveDC: ',idc,dir,'  ',site0,'  ',step
do k = 2,NDCsites
    if (all(DCoffset(:,k) == step)) then       ! this is in the direction of step
        ! move site contents by 3 sites in opposite direction to step
        site1 = site0 - step       ! old DC site
        site2 = site0 + 2*step     ! new DC site
!        write(*,'(a,7i4)') 'step dir: ',k,site1,site2
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
!        write(*,'(a,7i4)') 'other dir: ',k,site1,site2
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
!write(*,*) 'DC proximity counts: ',cnt

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

!!-----------------------------------------------------------------------------------------
!! Is site near a DC?
!! The criterion for a DC site might be different from an exit site.
!! prox = DC_DCprox*DCRadius for DC - DC
!!-----------------------------------------------------------------------------------------
!logical function tooNearDC(site,kdc,prox)
!integer :: site(3), kdc
!real :: prox
!integer :: idc
!real :: r(3), d
!
!if (kdc == 0) then
!    tooNearDC = .false.
!    return
!endif
!do idc = 1,kdc
!    if (.not.DClist(idc)%alive) cycle
!    r = site - DClist(idc)%site
!    d = norm(r)     ! units sites
!    if (d < prox) then
!        tooNearDC = .true.
!        return
!    endif
!enddo
!tooNearDC = .false.
!end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AssignFDCsites(ifdc,nassigned,ok)
integer :: ifdc, nassigned
logical :: ok
integer :: k, site(3), site1(3)

ok = .false.
site = FDClist(ifdc)%site
occupancy(site(1),site(2),site(3))%indx = -ifdc     ! This site holds an FDC (centre)
!occupancy(site(1),site(2),site(3))%DC(0) = -ifdc    ! This site holds a DC
nassigned = 1
!write(*,*) 'assignDCsites: ',site
do k = 2,NDCsites
    site1 = site + DCoffset(:,k)
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == -ifdc) then
!		write(*,'(a,6i4)') 'AssignDCsites: ???: site,site1: ',site,site1 
		nassigned = nassigned + 1
		cycle
	endif
    ! Reduce available blob sites only if the site is within the blob.
    ! A DC site should never fall outside, I think.  Check this.
    if (.not.inside_xyz(site1)) then
!		write(*,*) 'Error: AssignFDCsites: site outside grid: '
!		write(*,*) 'FDC site: ',site
!		write(*,*) 'DCoffset: ',k,DCoffset(:,k)
!		write(*,*) 'site1: ',site1
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
!    occupancy(site1(1),site1(2),site1(3))%DC(0) = -ifdc    ! This site holds a FDC
enddo
ok = .true.
end subroutine

end module
