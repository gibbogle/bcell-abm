! To handle FDCs

module FDC

use global
implicit none
save

integer, parameter :: FDC_TAG = 1
integer, parameter :: MRC_TAG = 2

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
real :: R, separation_limit(2,2), GC_centre(3)
real :: FDC_SEPARATION, FDC_MRC_SEPARATION
type(vector3_type) :: eradius
integer, parameter :: method = 3	! 1, 2 or 3

if (.not.use_FDCs) then
	NFDC = 0
	ok = .true.
	return
endif
call logger('PlaceFDCs')
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
    FDC_MRC_SEPARATION = 1.8
elseif (method == 3) then
	! The FDCs are located in an ellisoidal region above the blob centre.  This region is more flattened
	! than the blob ellipsoid (the factor 1.5 in x and z radii)
    GC_centre = Centre + (/0.,Radius%y/3, 0./)
    eradius%y = (2./3.)*Radius%y
    eradius%x = 1.5*ELLIPSE_RATIO*eradius%y
    eradius%z = eradius%x
    FDC_SEPARATION = 2.8
    FDC_MRC_SEPARATION = 2.0
endif
write(logmsg,*) 'GC_centre: ',GC_centre
call logger(logmsg)
FDC_list(1)%ID = 1
FDC_list(1)%site = GC_centre
FDC_list(1)%alive = .true.
call AssignFDCsites(FDC_TAG,1,nassigned,ok)
if (.not.ok) then
	return
endif
FDC_list(1)%nsites = nassigned
NFDC = 1
nBlist = nBlist - nassigned
write(logmsg,*) 'FDC site: ',1,FDC_list(1)%site,nassigned
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
					separation_limit(FDC_TAG,FDC_TAG) = FDC_SEPARATION*(0.8 + 0.3*R)
					separation_limit(FDC_TAG,MRC_TAG) = FDC_MRC_SEPARATION*(0.8 + 0.3*R)
					if (FDCSiteAllowed(FDC_TAG,site,separation_limit)) then
						NFDC = NFDC + 1
						FDC_list(NFDC)%ID = NFDC
						FDC_list(NFDC)%site = site
						FDC_list(NFDC)%alive = .true.
						call AssignFDCsites(FDC_TAG,NFDC,nassigned,ok)
						if (.not.ok) then
							return
						endif
						FDC_list(NFDC)%nsites = nassigned
						nBlist = nBlist - nassigned
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
			separation_limit(FDC_TAG,FDC_TAG) = FDC_SEPARATION*(0.8 + 0.3*R)
			separation_limit(FDC_TAG,MRC_TAG) = FDC_MRC_SEPARATION*(0.8 + 0.3*R)
			if (FDCSiteAllowed(FDC_TAG,site,separation_limit)) then
				NFDC = NFDC + 1
				FDC_list(NFDC)%ID = NFDC
				FDC_list(NFDC)%site = site
				FDC_list(NFDC)%alive = .true.
				call AssignFDCsites(FDC_TAG,NFDC,nassigned,ok)
				if (.not.ok) then
					return
				endif
				FDC_list(NFDC)%nsites = nassigned
				nBlist = nBlist - nassigned
				checked(dx,dy,dz) = .true.
				write(logmsg,'(a,10i4)') 'FDC site: ilim: ',NFDC,ilim,dx,dy,dz,NFDC,FDC_list(NFDC)%site,nassigned
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
			separation_limit(FDC_TAG,FDC_TAG) = FDC_SEPARATION*(0.8 + 0.3*R)
			separation_limit(FDC_TAG,MRC_TAG) = FDC_MRC_SEPARATION*(0.8 + 0.3*R)
			if (FDCSiteAllowed(FDC_TAG,site,separation_limit)) then
				NFDC = NFDC + 1
				FDC_list(NFDC)%ID = NFDC
				FDC_list(NFDC)%site = site
				FDC_list(NFDC)%alive = .true.
				call AssignFDCsites(FDC_TAG,NFDC,nassigned,ok)
				if (.not.ok) then
					return
				endif
				FDC_list(NFDC)%nsites = nassigned
				nBlist = nBlist - nassigned
				checked(dx,dy,dz) = .true.
!				write(logmsg,'(a,9i4)') 'FDC site: ilim: ',ilim,dx,dy,dz,NFDC,FDC_list(NFDC)%site,nassigned
!				call logger(logmsg)
				if (NFDC == NFDCrequired) then
					done = .true.
					exit
				endif
			endif				
		enddo
		if (done) exit
	enddo
endif

call AssignCellBdrySites(FDC_TAG)

write(logmsg,'(a,2i6)') 'Number of FDCs, B cells: ',NFDC,nBlist
call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
subroutine placeMRCs(ok)
logical :: ok
integer :: NMRCrequired, site(3), nassigned, dx, dy, dz, x, y, z, k, kpar=0
integer, parameter :: MAX_LIM = 30
real :: R, separation_limit(2,2)
real :: MRC_SEPARATION, FDC_MRC_SEPARATION
real :: MRC_LAYER_LIMIT, MRC_YFRAC_MIN
type(vector3_type) :: eradius
integer, parameter :: method = 1	
integer, parameter :: kmax = 100000

if (.not.use_MRCs) then
	NMRC = 0
	ok = .true.
	return
endif
call logger('PlaceMRCs')
DCoffset(:,1) = 0
DCoffset(:,2) = (/-1,0,0/)
DCoffset(:,3) = (/1,0,0/)
DCoffset(:,4) = (/0,-1,0/)
DCoffset(:,5) = (/0,1,0/)
DCoffset(:,6) = (/0,0,-1/)
DCoffset(:,7) = (/0,0,1/)
NMRCrequired = BASE_NMRC
if (method == 1) then
	! The MRCs are located in a layer adjacent to the capsule boundary of the blob,
	! defined by x^2/a^2 + y^2/b^2 + z^2/c^2 > d^2, where d < 1,
	! and y > e
	MRC_YFRAC_MIN = 0.4		! e
	MRC_LAYER_LIMIT = 0.6	! d
	eradius%x = MRC_LAYER_LIMIT*Radius%x
	eradius%y = MRC_LAYER_LIMIT*Radius%y
	eradius%z = MRC_LAYER_LIMIT*Radius%z
    MRC_SEPARATION = 3.0
    FDC_MRC_SEPARATION = 2.0
endif
NMRC = 0
if (method == 1) then
	k = 0
	do
		k = k+1
		if (k > kmax) then
			write(logmsg,*) 'Error: PlaceMRCs: max iterations exceeded: NMRC: ',NMRC
			call logger(logmsg)
			ok = .false.
			return
		endif
		R = par_uni(kpar)
		dx = (2*R-1)*Radius%x
		R = par_uni(kpar)
		dy = MRC_YFRAC_MIN*Radius%y + R*(1 - MRC_YFRAC_MIN)*Radius%y
		R = par_uni(kpar)
		dz = (2*R-1)*Radius%z
		site = Centre + (/dx,dy,dz/)
		if (.not.InsideEllipsoid(site,Centre,Radius)) cycle
		if (InsideEllipsoid(site,Centre,eradius)) cycle
		R = par_uni(kpar)
		separation_limit(MRC_TAG,MRC_TAG) = MRC_SEPARATION*(0.8 + 0.3*R)
		separation_limit(MRC_TAG,FDC_TAG) = FDC_MRC_SEPARATION*(0.8 + 0.3*R)
		if (FDCSiteAllowed(MRC_TAG,site,separation_limit)) then
			NMRC = NMRC + 1
			MRC_list(NMRC)%ID = NMRC
			MRC_list(NMRC)%site = site
			MRC_list(NMRC)%alive = .true.
			call AssignFDCsites(MRC_TAG,NMRC,nassigned,ok)
			if (.not.ok) then
				return
			endif
			MRC_list(NMRC)%nsites = nassigned
			nBlist = nBlist - nassigned
			if (NMRC == NMRCrequired) exit
		endif
	enddo
endif

call AssignCellBdrySites(MRC_TAG)

write(logmsg,'(a,3i6)') 'Number of MRCs, B cells, iterations: ',NMRC,nBlist,k
call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Check that a proposed FDC/MRC site is far enough from all existing FDCs and MRCs
! rlimit(tag1,tag2) = limit separation for tag1 and tag2 cells,
! where tag1,2 is FDC_TAG or MRC_TAG
!--------------------------------------------------------------------------------
logical function FDCSiteAllowed(tag,site0,rlimit)
integer :: tag, site0(3)
real :: rlimit(2,2)
integer :: ifdc, imrc, site(3), k
real :: r(3), d

do ifdc = 1,NFDC
	r = FDC_list(ifdc)%site - site0
	d = norm(r)
	if (d < rlimit(tag,FDC_TAG)) then
		FDCSiteAllowed = .false.
		return
	endif
enddo
do imrc = 1,NMRC
	r = MRC_list(imrc)%site - site0
	d = norm(r)
	if (d < rlimit(tag,MRC_TAG)) then
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
! Determine the FDC and MRC boundary sites that have
! specified chemokine concentration in SolveSteadystate_A(),
! or sites that receive FDC chemokine secretion at specified rates in
! SolveSteadystate_B().
! USE_CELL_SITES = true => actual FDC/MRC sites have specified conc.
!                = false => sites adjacent to an FDC/MRC are flagged,
!                  occupancy(:,:,:)%FDC_nbdry records the number of adjacent FDCs
!--------------------------------------------------------------------------------
subroutine AssignCellBdrySites(tag)
integer :: tag
integer :: i, j, ic, site(3), fsite(3), bsite(3), dx, dy, dz, NC
logical :: bdry(-2:2,-2:2,-2:2)

if (tag == FDC_TAG) then
	NC = NFDC
	occupancy(:,:,:)%FDC_nbdry = 0
elseif (tag == MRC_TAG) then
	NC = NMRC
	occupancy(:,:,:)%MRC_nbdry = 0
endif
if (USE_CELL_SITES) then
	do ic = 1,NC
		bdry = .false.
		do i = 2,7
			if (tag == FDC_TAG) then
				fsite = FDC_list(ic)%site + DCoffset(:,i)
				occupancy(fsite(1),fsite(2),fsite(3))%FDC_nbdry = 1
			elseif (tag == MRC_TAG) then
				fsite = MRC_list(ic)%site + DCoffset(:,i)
				occupancy(fsite(1),fsite(2),fsite(3))%MRC_nbdry = 1
			endif
		enddo
	enddo
else
	do ic = 1,NC
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
						if (tag == FDC_TAG) then
							site = FDC_list(ic)%site + (/dx,dy,dz/)
							occupancy(site(1),site(2),site(3))%FDC_nbdry = occupancy(site(1),site(2),site(3))%FDC_nbdry + 1
						elseif (tag == MRC_TAG) then
							site = MRC_list(ic)%site + (/dx,dy,dz/)
							occupancy(site(1),site(2),site(3))%MRC_nbdry = occupancy(site(1),site(2),site(3))%MRC_nbdry + 1
						endif
					endif
				enddo
			enddo
		enddo
	enddo
endif
end subroutine

!-----------------------------------------------------------------------------------------
! All associations of sites with DC are recomputed.
! That is for every site (x,y,z), create the list of DC that are near this site:
! occupancy(x,y,z)%DC(1:3), where occupancy(x,y,z)%DC(0) = number of nearby DC (if >= 0)
! To avoid having to explicitly select the closest DC (when there are more than DCDIM-1 near
! a site), the order of scanning the DC_list is randomized.
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
    if (.not.DC_list(idc)%alive) cycle
	if (DC_list(idc)%nsites < NDCsites) then
		call AssignFDCsites(0,idc,nassigned,ok)
		if (.not.ok) return
		DC_list(idc)%nsites = nassigned
	endif
    NDCalive = NDCalive + 1
    site = DC_list(idc)%site
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
    if (DC_list(idc)%alive) exit
enddo
dir = random_int(1,6,kpar)
call moveDC(idc,dir)
end subroutine

!-----------------------------------------------------------------------------------------
! A DC is allowed to move only if it has grown to its full extent,
! i.e. if DC_list()%nsites = NDCsites
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

if (DC_list(idc)%nsites /= NDCsites) return
step = neumann(:,dir)
site0 = DC_list(idc)%site
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
                Bcell_list(kcell)%site = site1
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
                Bcell_list(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    endif
enddo
DC_list(idc)%site = site0 + step
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
! tag = FDC_TAG, MRC_TAG
! ic = ifdc, imrc
!-----------------------------------------------------------------------------------------
subroutine AssignFDCsites(tag,ic,nassigned,ok)
integer :: tag, ic, nassigned
logical :: ok
integer :: ifdc, imrc, k, site(3), site1(3), indx(2)

ok = .false.
if (tag == FDC_TAG) then
	ifdc = ic
	site = FDC_list(ifdc)%site
	occupancy(site(1),site(2),site(3))%indx = -ifdc				! This site holds an FDC (centre)
elseif (tag == MRC_TAG) then
	imrc = ic
	site = MRC_list(imrc)%site
	occupancy(site(1),site(2),site(3))%indx = -(1000 + imrc)	! This site holds an MRC (centre)
endif
nassigned = 1
do k = 2,NDCsites
    site1 = site + DCoffset(:,k)
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if ((tag == FDC_TAG .and. indx(1) == -ifdc) .or. (tag == MRC_TAG .and. indx(1) == -imrc)) then
		nassigned = nassigned + 1
		cycle
	endif
    ! Reduce available blob sites only if the site is within the blob.
    ! A DC site should never fall outside, I think.  Check this.
    if (.not.inside_xyz(site1)) then
		call logger("Error: AssignFDCsites: site outside grid")
		return
	endif
    if (indx(1) == OUTSIDE_TAG) then
		call logger('AssignFDCsites: cell hits boundary')
		return
	endif
	if (indx(1) /= 0) cycle
	if (indx(2) /= 0) cycle
	nassigned = nassigned + 1
    occupancy(site1(1),site1(2),site1(3))%indx = -ic     ! This site holds a FDC or MRC (soma)
enddo
ok = .true.
end subroutine

end module
