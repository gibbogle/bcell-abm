! omp_diffuse.f90
!----------------------------------------------------
! This module handles diffusion of cytokines.
!----------------------------------------------------

module diffuse

use behaviour

implicit none

real, allocatable :: cytp(:,:,:,:)
logical, allocatable :: inblob(:,:,:)
integer, allocatable :: xminmax(:,:,:)
integer :: yminmax(2),zminmax(2)
integer, allocatable :: sitelist(:,:,:)
integer, allocatable :: neighbours(:,:,:)
!integer, parameter :: MAXX = 200000     ! = max sites in a slice = 1.333*pi*R^3/(2*Mnodes)

!integer :: sitelist(MAXX,3,8)
!integer :: neighbours(0:6,MAXX,8)
!logical :: inblob(NX,NY,NZ)
!integer :: xminmax(NY,NZ,2)
!real :: cytp(NX,NY,NZ,N_CYT)

contains

!-----------------------------------------------------------------------------------------
! This is just for testing to see how useful the idea is
!-----------------------------------------------------------------------------------------
subroutine setup_minmax
integer :: x,y,z,xmin,xmax,ymin,ymax,zmin,zmax

inblob = .false.
ymin = 99999
ymax = -99999
zmin = 99999
zmax = -99999
do z = 1,NZ
    do y = 1,NY
        xmin = 99999
        xmax = -99999
        do x = 1,NX
            if (occupancy(x,y,z)%indx(1) >= 0) then
                xmin = min(x,xmin)
                xmax = max(x,xmax)
                ymin = min(y,ymin)
                ymax = max(y,ymax)
                zmin = min(z,zmin)
                zmax = max(z,zmax)
                inblob(x,y,z) = .true.
            endif
        enddo
        xminmax(y,z,:) = (/xmin,xmax/)
    enddo
enddo
yminmax = (/ymin,ymax/)
zminmax = (/zmin,zmax/)
end subroutine

!-----------------------------------------------------------------------------------------
! Parallelized
! This is slow.  It should help to precompute array(s) conveying which neighbour sites
! need to be considered for each (x,y,z) DONE!
!-----------------------------------------------------------------------------------------
subroutine diffuser
integer :: idstep,kIL2,sweep,kpar
integer :: x,y,z,xx,yy,zz,dir,kt,z_lo,z_hi,ndsites(8),slice,nslices,nb(0:6),xmin,xmax,icyt
real(8) :: asum
logical, save :: first = .true.
!integer, allocatable :: sitelist(:,:,:), neighbours(:,:,:)
!real, allocatable :: cytp(:,:,:,:)

!allocate(sitelist(MAXX,3,2*Mnodes))
!allocate(neighbours(0:6,MAXX,2*Mnodes))
!allocate(cytp(NX,NY,NZ,Ncytokines))

if (Ncytokines > 1) then
    write(*,*) 'ERROR: diffuser: par_diffuser expects Ncytokines = 1, not: ',Ncytokines
    stop
endif

if (first) then
    call setup_minmax
    first = .false.
endif

if (IL2_DECAY_RATE /= 0) then
    kIL2 = cyt_seq(IL2_TAG)
    cyt(:,:,:,kIL2) = cyt(:,:,:,kIL2)*(1 - IL2_DECAY_RATE*DELTA_T)  ! decay of IL2
endif

if (Mnodes == 1) then
    nslices = 1
else
    nslices = 2*Mnodes
endif

do slice = 1,nslices
    if (Mnodes == 1) then
        z_lo = 1
        z_hi = NZ
    else
        z_lo = zoffset(slice-1) + 1
        z_hi = zoffset(slice)
    endif
    z_lo = max(z_lo,zminmax(1))
    z_hi = min(z_hi,zminmax(2))
    if (z_lo > z_hi) cycle
    ndsites(slice) = 0
    do z = z_lo,z_hi
        do y = yminmax(1),yminmax(2)
!        do y = 1,NY
!            do x = 1,NX
            xmin = xminmax(y,z,1)
            xmax = xminmax(y,z,2)
            if (xmin > xmax) cycle
            do x = xmin,xmax
!                if (occupancy(x,y,z)%indx(1) < 0) cycle
                if (.not.inblob(x,y,z)) cycle
                nb(0:6) = 0
                kt = 0
                do dir = 1,6
                    xx = x + neumann(1,dir)
                    yy = y + neumann(2,dir)
                    zz = z + neumann(3,dir)
                    if (xx < 1 .or. xx > NX) cycle
                    if (yy < 1 .or. yy > NY) cycle
                    if (zz < 1 .or. zz > NZ) cycle
!                    if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
                    if (.not.inblob(xx,yy,zz)) cycle
                    kt = kt+1
                    nb(kt) = dir
                enddo
                if (kt == 0) cycle
                nb(0) = kt
                ndsites(slice) = ndsites(slice) + 1       ! this is the number of sites in the slice
                sitelist(ndsites(slice),:,slice) = (/x,y,z/)
                neighbours(:,ndsites(slice),slice) = nb
            enddo
        enddo
    enddo
enddo

do idstep = 1,NDIFFSTEPS
    if (Mnodes > 1) then
        do sweep = 0,1
            !$omp parallel do private(slice)
            do kpar = 0,Mnodes-1
                slice = sweep+2*kpar+1
                if (mod(idstep,2) == 1) then
                    call par_diffuser(cyt,cytp,ndsites(slice),sitelist(:,:,slice),neighbours(:,:,slice))
                else
                    call par_diffuser(cytp,cyt,ndsites(slice),sitelist(:,:,slice),neighbours(:,:,slice))
                endif
            enddo
            !$omp end parallel do
        enddo
    else
        slice = 1
        if (mod(idstep,2) == 1) then
            call par_diffuser(cyt,cytp,ndsites(slice),sitelist(:,:,slice),neighbours(:,:,slice))
        else
            call par_diffuser(cytp,cyt,ndsites(slice),sitelist(:,:,slice),neighbours(:,:,slice))
        endif
    endif
enddo
if (mod(NDIFFSTEPS,2) == 1) then
    cyt = cytp
endif
do icyt = 1,Ncytokines
!    csum = 0
!    nsum = 0
!    cmin = 99999
!    cmax = -99999
!    do slice = 1,nslices
!        do isite = 1,ndsites(slice)
!            x = sitelist(isite,1,slice)
!            y = sitelist(isite,2,slice)
!            z = sitelist(isite,3,slice)
!            nsum = nsum + 1
!            csum = csum + cyt(x,y,z,icyt)
!            if (cyt(x,y,z,icyt) < cmin) then
!                cmin = cyt(x,y,z,icyt)
!                ix = x
!                iy = y
!                iz = z
!            endif
!            cmax = max(cmax,cyt(x,y,z,icyt))
!        enddo
!    enddo
    asum = sum(cyt(:,:,:,icyt))
    cyt_mean(icyt) = asum/NBcells
enddo
!write(*,*) 'nsum,NBcells: ',nsum,NBcells
!write(*,*) 'asum,csum: ',asum,csum
!write(*,*) 'cmin,cmax: ',cmin,cmax,ix,iy,iz
!stop
!deallocate(sitelist)
!deallocate(neighbours)
!deallocate(cytp)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_diffuser(c1,c2,ndsites,sitelist,neighbours)
integer :: ndsites,sitelist(:,:),neighbours(0:,:)
real :: c1(:,:,:,:), c2(:,:,:,:)
integer :: x,y,z,xx,yy,zz,kt,dir,k,isite,nb(0:6)
!real csum(MAX_CYT), local_delta_diff(MAX_CYT)
!real :: constit(MAX_CYT)
real :: csum,local_delta_diff,constit   ! test for a single cytokine case
integer, parameter :: local_neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))
logical :: flag

!local_delta_diff(1:Ncytokines) = delta_diff(1:Ncytokines)
local_delta_diff = delta_diff(1)

!constit(1:Ncytokines) = cyt_constit(1:Ncytokines)*(DELTA_T/NDIFFSTEPS)*M_pM*L_um3/Navo      ! constitutive production of cytokines
constit = cyt_constit(1)*(DELTA_T/NDIFFSTEPS)*M_pM*L_um3/Navo      ! constitutive production of cytokines

!nhits = 0
!!do z = z_lo,z_hi
!!    do y = 1,NY
!!        do x = 1,NX
do isite = 1,ndsites
    flag = .false.
    x = sitelist(isite,1)
    y = sitelist(isite,2)
    z = sitelist(isite,3)
    nb = neighbours(0:6,isite)
    kt = nb(0)
    csum = constit
    do k = 1,min(kt,6)
        dir = nb(k)
        xx = x + local_neumann(1,dir)
        yy = y + local_neumann(2,dir)
        zz = z + local_neumann(3,dir)
!        csum(1:Ncytokines) = csum(1:Ncytokines) + delta_diff(1:Ncytokines)*c1(xx,yy,zz,1:Ncytokines)
        csum = csum + local_delta_diff*c1(xx,yy,zz,1)

!        kt = 0
!        csum(1:Ncytokines) = constit(1:Ncytokines)
!        do dir = 1,6
!            xx = x + local_neumann(1,dir)
!            yy = y + local_neumann(2,dir)
!            zz = z + local_neumann(3,dir)
!            if (xx < 1 .or. xx > NX) cycle
!            if (yy < 1 .or. yy > NY) cycle
!            if (zz < 1 .or. zz > NZ) cycle
!            if (occupancy(xx,yy,zz)%indx(1) < 0) cycle
!            csum(1:Ncytokines) = csum(1:Ncytokines) + local_delta_diff(1:Ncytokines)*c1(xx,yy,zz,1:Ncytokines)
!            kt = kt + 1

    enddo
!    if (kt /= 0) nhits = nhits + 1
!    c2(x,y,z,1:Ncytokines) = csum(1:Ncytokines) + (1 - kt*delta_diff(1:Ncytokines))*c1(x,y,z,1:Ncytokines)
    c2(x,y,z,1) = csum + (1 - kt*local_delta_diff)*c1(x,y,z,1)
enddo
!write(*,*) 'par_diffuser: ',kpar,sweep,nhits
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_tester(aa,bb,sweep,kpar)
integer :: sweep, kpar
real, allocatable :: aa(:,:,:),bb(:,:,:)
integer :: x,y,z,xx,yy,zz,kt,dir,z_lo,z_hi
real :: c(6,1),csum(1)
logical :: included(6)

allocate(aa(NX,NY,NZ))
allocate(bb(NX,NY,NZ))

if (Mnodes == 1) then
    z_lo = 1
    z_hi = NX
else
    z_lo = zoffset(sweep+2*kpar) + 1
    z_hi = zoffset(sweep+2*kpar+1)
endif
c = 0
do z = 1,z_lo,z_hi
    do y = 1,NY
        do x = 1,NX
            kt = 0
            included = .false.
            do dir = 1,6
                xx = x + neumann(1,dir)
                yy = y + neumann(2,dir)
                zz = z + neumann(3,dir)
                if (xx < 1 .or. xx > NX) cycle
                if (yy < 1 .or. yy > NY) cycle
                if (zz < 1 .or. zz > NZ) cycle
                included(dir) = .true.
                c(dir,1) = aa(xx,yy,zz)
                kt = kt + 1
            enddo
            csum(1) = (1 - kt*delta_diff(1))*aa(x,y,z)
            do dir = 1,6
                if (included(dir)) then
                    csum(1) = csum(1) + delta_diff(1)*c(dir,1)
                endif
            enddo
            bb(x,y,z) = csum(1)
        enddo
    enddo
enddo

deallocate(aa)
deallocate(bb)

end subroutine

!-----------------------------------------------------------------------------------------
! If "mean field" concentrations are used for cytokines (i.e. not use_diffusion)
! the net number of molecules released into the blob in one time step (in updater()) is
! accumulated in dcyt_mols(:)
! The overall total count of molecules in cyt_mols(:) is updated every time step and
! broadcast to the nodes.
!-----------------------------------------------------------------------------------------
subroutine molsynch
integer :: k, kIL2  !, source, dest, ierr, status(MPI_STATUS_SIZE)
real :: mols_pM
!real, allocatable :: dmols(:)

!allocate(dmols(Ncytokines))
!if (me == 0) then
    kIL2 = cyt_seq(IL2_TAG)
    cyt_mols(kIL2) = cyt_mols(kIL2)*(1 - IL2_DECAY_RATE*DELTA_T)                ! decay of IL2
    mols_pM = L_um3*M_pM/(NBcells*Vc*Navo)
!    write(*,'(a,4f8.1)') 'molsynch: ',cyt_mols(kIL2),dcyt_mols(kIL2), &
!        cyt_constit(kIL2)*DELTA_T*NBcells*Vc,mols_pM*cyt_mols(kIL2)

    cyt_mols = cyt_mols + dcyt_mols + cyt_constit*DELTA_T*NBcells*Vc  ! constitutive production of cytokines
    dcyt_mols = 0
!    if (Mnodes > 1) then
!        do source = 1,Mnodes-1
!            call MPI_RECV(dmols, Ncytokines, MPI_REAL, source, MOL1_TAG, MPI_COMM_WORLD, status, ierr)
!            cyt_mols = cyt_mols + dmols
!        enddo
!        do dest = 1,Mnodes-1
!            call MPI_SEND(cyt_mols, Ncytokines, MPI_REAL, dest, MOL2_TAG, MPI_COMM_WORLD, ierr)
!        enddo
!    endif
!elseif (Mnodes > 1) then
!    dest = 0
!    call MPI_SEND(dcyt_mols, Ncytokines, MPI_REAL, dest, MOL1_TAG, MPI_COMM_WORLD, ierr)
!    dcyt_mols = 0
!    source = 0
!    call MPI_RECV(cyt_mols, Ncytokines, MPI_REAL, source, MOL2_TAG, MPI_COMM_WORLD, status, ierr)
!endif
!deallocate(dmols)
do k = 1,Ncytokines
    if (cyt_mols(k) < 0) then
        write(*,*) 'molsynch: cyt_mols < 0: ',k,cyt_mols(k)
        stop
    endif
enddo
end subroutine


end module
