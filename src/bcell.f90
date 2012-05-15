!-----------------------------------------------------------------------------------------
! Simplifying assumptions: 
! (1) The follicle is assumed always to be approximately bounded by a spheroid, 
!     the surface created by rotating an ellipse about the X axis.  The ellipse
!     is defined by the major axis radius (in the X direction) and the eccentricity e.
!     If the major axis radius = a, minor = b, b/a = (1-e^2)^1/2
! (2) The partitioning of the blob is always created by slice boundaries that
!     are perpendicular to the X axis, i.e. parallel to the YZ plane.
! (3) At any x value in every slice the available sites form a convex set in 2D.
!     I.e. for any y value there are either no available sites, or all the
!     sites in the range y1 <= y <= y2 are available.  Similarly for any z.
!     This means that the blob shape is fully described by 2 sets of spans:
!     span(x,y) x=1,NX, y=1,NY => available sites are for span%lo <= z <= span%hi
!     span(z,x) x=1,NX, z=1,NZ => available sites are for span%lo <= y <= span%hi
!
! Implementing a list of cells:
! In this version the cells in the domain are stored in a list, while the
! occupancy array holds the indices of cells in the list.  When a cell
! leaves the domain (i.e. crosses domain boundary, leaves the paracortex,
! or dies) a gap is created in the list.  The locations of such gaps are stored
! in the gaplist, the total number of such gaps is ngaps.  A cell entering the
! domain is allocated an index from the tail of this list, if ngaps > 0, or
! else it is added to the end of the cell list.
!
! Cell trafficking:
! B cells enter the follicle at a rate that is determined by the vascularity,
! which in turn is controlled by the inflammation signal.  Cells actually 
! transmigrate from the blood to the LN paracortex via HEVs, and then migrate
! to the lower follicle boundary, presumably using chemotactic cues.  For this
! model we assume that B cells arrive at the T zone-follicle interface at a
! computed rate.  The entry locations will be randomly (uniformly) distributed
! on a portion of the lower surface of the ellipsoid.  If the minor axis radius
! of the ellipse is Rb (= Ra/ELLIPSE_RATIO), the entry locations are given by
! surface sites with y < -ENTRY_ALPHA*Rb, where ENTRY_ALPHA is a specified parameter,
! e.g. ENTRY_ALPHA = 0.5.  The intersection of the plane y = -ENTRY_ALPHA*Rb with
! the ellipsoid surface defines an ellipse with major axis radius given by:
!  Ra*(1-ENTRY_ALPHA^2)^1/2
! which equals Ra/2 when ENTRY_ALPHA = 0.75^1/2 = 0.866
! B cell egress also occurs on lower surface, in a portion parametrised by EXIT_ALPHA.
!-----------------------------------------------------------------------------------------

module main_mod
use global
use behaviour
use ode_diffuse
use FDC
use fields
use winsock
!use aviewer 

IMPLICIT NONE

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
if (Mnodes == 1) return
!!DEC$ IF ( DEFINED (_OPENMP) .OR. DEFINED (IBM)) 
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine array_initialisation(ok)
logical :: ok
integer :: x,y,z,k
integer :: MAXX, z1, z2
integer :: cog_size !, noncog_size
real :: d, rr(3)
type(cog_type) :: cog

ok = .false.
call logger("call rng_initialisation")
call rng_initialisation
call logger("did rng_initialisation")
!call check_rng

cog_size = (sizeof(cog) + 2)/4
nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
ngaps = 0
max_ngaps = 5*NY*NZ
ID_offset = BIG_INT/Mnodes
nlist = 0
MAX_COG = 0.5*NX*NY*NZ

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
allocate(xoffset(0:2*Mnodes))
allocate(xdomain(NX))
x0 = (NX + 1.0)/2.        ! global value
y0 = (NY + 1.0)/2.
z0 = (NZ + 1.0)/2.
aRadius = BLOB_RADIUS     ! starting value
if (2*aRadius > NX-2) then
    write(logmsg,'(a,i6,f7.1)') 'ERROR: NX too small for RADIUS: ',NX,aRadius
	call logger(logmsg)
	logmsg = 'Fix: Reduce BLOB_RADIUS or increase NX'
	call logger(logmsg)
    return
endif
bRadius = aRadius/ELLIPSE_RATIO

max_nlist = 1.5*NX*NY*NZ

allocate(DClist(MAX_DC))
allocate(FDClist(MAX_FDC))
allocate(occupancy(NX,NY,NZ))
allocate(cellist(max_nlist))
allocate(gaplist(max_ngaps))
allocate(nz_sites(NZ))
allocate(nz_totsites(NZ))
allocate(nz_cells(NZ))
allocate(nz_excess(NZ))
allocate(cognate_list(MAX_COG))

do k = 1,max_nlist
	nullify(cellist(k)%cptr)
enddo


call make_reldir

nz_excess = 0
Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
ncogseed = 0
lastcogID = 0
lastID = 0
k_nonrandom = 0
lastNBcells = 0
nadd_sites = 0
lastbalancetime = 0
localres%dN_EffCogTC = 0
localres%dN_EffCogTCGen = 0
localres%N_EffCogTC = 0
localres%N_EffCogTCGen = 0
totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%N_EffCogTC = 0
totalres%N_EffCogTCGen = 0
totalres%N_dead = 0
totalres%dN_dead = 0

if (evaluate_residence_time) then
    allocate(Tres_dist(int(days*24)))
    Tres_dist = 0
endif

!if (use_cytokines) then
!    allocate(cytp(NX,NY,NZ,N_CYT))
!endif
!if (use_diffusion) then
!    if (.not.use_cytokines) then
!        write(logmsg,*) 'Cannot use_diffusion without use_cytokines'
!	    call logger(logmsg)
!        stop
!    endif
!    allocate(xminmax(NY,NZ,2))
!    allocate(inblob(NX,NY,NZ))
!    MAXX = 1.5*PI*(NX/2)**3/(2*Mnodes)
!    allocate(sitelist(MAXX,3,8))
!    allocate(neighbours(0:6,MAXX,8))
!endif

ok = .true.

end subroutine

!--------------------------------------------------------------------------------
! Generates the arrays wz(), zoffset() and zdomain().
! The domains (slices) are numbered 0,...,2*Mnodes-1
! wz(k) = width of the slice for kth domain
! zoffset(k) = offset of kth domain occupancy array in the occupancy array.
! zdomain(x) = domain that global z lies in.
! The kth domain (slice) extends from z = zoffset(k)+1 to z = zoffset(k+1)
! The idea is to set the domain boundaries such that each domain has roughly the
! same number of available sites.
! This is the initial split, which will continue to be OK if:
! not using a blob, or Mnodes <= 2
!--------------------------------------------------------------------------------
subroutine make_split
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real :: dNT, diff1, diff2
!integer, save :: lastNBcells = 0
!type(cell_type), pointer :: cell
logical :: show = .false.

!write(*,*) 'make_split: Mnodes: ',Mnodes,use_blob
if (Mnodes == 1) then
    Mslices = 1
    zdomain = 0
    return
endif
Mslices = 2*Mnodes
allocate(wz(0:Mslices))
allocate(ztotal(0:Mslices))
allocate(scount(NX))
if (use_blob) then
    dNT = abs(NBcells - lastNBcells)/real(lastNBcells+1)
    if (dNT < 0.03) then
!       write(*,*) 'debugging make_split: ',NBcells,lastNBcells,dNT
        return
    endif
    lastNBcells = NBcells
    if (show) write(*,*) 'make_split: dNT: ',dNT
    nsum = 0
    do z = 1,NZ
        k = 0
        do y = 1,NY
            do x = 1,NX
                if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
                    k = k + 1
                endif
            enddo
        enddo
        scount(z) = k
        nsum = nsum + scount(z)
    enddo
    Ntot = nsum
    N = Ntot/Mslices
    nsum = 0
    last = 0
    k = 0
    do z = 1,NZ
        nsum = nsum + scount(z)
        if (nsum >= (k+1)*N) then
            diff1 = nsum - (k+1)*N
            diff2 = diff1 - scount(z)
            if (abs(diff1) < abs(diff2)) then
                wz(k) = z - last
                last = z
            else
                wz(k) = z - last - 1
                last = z - 1
            endif
            k = k+1
            if (k == Mslices-1) exit
        endif
    enddo
    wz(Mslices-1) = NZ - last
    if (show) then
        write(*,*) 'Ntot, N: ',Ntot,N
        write(*,'(10i6)') scount
    endif
else
    wz = NZ/Mslices
    wsum = 0
    do k = 0,Mslices-1
        wsum = wsum + wz(k)
    enddo
    do k = 0,Mslices-1
        if (wsum < NZ) then
            wz(k) = wz(k) + 1
            wsum = wsum + 1
        endif
    enddo
endif
zoffset(0) = 0
do k = 1,Mslices-1
    zoffset(k) = zoffset(k-1) + wz(k-1)
enddo
zoffset(Mslices) = NZ
z = 0
do kdomain = 0,Mslices-1
    do k = 1,wz(kdomain)
        z = z+1
        zdomain(z) = kdomain      ! = kpar with two sweeps
    enddo
enddo
if (show) then
    write(*,*) 'zoffset: ',zoffset
    write(*,*) 'wz:      ',wz
    write(*,*) 'zdomain: '
    write(*,'(10i4)') zdomain
endif
ztotal = 0
do k = 0,2*Mnodes-1
    do z = zoffset(k)+1,zoffset(k+1)
        ztotal(k) = ztotal(k) + scount(z)
    enddo
    if (show) write(*,*) k,ztotal(k)
enddo
deallocate(wz)
deallocate(ztotal)
deallocate(scount)
end subroutine

!--------------------------------------------------------------------------------
! Makes an approximate count of the number of sites of the spherical blob that
! are in the xth slice.  Uses the area of the slice.
! The blob centre is at (x0,y0,z0), and the blob radius is R = aRadius
!--------------------------------------------------------------------------------
integer function slice_count(x)
integer :: x
real :: r2

r2 = aRadius**2 - (x-x0)**2
if (r2 < 0) then
    slice_count = 0
else
    slice_count = PI*r2
endif
end function


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine motility_calibration
integer :: ic,id,k,ds(3), imin, isub
integer :: NBcells, nvar0, nvar, ntime2
integer :: ns = 10000
integer :: npaths = 20
integer :: npos = 301
integer :: nbeta = 30, nrho = 50
integer :: ibeta, irho, kpath
real :: dt, time1 = 5, time2 = 30   !time2 = 15   ! minutes (used for ICB paper results)
real :: tagrad = 10
real :: dbeta, drho
real :: betamin = 0.25, betamax = 0.9		! 0.25 - 0.90 for Model_N, 0.15 - 0.80 for Model_M
real :: rhomin = 0.20, rhomax = 0.85			! 0.20 - 0.85 for Model_N, 0.20 - 0.85 for Model_M
real :: Cm,speed,ssum,d
integer, allocatable :: tagid(:), tagseq(:), tagsite(:,:,:), pathcell(:)
integer, allocatable :: prevsite(:,:)   ! for mean speed computation
real, allocatable :: Cm_array(:,:), S_array(:,:)
type(cell_type) :: cell
logical :: ok

write(*,*) 'motility_calibration'

use_chemotaxis = .false.
chemo_K_exit = 0
NBcells = NX*NY*NZ
if (motility_save_paths) then
    allocate(pathcell(npaths))
endif

write(*,*) 'rho,beta: ',rho,beta

ntime2 = time2
nvar0 = time1
nvar = time2    ! number of minutes to simulate
dt = DELTA_T*nsteps_per_min
TagRadius = tagrad

if (motility_param_range) then
	drho = (rhomax-rhomin)/(nrho-1)
	dbeta = (betamax-betamin)/(nbeta-1)
elseif (n_multiple_runs > 1) then
	nbeta = n_multiple_runs
	nrho = 1
	betamin = beta
	rhomin = rho
	drho = 0
	dbeta = 0
else
	nbeta = 1
	nrho = 1
	betamin = beta
	rhomin = rho
	drho = 0
	dbeta = 0
endif

allocate(Cm_array(nrho,nbeta))
allocate(S_array(nrho,nbeta))

write(*,*) 'nbeta,nrho: ',nbeta,nrho
do ibeta = 1,nbeta
    do irho = 1,nrho
	    rho = rhomin + (irho-1)*drho
	    beta = betamin + (ibeta-1)*dbeta
	    write(*,'(a,2i3,2f6.2)') ' beta, rho: ',ibeta,irho,BETA,RHO
	    call compute_dirprobs
	    write(*,*) 'dirprob: '
	    write(*,'(10f7.3)') dirprob(0:nreldir)
	    call placeCells(ok)
	    if (.not.ok) stop
        call make_split
        if (nlist > 0) then
	        write(*,*) 'make tag list: NBcells,nlist,ntagged: ',NBcells,nlist,ntagged

            allocate(tagseq(NBcells))
            allocate(tagid(ntagged))
            allocate(tagsite(3,ntagged,0:nvar))
            tagseq = 0
            k = 0
	        kpath = 0
            do ic = 1,nlist
                if (cellist(ic)%ctype == TAGGED_CELL) then
                    id = cellist(ic)%ID
                    k = k+1
                    tagid(k) = id
                    tagseq(id) = k
                    tagsite(:,k,0) = cellist(ic)%site
					if (motility_save_paths) then
						if (kpath < npaths) then
							kpath = kpath + 1
							pathcell(kpath) = ic
						endif
					endif
                endif
            enddo
        endif

	    if (checking > 0) call checker
	    if (checking > 0) call checker
        ns = min(ns,nlist)
        allocate(prevsite(3,ns))
        do ic = 1,ns
            prevsite(:,ic) = cellist(ic)%site
        enddo
        ssum = 0

        !
        ! Now we are ready to run the simulation
        !
        if (motility_save_paths) then
            open(nfpath,file='path.out',status='replace')
            write(nfpath,'(i3,a)') npaths,' paths'
        endif

        write(*,*) 'nvar, nsteps_per_min: ',nvar, nsteps_per_min
        istep = 0
        do imin = 1,nvar
            do isub = 1,nsteps_per_min
                istep = istep + 1
                call mover(ok)
                if (.not.ok) stop
                call squeezer(.false.)
!                write(*,*) 'Cell 1: ',cellist(1)%site
                do ic = 1,ns
                    ds = cellist(ic)%site - prevsite(:,ic)
                    prevsite(:,ic) = cellist(ic)%site
                    d = sqrt(real(ds(1)*ds(1) + ds(2)*ds(2) + ds(3)*ds(3)))
                    ssum = ssum + d*DELTA_X/DELTA_T
                enddo
                if (motility_save_paths) then
                    k = (imin-1)*nsteps_per_min + isub
                    if (k >= nvar0*nsteps_per_min .and. k < nvar0*nsteps_per_min + npos) then
                        write(nfpath,'(160i4)') (cellist(pathcell(kpath))%site(1:2),kpath=1,npaths)
                    endif
                endif
            enddo
            write(*,*) 'speed: ',ssum/(ns*nsteps_per_min*imin)

            do ic = 1,nlist
                cell = cellist(ic)
                if (cell%ctype == TAGGED_CELL) then
                    id = cell%ID
                    k = tagseq(id)
                    tagsite(:,k,imin) = cell%site
!                   if (k < 20) write(*,*) 'site: ',cell%ID,cell%site
                endif
            enddo
        enddo

        call compute_Cm(tagsite,ntagged,nvar0,nvar,dt,Cm)
        speed = ssum/(ns*nvar*nsteps_per_min)
        write(*,'(a,2f8.2)') 'speed, Cm: ',speed,Cm
        if (allocated(tagid))   deallocate(tagid)
        if (allocated(tagseq))  deallocate(tagseq)
        if (allocated(tagsite)) deallocate(tagsite)
	enddo
enddo
if (allocated(pathcell)) deallocate(pathcell)
deallocate(Cm_array)
deallocate(S_array)
if (motility_save_paths) then
    close(nfpath)
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine diffusion_calibration
!!!use ifport
integer :: x, y, z, istep
integer, parameter :: Ntimes = 50   !500
real(DP) :: t1, t2, Cnumeric(Ntimes),Canalytic(Ntimes)
logical :: ok

write(*,*) 'diffusion_calibration: '

x = NX/4
!call analytical_soln(x,Canalytic,Ntimes)
call placeCells(ok)
if (.not.ok) stop
call make_split

!t1 = timef()
call cpu_time(t1)
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            occupancy(x,y,z)%indx = 1
        enddo
    enddo
enddo

!do istep = 1,Ntimes
!    call diffuser
!    write(*,*) istep,cyt(NX/4,NY/2,NZ/2,1)
!    Cnumeric(istep) = cyt(NX/4,NY/2,NZ/2,1)
!enddo
!t2 = timef()
call cpu_time(t2)
write(*,'(a,f10.2)') 'Time: ',t2-t1
write(nfout,*) 'NX, x, NDIFFSTEPS: ',NX,x,NDIFFSTEPS
do istep = 1,Ntimes
    write(nfout,'(i6,f8.2,2f8.4)') istep,istep*DELTA_T,Cnumeric(istep),Canalytic(istep)
enddo

call wrapup
stop

end subroutine


!-----------------------------------------------------------------------------------------
! To test AddBcell()
!-----------------------------------------------------------------------------------------
subroutine add_random_cells(n, ctype, gen, stage, region)
integer :: n, ctype, gen, stage, region
integer :: k, x, y, z, site(3), slots, kcell
integer :: kpar=0
logical :: ok

k = 0
do while (k < n)
    x = random_int(1,NX,kpar)
    y = random_int(1,NY,kpar)
    z = random_int(1,NZ,kpar)
    site = (/x,y,z/)
    slots = getslots(site)
    if (occupancy(x,y,z)%indx(1) >= 0 .and. slots < BOTH) then
        if (dbug) write(*,'(a,7i6)') 'add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call AddBcell(site,ctype,gen,stage,region,kcell,ok)
        if (dbug) write(*,'(a,7i6)') 'after add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call checkslots('add_random_cells: ',site)
        k = k+1
    endif
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Scans the cell list and builds the counts nz_sites() and nz_excess().
! nz_sites(k)  = number of available sites at the slice z = k
! nz_excess(k) = total excess of cells over sites in the paracortex zone with z >= k
! These values are used in jumper() to adjust the jump probabilities.  The probability
! of jumps in the -z direction is increased by an increment that is proportional to
! nz_excess(k)/nz_sites(k)
! This version is to quantify the cell distribution in the blob.
!-----------------------------------------------------------------------------------------
subroutine scanner
integer :: x, y, z, ns, nc, nst, nct, indx(2), k
integer :: yfraction(NY), ynsb(NY), yncb(NY)
real, allocatable, save :: ysum(:)
integer, save :: nh = 0
integer :: excess, nextra, i, idn, imin=0, nz1, nz2
real :: eratio(NZ), nz_sites0(NZ)
real :: df, dfmin
integer :: nsb, ncb, nsbt, ncbt
type (boundary_type), pointer :: bdry

if (nh == 0) then
	allocate(ysum(NY))
	ysum = 0
endif
yfraction = 0
ynsb = 0
yncb = 0
nct = 0
nst = 0
nsbt = 0
ncbt = 0
k = 0
do y = 1,NY
    ns = 0
    nc = 0
	nsb = 0
	ncb = 0
    do z = 1,NZ
        do x = 1,NX
            indx = occupancy(x,y,z)%indx
            if (indx(1) < 0) cycle       ! OUTSIDE_TAG or DC
            ns = ns+1
            if (indx(1) > 0) nc = nc+1
            if (indx(2) > 0) nc = nc+1
            bdry => occupancy(x,y,z)%bdry
            if (associated(bdry)) then
				if (bdry%exit_OK) then
		            nsb = nsb + 1
			        if (indx(1) > 0) ncb = ncb+1
				    if (indx(2) > 0) ncb = ncb+1
				endif
			endif	
        enddo
    enddo
    nst = nst + ns
    nct = nct + nc
    nsbt = nsbt + nsb
    ncbt = ncbt + ncb
    if (ns > 0) then
		k = k+1
		yfraction(k) = (100.*nc)/(2.*ns) + 0.5
		ynsb(k) = nsb
		yncb(k) = ncb
	endif
!    nz_sites(z) = ns
!    nz_cells(z) = nc
!    excess = excess + nc - ns
!    nz_excess(z) = excess
enddo
nh = nh + 1
ysum = ysum + yfraction
write(nfout,'(a,f8.2)') 'Hour: ',istep*DELTA_T/60
write(nfout,'(25i3)') int(ysum(1:k)/nh)
write(nfout,'(25i4)') int(ynsb(1:k/2))
write(nfout,'(25i4)') int(yncb(1:k/2))
write(nfout,'(a,f8.4)') 'Fraction of exit sites occupied: ',real(ncbt)/(2*nsbt)
write(nfout,*) 'Total exit sites, exit cells: ',nsbt,ncbt
NXcells = ncbt
return

nz_excess = 0
excess = 0
nextra = nst - nct
! This is the imbalance between number of sites and number of T cells
! resulting from (a) DC sites and (b) sites to be added (nadd_sites)
! We need to adjust either nz_sites(:) or nz_cells(:) to bring them into balance,
! so that eratio(:) can be computed correctly to generate a drift.
! The complication arises because we want to spread the adjustment over the slices
! in a way that is proportionate to the number of sites in the slice.
! Choose to adjust nz_sites(:) (arbitrarily).

if (nextra /= 0) then
    if (nextra > 0) then
        idn = -1
    else
        idn = 1
    endif
    nz_sites0 = nz_sites    ! This conveys approx the shape of the blob - we want to maintain this
    do k = 1,abs(nextra)    ! we need to remove/add this many sites from/to nz_sites(:)
        dfmin = 1.0e10
        do i = 1,NZ
            if (nz_sites(i) == 0) cycle
            df = abs(real(nz_sites(i) + idn)/nz_sites0(i) - 1)
            if (df < dfmin) then
                dfmin = df
                imin = i
            endif
        enddo
        nz_sites(imin) = nz_sites(imin) + idn
        nst = nst + idn
    enddo
    nz_excess = 0
    excess = 0
    do z = NZ,1,-1
        excess = excess + nz_cells(z) - nz_sites(z)
        nz_excess(z) = excess
    enddo
endif
nz1 = NZ
nz2 = 1
do z = NZ,1,-1
    if (z == NZ) then
        nz_totsites(z) = nz_sites(z)
    else
        nz_totsites(z) = nz_sites(z) + nz_totsites(z+1)
    endif
    if (nz_sites(z) > 0) then
        eratio(z) = 100*nz_excess(z)/real(nz_totsites(z))
        nz1 = min(z,nz1)
        nz2 = max(z,nz2)
    else
        eratio(z) = 0
    endif
enddo

end subroutine


!-----------------------------------------------------------------------------------------
! Instead of adding a site immediately when a B cell enters the follicle, the adjustment
! to the number of sites in the follicle is carried out in balancer(), when the net
! effect of inflows and outflows is accounted for.
!-----------------------------------------------------------------------------------------
subroutine CellInflux(ninflow,ok)
integer :: ninflow
logical :: ok
integer :: k, x, y, z, indx(2), site(3), kcell, ctype, gen, region, kpar=0
real :: R

region = FOLLICLE
gen = 1
k = 0
do while (k < ninflow)
    call getEntrySite(site,ok)
    if (.not.ok) then
		call logger('Error: CellInflux: no entry site found')
		return
	endif
    indx = occupancy(site(1),site(2),site(3))%indx
    if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC
    if (indx(1) == 0 .or. indx(2) == 0) then
        if (evaluate_residence_time) then
            ! This is for measuring residence time
            if (istep > istep_res1 .and. istep <= istep_res2) then
                ctype = RES_TAGGED_CELL   ! Use ctype = 2 = COG_TYPE_TAG for residence time tagging, to have cptr%entrytime
                ninflow_tag = ninflow_tag + 1
            else
                ctype = 1
            endif
        else
            ctype = select_cell_type(kpar)
        endif
        if (ctype /= NONCOG_TYPE_TAG) then
            ncogseed = ncogseed + 1
        endif
        call AddBcell(site,ctype,gen,NAIVE,region,kcell,ok)
        if (.not.ok) then
			call logger('Error: CellInflux: failed to add B cell')
			return
		endif
        k = k+1
        cycle
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Under steady-state conditions (no inflammation) we'd like the number of cells that leave 
! to match approximately the number that enter.
! Egress in this case is possible for any cell at a boundary site, i.e. a site with an
! 'outside' neighbour, in the patch specified by EGRESS_ALPHA.
!-----------------------------------------------------------------------------------------
subroutine traffic(ok)
logical :: ok
integer :: x, y, z, k, kcell, indx(2), ctype, gen, stage, region, status, site(3), n, slot
integer :: node_inflow, node_outflow, add(3), net_inflow, ihr
real(DP) :: R, df, prob
real :: tnow, noncog_exit_prob
logical :: can_leave, left
integer :: kpar=0
integer :: nb, nc, nx
type (boundary_type), pointer :: bdry
type (cog_type), pointer :: p

ok = .true.
!write(*,*) 'traffic'
tnow = istep*DELTA_T
node_inflow = InflowTotal
df = InflowTotal - node_inflow
R = par_uni(kpar)
if (R < df) then
    node_inflow = node_inflow + 1
endif
noncog_exit_prob = base_exit_prob*12/residence_time
! To speed attainment of steady state when inflammation_signal = 0
if (NBcells < NBcells0) then
	noncog_exit_prob = noncog_exit_prob*(real(NBcells)/NBcells0)**3
endif

! Inflow
call cellInflux(node_inflow,ok)
if (.not.ok) return
nadd_sites = nadd_sites + node_inflow

! Outflow
! Traverse the boundary list to locate cells at bdry sites with %exit_ok, i.e. near the lower boundary egress surface,
! and allow some to exit.
node_outflow = 0
nb = 0
nc = 0
nx = 0
bdry => bdrylist
do while ( associated ( bdry )) 
    if (bdry%exit_ok) then
		nb = nb + 1
		site = bdry%site
		indx = occupancy(site(1),site(2),site(3))%indx
		do slot = 2,1,-1
			kcell = indx(slot)
			if (kcell > 0) then
				can_leave = .false.
				p => cellist(kcell)%cptr
				if (associated(p)) then
					status = get_status(p)
					if (status == PLASMA) then
						can_leave = .true.
					endif
				else
					R = par_uni(kpar)
					if (R < noncog_exit_prob) then    ! this cell can leave
						can_leave = .true.
					endif
				endif
				if (can_leave) then
					call CellExit(kcell,slot,site)
					nx = nx + 1
					node_outflow = node_outflow + 1
					if (evaluate_residence_time) then
						if (cellist(kcell)%ctype == RES_TAGGED_CELL) then
							noutflow_tag = noutflow_tag + 1
							restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
							ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
							Tres_dist(ihr) = Tres_dist(ihr) + 1
						endif
					endif
				endif
			endif
		enddo
	endif
    bdry => bdry%next
enddo
nadd_sites = nadd_sites - node_outflow
if (mod(istep,4*60) == 0) then
	write(nfout,*) nb,nx,noncog_exit_prob
endif
return

! Unused code follows
k = 0
do while (k < node_outflow)
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out x R: ',R
    x = 1 + R*NX
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out y R: ',R
    y = 1 + R*NY
    R = par_uni(kpar)
    if (dbug) write(nfres,'(a,f10.6)') 'out z R: ',R
    z = 1 + R*NZ        ! any z is OK to exit
    indx = occupancy(x,y,z)%indx
    if (indx(1) < 0) cycle      ! OUTSIDE_TAG or DC 
    if (indx(2) > 0) then
        slot = 2
    elseif (indx(1) > 0) then
        slot = 1
    else
        cycle
    endif
    kcell = indx(slot)
    site = (/x,y,z/)
    call CellExit(kcell,slot,site)
!    if (.not.left) cycle
    k = k+1

    if (evaluate_residence_time) then
        if (cellist(kcell)%ctype == RES_TAGGED_CELL) then
            noutflow_tag = noutflow_tag + 1
            restime_tot = restime_tot + tnow - cellist(kcell)%entrytime
            ihr = (tnow - cellist(kcell)%entrytime)/60. + 1
            Tres_dist(ihr) = Tres_dist(ihr) + 1
        endif
    endif
enddo
nadd_sites = nadd_sites - node_outflow
!write(*,'(a,4i6)') 'istep,in,out: ',istep,node_inflow,node_outflow,NBcells
end subroutine

!-----------------------------------------------------------------------------------------
! Determine whether a cognate T cell is licensed to exit the DCU.
! Possible exit rules:
! (1) gen >= NGEN_EXIT
! (2) gen > 1 and act < EXIT_THRESHOLD
! (3) Allow exit to any cognate cell that gets there (was gen > 1 and CD69 < CD69_threshold)
!-----------------------------------------------------------------------------------------
logical function exitOK(p,ctype)
type(cog_type), pointer :: p
integer :: ctype

exitOK = .true.
if (ctype == COG_TYPE_TAG) then
    exitOK = .false.
endif
end function



!-----------------------------------------------------------------------------------------
! The cell kcell in slot at esite(:) is a candidate to exit.  If it meets the criteria,
! left = .true.
! Currently, when a T cell exits the LN its ID is set to 0, and a gap is recorded in cellist(:).
!-----------------------------------------------------------------------------------------
subroutine CellExit(kcell,slot,esite)
integer :: kcell, slot, esite(3)
!logical :: left
integer :: x, y, z, ctype, gen, stage, region, status
real :: tnow
logical :: cognate, activated
type(cog_type), pointer :: p

tnow = istep*DELTA_T
!left = .false.
if (evaluate_residence_time) then
    cognate = .false.
elseif (associated(cellist(kcell)%cptr)) then
    cognate = .true.
    p => cellist(kcell)%cptr
!	call get_stage(p,stage,region)
	stage = get_stage(p)
	status = get_status(p)
    gen = get_generation(p)
    ctype = cellist(kcell)%ctype
!    if (.not.exitOK(p,ctype)) then
!        write(*,*) 'cell_exit: cognate exit suppressed: ',kcell
!        return
!    endif
else
    cognate = .false.
endif
! For initial testing, remove cells that leave the LN
x = esite(1)
y = esite(2)
z = esite(3)
if (slot == 2) then
    occupancy(x,y,z)%indx(2) = 0
else
    if (occupancy(x,y,z)%indx(2) == 0) then
        occupancy(x,y,z)%indx(1) = 0
    else    ! shift cell in slot 2 to slot 1 (indx(2) to indx(1))
        occupancy(x,y,z)%indx(1) = occupancy(x,y,z)%indx(2)
        occupancy(x,y,z)%indx(2) = 0
    endif
endif
NBcells = NBcells - 1
if (cognate) then
!    if (.not.evaluate_residence_time .and. activated) then
!		call efferent(p,ctype)
!	endif
	ngaps = ngaps + 1
	gaplist(ngaps) = kcell
	cellist(kcell)%ID = 0
    cognate_list(p%cogID) = 0
    write(logmsg,'(a,3i6)') 'CellExit: cognate cell left: status,stage: ',kcell,status,stage
    call logger(logmsg)
else
	ngaps = ngaps + 1
	gaplist(ngaps) = kcell
	cellist(kcell)%ID = 0
endif
!left = .true.
end subroutine

!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------- 
subroutine initialise_vascularity
if (VEGF_MODEL == 2) then	! Not used
	VEGF_beta = 0
	VEGF_baserate = 0
    VEGF_decayrate = 0.0012
    vasc_maxrate = 0.001
    VEGFmass = 0
    vasc_beta = 0.00001
    vasc_decayrate = 0.001
else	! VEGF_MODEL = 1
!	VEGF_beta = 4.0e-8
	VEGF_baserate = VEGF_beta*NBcells0
!    VEGF_decayrate = 0.002		! delta_G
!    vasc_maxrate = 0.0006		! alpha_V
    VEGFmass = VEGF_baserate/VEGF_decayrate    ! steady-state VEGF level M_G0
    Cvegf0 = VEGFmass/NBcells0	! taking K_V = 1.0
!    vasc_beta = 1.5				! beta_V
    vasc_decayrate = vasc_maxrate*hill(Cvegf0,vasc_beta*Cvegf0,vasc_n)	! delta_V
 !   write(*,*) 'Vascularity parameters:'
 !   write(*,*) 'alpha_G,beta_G,delta_G: ',VEGF_alpha, VEGF_beta, VEGF_decayrate
 !   write(*,*) 'alpha_V,beta_V,delta_V: ',vasc_maxrate,vasc_beta,vasc_decayrate
 !   write(*,*) 'Cvegf0,VEGF0,VEGF_baserate: ',Cvegf0,VEGFmass,VEGF_baserate
 !   write(*,*) 'vasc_beta*Cvegf0: ',vasc_beta*Cvegf0
endif
Vascularity = 1.00
!write(*,*) 'VEGF_MODEL, VEGF_baserate: ',VEGF_MODEL, VEGF_baserate
end subroutine

!-----------------------------------------------------------------------------------------
! Vascularity responds to the VEGF level.  The rate of production of VEGF is proportional to
! either:
! (VEGF_MODEL_1) a specified inflammation signal
! or
! (VEGF_MODEL_2) the total DC activity level (i.e. antigen load).
!
! In Model 1
!	VEGFsignal depends on the inflammation level
! In Model 2
!   VEGFsignal depends on total DC activity = total antigen density (normalized)
!
!   VEGF_baserate = constitutive rate of VEGF production (per LN volume, i.e. per T cell)
!   dVEGFdt = current rate of VEGF secretion (artificial units)
!   VEGF = current total mass of VEGF (artificial units)
!   Cvegf = current VEGF concentration (artificial units)
!   vasc_beta, vasc_n = Hill function parameters for dependence of vascularity growth
!                    on VEGF concentration
!   vasc_maxrate = maximum rate constant for growth of relative vascularity
!   vasc_decayrate = rate constant for decline of relative vascularity
!   Vascularity = current relative vascularity
!   Note: if inflammation = 0, i.e. VEGFsignal = 0, there should be no change in vascularity,
!   i.e. we should have dVdt = 0.
!-----------------------------------------------------------------------------------------
subroutine vascular
real :: VEGFsignal=0, dVEGFdt, Nfactor
real :: Cck1 = 1.5, Cck2 = 0.1

!write(*,*) 'vascular'
if (.not.vary_vascularity) then
    Vascularity = 1.0
    return
endif
VEGFsignal = get_inflammation() ! Rate of secretion of VEGF is proportional to inflammation  
! This is a measure to speed up attainment of steady-state
if (VEGFsignal == 0 .and. NBcells < NBcells0) then
	Cvegf = Cvegf0
	VEGFmass = NBcells0*Cvegf
    Vascularity = real(NBcells0)/NBcells
    return
endif
!Nfactor = real(NBcells)/NBcells0
dVEGFdt = VEGFsignal*VEGF_alpha + VEGF_baserate - VEGF_decayrate*VEGFmass
! Mass of VEGF is augmented by rate, and is subject to decay
VEGFmass = VEGFmass + dVEGFdt*DELTA_T
Cvegf = VEGFmass/NBcells   ! concentration (proportional to, anyway) 
if (VEGF_MODEL == 2) then ! not used
    dVdt = vasc_maxrate*hill(Cvegf,vasc_beta,vasc_n)*Vascularity - vasc_decayrate*(Vascularity - 1)
else	! VEGF_MODEL = 1
    dVdt = vasc_maxrate*hill(Cvegf,vasc_beta*Cvegf0,vasc_n)*Vascularity - vasc_decayrate*Vascularity	!  this works!
endif
Vascularity = max(Vascularity + dVdt*DELTA_T, 1.0)
!write(*,*) 'dVEGFdt, Cvegf, dVdt: ',dVEGFdt, Cvegf, dVdt, Vascularity
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function hill(x,b,n)
real :: x, b
integer :: n
hill = x**n/(x**n + b**n)
end function

!-----------------------------------------------------------------------------------------
! The total number of B cells added to the blob since the last balancing is nadd_sites.
! (Note that this can be negative, indicating a net loss of B cells). 
! A balancing is triggered either when this count exceeds a limit nadd_limit, which is a
! specified fraction of the original B cell population NBcells0, or when the time since
! the last balancing exceeds BALANCER_INTERVAL, and an adjustment to site count is needed.
! Sites needed or made available by a change to the DC population are accounted for.
! Either new sites are added (made available) or existing sites are removed (made OUTSIDE).
!-----------------------------------------------------------------------------------------
subroutine balancer(ok)
logical :: ok
integer :: nadd_total, nadd_limit, n, nadded
integer :: k, idc, naddDC, naddex, nremex, Nexits0, dexit
real :: tnow
integer :: kpar = 0
logical :: blob_changed

ok = .true.
!write(*,*) 'balancer: ',istep

blob_changed = .false.
tnow = istep*DELTA_T
nadd_limit = 0.01*NBcells0
nadd_total = nadd_sites

if ((abs(nadd_total) > nadd_limit) .or. (tnow > lastbalancetime + BALANCER_INTERVAL)) then
    if (dbug) write(nflog,*) 'balancer: nadd_total: ',nadd_total,nadd_limit,lastbalancetime,BALANCER_INTERVAL
    if (dbug) write(nflog,*) 'call squeezer'
	call squeezer(.false.)
	if (dbug) write(nflog,*) 'did squeezer'
    if (dbug) write(nflog,*) 'nadd_total: ',nadd_total
    if (nadd_total > 0) then
        n = nadd_total
	    if (dbug) write(nflog,*) 'call add_sites'
        call AddSites(n,ok)
        if (.not.ok) return
        n = 0
	    if (dbug) write(nflog,*) 'did add_sites'
        blob_changed = .true.
    elseif (nadd_total < 0) then
        n = -nadd_total
	    if (dbug) write(nflog,*) 'call removeSites'
        call RemoveSites(n,ok)
        if (.not.ok) return
        n = 0
	    if (dbug) write(nflog,*) 'did removeSites'
!	    write(logmsg,*) 'did removeSites: exit #10: ',exitlist(10)%site,exitlist(5)%site
!	    call logger(logmsg) 
!	    call checkexits("after removeSites")
        blob_changed = .true.
    else
        n = 0
    endif
    if (n /= 0) then
        write(logmsg,'(a,4i6)') 'Error: balancer: remaining sites: ',n, &
			NBcells,Nsites
		call logger(logmsg)
        ok = .false.
        return
    endif
! The cognate list is maintained at the time that a cell arrives or leaves
    if (NBcells /= Nsites) then
	    write(logmsg,'(a,3i10)') 'Error: balancer: cells /= sites: ', &
			NBcells,Nsites
	    call logger(logmsg)
	    ok = .false.
	    return
	endif
    lastbalancetime = tnow
    nadd_sites = 0
!    if (blob_changed) then
!		if (dbug) write(nflog,*) 'call make_split'
!        call make_split
!        if (use_diffusion) then
!            call setup_minmax
!        endif
!        if (dbug) write(nflog,'(a,2i6)') 'balancer: nadd_total, radius: ',nadd_total,int(aRadius)
!    endif
!    if (use_portal_egress) then
!		call adjustExits
!	endif
else
    call set_globalvar
endif
return


end subroutine




!--------------------------------------------------------------------------------
! Counts efferent cognate cells, and records their generation distribution.
! Only activated cells (stage >= CLUSTERS) are counted
!--------------------------------------------------------------------------------
subroutine efferent(p,ctype)
type (cog_type), pointer :: p
integer :: ctype, gen, i
real :: avid

gen = get_generation(p)
if (ctype > NCTYPES) then
    write(*,*) 'efferent: bad cell type:', ctype
    stop
endif
!localres%dN_EffCogTC(ctype)  = localres%dN_EffCogTC(ctype) + 1
!localres%dN_EffCogTCGen(gen) = localres%dN_EffCogTCGen(gen) + 1
!localres%N_EffCogTC(ctype)   = localres%N_EffCogTC(ctype) + 1
!localres%N_EffCogTCGen(gen)  = localres%N_EffCogTCGen(gen) + 1
totalres%dN_EffCogTC(ctype)  = totalres%dN_EffCogTC(ctype) + 1
totalres%dN_EffCogTCGen(gen) = totalres%dN_EffCogTCGen(gen) + 1
totalres%N_EffCogTC(ctype)   = totalres%N_EffCogTC(ctype) + 1
totalres%N_EffCogTCGen(gen)  = totalres%N_EffCogTCGen(gen) + 1

if (log_results) then
    ! Record avidity statistics for exiting cells
    avid = p%avidity
    if (avid_count%logscale) then
        avid = log10(avid)
    endif
    if (avid_count%nbins == 1) then
        i = 1
    else
        !i = (avid-avidity_min)*1.01/avidity_step + 1
        i = (avid-avid_count%binmin)/avid_count%binstep + 1.5
        i = max(i,1)
        !i = min(i,avidity_nlevels)
        i = min(i,avid_count%nbins)
    endif
    avid_count%ndist(i) = avid_count%ndist(i) + 1
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine show_snapshot(ok)
logical :: ok
integer :: kcell, ctype, stype, ncog, noncog, ntot, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead
real :: stim(FINISHED), IL2sig(FINISHED), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

ok = .true.
return

allocate(gendist(BC_MAX_GEN))
allocate(div_gendist(BC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
!		call get_stage(p,stage,region)
        stage = get_stage(p)
        nst(stage) = nst(stage) + 1
        stim(stage) = stim(stage) + p%stimulation
!        IL2sig(stage) = IL2sig(stage) + get_IL2store(p)
        gen = get_generation(p)
        if (gen < 0 .or. gen > BC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'show_snapshot: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
        if ((gen == 1 .and. p%stimulation > FIRST_DIVISION_THRESHOLD(1)) .or. &
			(gen > 1 .and. p%stimulation > DIVISION_THRESHOLD(1))) then
			div_gendist(gen) = div_gendist(gen) + 1
        endif
        max_TCR = max(p%stimulation,max_TCR)
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad stype: ',ctype,stype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
        IL2sig(i) = IL2sig(i)/nst(i)
    else
        stim(i) = 0
        IL2sig(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = BC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:BC_MAX_GEN))
do i = BC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
!mols_pM = L_um3*M_pM/(NBcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif
if (.not.use_TCP .and. use_cognate) then
write(*,'(a)') '----------------------------------------------------------------------'
write(*,*) 'use_cognate: ',use_cognate
write(*,'(a,i6,3i8,a,2i8)') 'snapshot: ',istep,ntot,ncogseed,ncog,'     dead: ',dNdead,Ndead
write(*,'(a,7i7)')   '# in stage:  ',nst
write(*,'(a,7f7.0)') 'stimulation: ',stim
write(*,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogTC(2:3),sum(totalres%dN_EffCogTC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogTC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogTCGen(1:neffgens)
!if (use_cytokines) then
!    do iseq = 1,Ncytokines
!        if (use_diffusion) then
!            cyt_conc = cyt_mean(iseq)
!        else
!            cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
!        endif
!        write(*,'(3a,f8.4)') 'Mean cytokine conc: ',cyt_name(cyt_tag(iseq)),'  ',cyt_conc
!    enddo
!endif
!write(*,'(a)') '----------------------------------------------------------------------'

!call check_cognate_list
!kcog = 1
!kcell = cognate_list(kcog)
!if (kcell > 0) then
!    write(*,'(2i6,f8.4,f8.1)') kcog,kcell,cellist(kcell)%cptr%stimrate,cellist(kcell)%cptr%CD69
!endif
!write(*,*) '========= Average time to IL2 threshold: ',nIL2thresh,tIL2thresh/max(1,nIL2thresh)
!write(*,'(a)') '----------------------------------------------------------------------'
endif

!call get_cognate_dist(ncog1,ncog2)

!write(nfout,'(i8,f8.0,i8,6i8,25f7.4)') istep,tnow/60,0,ntot,ncogseed,ncog,dNdead,Ndead,teffgen, &
!    fac*totalres%N_EffCogTCGen(1:TC_MAX_GEN)
if (use_tcp) then
!    if (.not.awp_1%is_open) then
!        call logger("in show_snapshot: awp_1 is not open")
!    endif
!    write(msg,'(2(i6,f8.0),5i8)') istep,tnow,NDCalive,act,ntot,ncogseed,ncog,Ndead,teffgen
!    call winsock_send(awp_1,msg,len_trim(msg),error)
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
    call logger(msg)
endif

! To plot outflow variation with time
!write(nfout,'(2f8.2)') tnow/60,OutflowTotal

!write(nfout,'(a,7f7.0)') 'stimulation: ',stim
!write(nfout,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
deallocate(gendist)
totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_stim_dist
integer :: kcell, ctype, stype, ncog, noncog
real :: s
type (cog_type), pointer :: p

ncog = 0
noncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        s = p%stimulation
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: compute_stim_dist: bad stype: ',ctype,stype
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_cognate_list
integer :: k, kcell

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > nlist) then
        write(*,*) 'check_cognate_list: bad cognate_list entry > nlist: ',lastcogID,k,kcell,nlist
        stop
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The aim is to display the distribution of cognate cells in the z direction, in order
! to investigate failure of cognate cell proliferation to scale with NBcells0.
! For now, just look at cognate fractions in upper and lower hemispheres.
!-----------------------------------------------------------------------------------------
subroutine get_cognate_dist(ncog1,ncog2)
integer :: z,kcell,ntot1,ntot2,ncog1,ncog2
logical :: cognate

ntot1 = 0
ntot2 = 0
ncog1 = 0
ncog2 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    z = cellist(kcell)%site(3)
    cognate = (associated(cellist(kcell)%cptr))
    if (z < z0) then
        ntot1 = ntot1 + 1
        if (cognate) ncog1 = ncog1 + 1
    else
        ntot2 = ntot2 + 1
        if (cognate) ncog2 = ncog2 + 1
    endif
enddo
write(*,*) 'Lower: ',ntot1,ncog1,real(ncog1)/ntot1
write(*,*) 'Upper: ',ntot2,ncog2,real(ncog2)/ntot2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testqsort(n)
integer :: n
integer :: i
integer :: kpar = 0
real, allocatable :: a(:)
integer, allocatable :: t(:)

call rng_initialisation

allocate(a(n))
allocate(t(n))
do i = 1,n
    a(i) = par_uni(kpar)
    t(i) = i
enddo
call qsort(a,n,t)
end subroutine



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine vascular_test
integer :: nsteps = 10.0*24*60/DELTA_T
real :: inflow0, act, tnow, nsum, Fin0, Fin, Fout, Tres

Tres = 24
NBcells0 = 100000
use_chemotaxis = .false.
call initialise_vascularity

NBcells = NBcells0
Fin0 = NBcells/(Tres*60)
write(*,*) nsteps,NBcells

nsum = 0
do istep = 1,nsteps
	tnow = istep*DELTA_T
    call vascular
!    call generate_traffic(inflow0)
	Fin = Fin0*Vascularity*DELTA_T
	Fout = NBcells*DELTA_T/(Tres*60)
    if (mod(istep,60) == 0) then
        write(nfout,'(i6,f8.3,5e14.5,i8)') istep,tnow, &
			VEGFmass, &
            dVdt, &
            Vascularity, &
            Fin, Fout, &
            NBcells
    endif
    nsum = nsum + (Fin - Fout)
    NBcells = NBcells0 + nsum
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testrnor
integer :: n = 100000
integer :: k, j
integer :: kpar = 0
real :: r, rmin=1.0e10, rmax = -1.0e10

do k = 1,n
    do j = 1,n
        r = par_rnor(kpar)
        rmin = min(r,rmin)
        rmax = max(r,rmax)
    enddo
    write(*,'(i12,2e12.4)') k,rmin,rmax
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testrandom
integer :: n = 100000
integer :: k, j, ncog
integer :: kpar = 0

ncog = 0
do k = 1,n
    j = select_cell_type(kpar)
    if (j == COG_CD4_TAG) ncog = ncog + 1
    if (mod(k,1000) == 0) then
        write(*,*) ncog,real(ncog)/k
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine old_process_command_line
integer :: i, cnt, len, status
character :: c*(64), b*(256)
character*(64) :: progname

call get_command (b, len, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b (1:len)
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c (1:len)
progname = c(1:len)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
!    Mnodes = 4      ! for profiling
!    write(*,*) 'Ruuning with Mnodes = ',Mnodes
endif
do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
!    write (*,*) 'command arg ', i, ' = ', c (1:len)
    if (i == 1) then
!!!        read(c(1:len),'(i)') Mnodes
        read(c(1:len),*) Mnodes
        write(*,*) 'Requested threads: ',Mnodes
    elseif (i == 2) then
        inputfile = c(1:len)
        write(*,*) 'Input file: ',inputfile
    elseif (i == 3) then
        outputfile = c(1:len)
        write(*,*) 'Output file: ',outputfile
!    elseif (i == 4) then
!        resultfile = c(1:len)
!        write(*,*) 'Result file: ',resultfile 
    endif
end do
write (*,*) 'command line processed'
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine write_Tres_dist
integer :: k, kmax
real(8) :: Tres

do k = int(days*24),1,-1
    if (Tres_dist(k) /= 0) then
        kmax = k
        exit
    endif
enddo
Tres = 0
do k = 1,kmax
    Tres = Tres + (k-0.5)*Tres_dist(k)/noutflow_tag
enddo
write(nfout,'(a)') 'Transit time distribution'
write(nfout,'(a)') 'Parameters: '
if (TAGGED_EXIT_CHEMOTAXIS) then
	write(nfout,'(a,L)') '  TAGGED_EXIT_CHEMOTAXIS:     ',TAGGED_EXIT_CHEMOTAXIS
	write(nfout,'(a,f6.3)') '  TAGGED_CHEMO_FRACTION: ',TAGGED_CHEMO_FRACTION 
	write(nfout,'(a,f6.3)') '  TAGGED_CHEMO_ACTIVITY: ',TAGGED_CHEMO_ACTIVITY
endif
write(nfout,'(a,f6.3)') 'chemo_K_exit:            ',chemo_K_exit
write(nfout,'(a,f6.3)') '  K1_S1P1:               ',K1_S1P1
write(nfout,'(a,3i8)') 'Results: noutflow_tag, ninflow_tag,kmax: ',noutflow_tag,ninflow_tag,kmax
write(nfout,'(a,f6.1)') 'Residence time: ',Tres
!write(nfout,'(a,f6.1)') 'Residence time from restime_tot: ',restime_tot/60.
write(nfout,'(a)') ' Hour   Probability'
do k = 1,kmax
    write(nfout,'(f6.1,e12.4)') k-0.5,Tres_dist(k)/noutflow_tag
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The distributions gendist, tcrdist and aviddist are all for the current DCU population.
! avid_count is for the recent efferent population.
!-----------------------------------------------------------------------------------------
subroutine write_results
integer :: kcell, ctype, stype, ncog, ntot, i
integer :: gen
real :: tcr, avid, dtcr, hour
type (cog_type), pointer :: p
integer :: gendist(BC_MAX_GEN),aviddist(MAX_AVID_LEVELS),tcrdist(tcr_nlevels)
character*(60) :: fmtstr = '(f6.2,2i8,4x,15f7.4,4x,10f7.4,4x,10f7.4,4x,10i7)'

write(fmtstr(14:15),'(i2)') BC_MAX_GEN
write(fmtstr(24:25),'(i2)') tcr_nlevels
write(fmtstr(34:35),'(i2)') avidity_nlevels
write(fmtstr(44:45),'(i2)') avidity_nlevels
hour = istep*DELTA_T/60
dtcr = TCR_limit/TCR_nlevels
gendist = 0
aviddist = 0
tcrdist = 0
ntot = 0
ncog = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    ntot = ntot + 1
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog = ncog + 1
        ! TCR stimulation distribution
        tcr = p%stimulation
        i = tcr/dtcr + 1
        i = min(i,TCR_nlevels)
        tcrdist(i) = tcrdist(i) + 1
        ! T cell generation distribution
        gen = get_generation(p)
        gendist(gen) = gendist(gen) + 1
        ! T cell avidity distribution
        avid = p%avidity
        if (avidity_logscale) then
            avid = log10(avid)
        endif
        if (avidity_nlevels == 1) then
            i = 1
        else
            i = (avid-avidity_min)/avidity_step + 1.5
!           write(nfout,'(a,2f8.4,3i7)') 'Count: ',p%avidity,avid,i,kcell,p%cogID
            i = max(i,1)
            i = min(i,avidity_nlevels)
        endif
        aviddist(i) = aviddist(i) + 1
    endif
enddo
if (fix_avidity) then
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog,real(aviddist)/ncog, &
        avid_count%ndist
else
    write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog,real(tcrdist)/ncog
endif
avid_count_total%ndist = avid_count_total%ndist + avid_count%ndist
avid_count%ndist = 0
!write(nfout,'(8i6)') aviddist
end subroutine

!-----------------------------------------------------------------------------------------
! A different approach to tagging cells for residence time computation
!-----------------------------------------------------------------------------------------
subroutine tag_cells
integer :: kcell
type(cell_type),pointer :: cell

do kcell = 1,nlist
    cell => cellist(kcell)
    if (cell%ID == 0) cycle
    cell%ctype = RES_TAGGED_CELL
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Various logging counters are initialized here.
!-----------------------------------------------------------------------------------------
subroutine init_counters

ninflow_tag = 0
noutflow_tag = 0
restime_tot = 0
if (log_results) then
    if (.not.use_cognate) then
        write(*,*) 'No use logging results with no cognate cells'
        stop
    endif
    avid_count%nbins = avidity_nlevels
    allocate(avid_count%ndist(avid_count%nbins))
    avid_count%period = ntres
    avid_count%logscale = avidity_logscale
    avid_count%binmin = avidity_min
    avid_count%binstep = avidity_step
    avid_count%ndist = 0
    avid_count%total = 0

    allocate(avid_count_total%ndist(avid_count%nbins))
    avid_count_total = avid_count
!    avid_count_total%nbins = avid_count%nbins
!    avid_count_total%period = avid_count%period
!    avid_count_total%logscale = avid_count%logscale
!    avid_count_total%binmin = avidity_min
!    avid_count_total%binstep = avidity_step
!    avid_count_total%ndist = 0
!    avid_count_total%total = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_chemoactivity(cave)
real :: cave
type (cell_type), pointer :: cell
integer :: kcell

cave = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    cell => cellist(kcell)
    cave = cave + chemo_active_exit(cell)
enddo
cave = cave/NBcells
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_cum_prob
real :: m = 30, s = 2.0
real :: p1, p2, a
integer :: i

p1 = log(m)
p2 = log(s)
do i = 1,30
    a = 2*i
    write(*,'(i4,2f8.3)') i,a,1-cum_prob_lognormal(a,p1,p2)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_pause
!!!use ifport
logical :: paused

inquire(file=pausefile,exist=paused)
if (paused) then
	call logger('Pause order received')
	do while (paused)
!!!		call sleepqq(100)
        call sleeper(1)   ! Too coarse!
		inquire(file=pausefile,exist=paused)
	enddo
	call logger('Resuming ...')
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine SaveGenDist
integer :: gendist(BC_MAX_GEN)
integer :: k, kcell, stage, region, gen, maxg
real :: fac

gendist = 0
maxg = 0
do k = 1,lastcogID
	kcell = cognate_list(k)
	if (kcell > 0) then
!		call get_stage(cellist(kcell)%cptr,stage,region)
!		if (region /= FOLLICLE) cycle 
		gen = get_generation(cellist(kcell)%cptr)
		maxg = max(maxg,gen)
		gendist(gen) = gendist(gen) + 1
	endif
enddo
fac = 1.0/sum(gendist)
write(nflog,*) 'gendist:'
write(nflog,'(20f7.4)') fac*gendist(1:maxg)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim,NY_dim,NZ_dim) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine get_summary_old(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)
logical :: ok
integer :: kcell, ctype, stype, ncog(2), noncog, ntot, nbnd, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead, nact
real :: stim(2*STAGELIMIT), IL2sig(2*STAGELIMIT), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

if (firstSummary) then
	write(nfout,'(a)') "==========================================================================="
	firstSummary = .false.
endif
ok = .true.

allocate(gendist(TC_MAX_GEN))
allocate(div_gendist(TC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nbnd = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    if (associated(p)) then
!		call get_stage(p,stage,region)
		stage = get_stage(p)
		region = get_region(p)
	else
		stage = 0
		region = FOLLICLE
	endif
	if (region == FOLLICLE) then
	    ntot = ntot + 1
	endif
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        ncog(region) = ncog(region) + 1
        nst(stage) = nst(stage) + 1
        stim(stage) = stim(stage) + p%stimulation
!        IL2sig(stage) = IL2sig(stage) + get_IL2store(p)
        gen = get_generation(p)
        if (gen < 0 .or. gen > TC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'get_summary: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
        if ((gen == 1 .and. p%stimulation > FIRST_DIVISION_THRESHOLD(1)) .or. &
			(gen > 1 .and. p%stimulation > DIVISION_THRESHOLD(1))) then
			div_gendist(gen) = div_gendist(gen) + 1
        endif
        max_TCR = max(p%stimulation,max_TCR)
    elseif (stype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad stype: ',ctype,stype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
        IL2sig(i) = IL2sig(i)/nst(i)
    else
        stim(i) = 0
        IL2sig(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = TC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:TC_MAX_GEN))
do i = TC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
!mols_pM = L_um3*M_pM/(NBcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif
if (.not.use_TCP .and. use_cognate) then
write(*,'(a)') '----------------------------------------------------------------------'
write(*,'(a,i6,4i8,a,2i8)') 'snapshot: ',istep,ntot,ncogseed,ncog,'     dead: ',dNdead,Ndead
write(*,'(a,7i7)')   '# in stage:  ',nst
write(*,'(a,7f7.0)') 'stimulation: ',stim
write(*,'(a,7f7.0)') 'IL2 signal:  ',IL2sig
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogTC(2:3),sum(totalres%dN_EffCogTC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogTC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogTCGen(1:neffgens)
write(*,'(a,f6.0)') 'max_TCR: ',max_TCR
!if (use_cytokines) then
!    do iseq = 1,Ncytokines
!        if (use_diffusion) then
!            cyt_conc = cyt_mean(iseq)
!        else
!            cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
!        endif
!        write(*,'(3a,f8.4)') 'Mean cytokine conc: ',cyt_name(cyt_tag(iseq)),'  ',cyt_conc
!    enddo
!endif
write(*,'(a)') '----------------------------------------------------------------------'

write(*,*) '========= Average time to IL2 threshold: ',nIL2thresh,tIL2thresh/max(1,nIL2thresh)
write(*,'(a)') '----------------------------------------------------------------------'
endif

if (use_tcp) then
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
    call logger(msg)
endif
deallocate(gendist)
deallocate(div_gendist)

totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

!write(nfout,'(i4,i8,i4,4i8,4i6,i8,25f7.4)') int(tnow/60),istep,0,ntot,ncogseed,ncog,Ndead, &
!	nbnd,int(InflowTotal),Nexits, teffgen, fac*totalres%N_EffCogTCGen(1:TC_MAX_GEN)
summaryData(1:13) = (/int(tnow/60),istep,0,0,ntot,ncogseed,ncog,Ndead, &
	nbnd,int(InflowTotal),Nexits, teffgen/)
write(nflog,*) 'ndivisions = ',ndivisions

write(logmsg,'(a,i5,a,2i4,i6,i8,100i4)') 'In: ',check_inflow,' Out: ',Nexits,Lastexit,sum(check_egress),NBcells
	!,(check_egress(i),i=1,Lastexit)
call logger(logmsg)
check_inflow = 0
check_egress = 0
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, cellist(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)
logical :: ok
integer :: kcell, ctype, stype, ncog, noncog, ntot, nbnd, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead, nvasc
real :: stim(2*STAGELIMIT), IL2sig(2*STAGELIMIT), tgen, tnow, fac, act, cyt_conc, mols_pM
type (cog_type), pointer :: p
integer :: nst(FINISHED)
integer, allocatable :: gendist(:)
integer, allocatable :: div_gendist(:)  ! cells that are capable of dividing
character*(6) :: numstr
character*(256) :: msg

if (firstSummary) then
	write(nfout,'(a)') "==========================================================================="
	firstSummary = .false.
endif
ok = .true.

allocate(gendist(BC_MAX_GEN))
allocate(div_gendist(BC_MAX_GEN))
tnow = istep*DELTA_T
noncog = 0
ncog = 0
ntot = 0
nbnd = 0
nst = 0
stim = 0
IL2sig = 0
gendist = 0
div_gendist = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    p => cellist(kcell)%cptr
    if (associated(p)) then
!		call get_stage(p,stage,region)
		stage = get_stage(p)
		region = get_region(p)
	else
		stage = 0
		region = FOLLICLE
	endif
	if (region == FOLLICLE) then
	    ntot = ntot + 1
	else
		write(logmsg,*) 'kcell, region: ',kcell,region
		call logger(logmsg)
		stop
	endif
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
    if (ctype == COG_TYPE_TAG) then
        ncog = ncog + 1
!        nst(stage) = nst(stage) + 1
!        stim(stage) = stim(stage) + p%stimulation
        gen = get_generation(p)
        if (gen < 0 .or. gen > BC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'get_summary: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
    elseif (ctype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    else
        write(logmsg,*) 'ERROR: get_summary: bad stype: ',ctype,stype
        call logger(logmsg)
        stop
    endif
enddo

tgen = sum(gendist)
do i = BC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogTCGen(1:BC_MAX_GEN))
do i = BC_MAX_GEN,1,-1
    if (totalres%N_EffCogTCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead
!mols_pM = L_um3*M_pM/(NBcells*Vc*Navo)

if (teffgen > 0) then
    fac = 1/real(teffgen)
else
    fac = 0
endif

if (use_tcp) then
    msg = ''
    do i = 1,ngens
		write(numstr,'(i6)') gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//'('
		write(numstr,'(i6)') div_gendist(i)
		msg = trim(msg)//trim(adjustl(numstr))
		msg = trim(msg)//')-'
	enddo
!    call logger(msg)	! msg sending is suppressed for now 
endif
deallocate(gendist)
deallocate(div_gendist)

totalres%dN_EffCogTC = 0
totalres%dN_EffCogTCGen = 0
totalres%dN_Dead = 0

summaryData(1:9) = (/int(tnow/60),istep,ntot,ncogseed,ncog,Ndead,int(InflowTotal*60/DELTA_T), int(100*vascularity), teffgen/)
check_inflow = 0
check_egress = 0
end subroutine

!-------------------------------------------------------------------------------- 
! The GUI calls this subroutine to fetch the cell info needed to identify and render 
! the cells:
!   id			the cell's sequence number
!   position	(x,y,z)
!   state       this is translated into a colour
!
! The info is stored in integer arrays, one for B cells, one for FDCs,
! and one for cell-cell bonds (not used).
! As a quick-and-dirty measure, the first 7 B cells in the list are actually 
! markers to provide a visual indication of the extent of the follicular blob.
!--------------------------------------------------------------------------------
subroutine get_scene(nBC_list,BC_list,nFDC_list,FDC_list,nbond_list,bond_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nFDC_list, nBC_list, nbond_list, FDC_list(*), BC_list(*), bond_list(*)
integer :: k, kc, kcell, site(3), j, jb, idc, fdcsite(3)
integer :: r, g, b
integer :: x, y, z
real :: dcstate
integer :: ifdcstate, ibcstate, stype, ctype, stage, region
real :: bcell_diam = 0.9
real :: FDC_diam = 1.8
integer :: gen, bnd(2)
integer, parameter :: axis_centre = -2	! identifies the ellipsoid centre
integer, parameter :: axis_end    = -3	! identifies the ellipsoid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the ellipsoid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 5			! the size of the info package for a cell (number of integers)

k = 0
! Need some markers to delineate the follicle extent.  These 7 "cells" are used to convey the follicle centre
! and the approximate ellipsoidal blob limits in the 3 axis directions.
do k = 1,7
	select case (k)
	case (1)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_centre
	case (2)
		x = Centre(1) - aRadius - 2
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (3)
		x = Centre(1) + aRadius + 2
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (4)
		x = Centre(1) + 0.5
		y = Centre(2) - bRadius - 2
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_bottom
	case (5)
		x = Centre(1) + 0.5
		y = Centre(2) + bRadius + 2
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (6)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
		z = Centre(3) - bRadius - 2
		site = (/x, y, z/)
		ibcstate = axis_end
	case (7)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
		z = Centre(3) + bRadius + 2
		site = (/x, y, z/)
		ibcstate = axis_end
	end select

	j = ninfo*(k-1)
	BC_list(j+1) = k-1
!	write(logmsg,*) 'cell list #: ',j+1,k-1,site
!	call logger(logmsg)
	BC_list(j+2:j+4) = site
	BC_list(j+5) = ibcstate
enddo
k = k-1

! B cell section
do kc = 1,lastcogID
	kcell = cognate_list(kc)
	if (kcell > 0) then
		region = get_region(cellist(kcell)%cptr)
		if (region /= FOLLICLE) cycle
		k = k+1
		j = ninfo*(k-1)
		site = cellist(kcell)%site
		call BcellColour(kcell,r,g,b)
!		gen = get_generation(cellist(kcell)%cptr)
!		ibcstate = gen
		! Need bcstate to convey non-activated status, i.e. 0 = non-activated
		BC_list(j+1) = kc-1 + 7
!	write(logmsg,*) 'cell list #: ',j+1,kc-1+7,site
!	call logger(logmsg)
		BC_list(j+2:j+4) = site
		BC_list(j+5) = rgb(r,g,b)
	endif
enddo
nBC_list = k

! FDC section
nFDC_list = 0
if (NFDC > 0) then
	k = 0
    do kcell = 1,NFDC
        if (FDClist(kcell)%alive) then
			k = k+1
			j = ninfo*(k-1)
            site = FDClist(kcell)%site
            ifdcstate = 100
			FDC_list(j+1) = kcell-1
			FDC_list(j+2:j+4) = site
			FDC_list(j+5) = ifdcstate
        endif
    enddo
    nFDC_list = k
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cognate B cell colour depends on stage, state, receptor expression level.
!-----------------------------------------------------------------------------------------
subroutine BcellColour(kcell,r,g,b)
integer :: kcell, r, g, b
integer :: stage
type(cog_type),pointer :: p

p => cellist(kcell)%cptr
stage = get_stage(p)
select case(stage)
case (NAIVE)
	r = 40
	g = 140
	b = 0
case (ANTIGEN_MET, CCR7_UP)
	r = 60
	g = 230
	b = 0
case (TCELL_MET, EBI2_UP)
	r = 30
	g = 170
	b = 170
case (DIVIDING)
	if (BCL6_HI) then
		r = 250
		g = 250
		b = 40
	else
		r = 90
		g = 230
		b = 230
	endif
case (GCC_COMMIT, BCL6_UP)
	r = 250
	g = 250
	b = 40
case (PLASMA)
	r = 240
	g = 100
	b = 60
case default
	r = 255
	g = 255
	b = 255
end select
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(r, g, b)
integer :: r, g, b

rgb = ishft(r,16) + ishft(g,8) + b
end function

!-----------------------------------------------------------------------------------------
! Now this is used only to set use_TCP = .false.
! The lines
!    call get_command (b, len, status)
!    call get_command_argument (0, c, len, status)
! were failing with gfortran (don't know why), but in any case there was no need to
! get the command line arguments in this way.
!-----------------------------------------------------------------------------------------
subroutine process_command_line(ncpu,infile,outfile)
!DEC$ ATTRIBUTES DLLEXPORT :: process_command_line
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"PROCESS_COMMAND_LINE" :: process_command_line
integer :: i, cnt, len, status
integer :: ncpu
character :: c*(64), b*(256)
character*(64) :: infile,outfile
character*(64) :: progname

!write(*,*) 'process_command_line'
use_TCP = .false.   ! because this is called from para_main()							! --> use_TCP

return

ncpu = 3
infile = 'omp_para.inp'
outfile = 'omp_para.out'
!resfile = 'result.out'
!runfile = ' '

call get_command (b, len, status)
if (status .ne. 0) then
    write (logmsg,'(a,i4)') 'get_command failed with status = ', status
    call logger(logmsg)
    stop
end if
call logger('command: ')
call logger(b)
c = ''
call get_command_argument (0, c, len, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    write(*,*) c
    stop
end if
progname = c(1:len)
cnt = command_argument_count ()
if (cnt < 1) then
    write(*,*) 'Use: ',trim(progname),' num_cpu'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, len, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:len),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:len)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:len)																! --> outfile
        write(*,*) 'Output file: ',outfile
!    elseif (i == 4) then
!        resfile = c(1:len)																! --> resfile
!        write(*,*) 'Result file: ',resfile
    endif
end do

end subroutine

!-----------------------------------------------------------------------------------------
! Until I find a time function in gfortran
!----------------------------------------------------------------------------------------- 
!real(8) function timef()
!timef = 0
!end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
real(DP) :: t1, t2, tmover=0, tfields=0
integer :: k, kcell
logical :: ok
logical, save :: first = .true.

if (first) then
	first = .false.
	tmover = 0
	tfields = 0
endif
res = 0
dbug = .false.
ok = .true.
istep = istep + 1
if (dbug) then
	write(logmsg,*) 'simulate_step: ',istep
	call logger(logmsg)
endif

if (test_chemotaxis) then
	if (istep == 1) then
		do k = 1,lastcogID
			kcell = cognate_list(k)
			if (kcell > 0) then
				call set_stage(cellist(kcell)%cptr,min(k,BCL6_UP))
				call ReceptorLevel(kcell,NAIVE_TAG,GCC_TAG,1.,0.,1.,cellist(kcell)%receptor_level)
!				write(logmsg,'(a,i6,5f8.3)') 'Cognate cell receptor levels: ',kcell,cellist(kcell)%receptor_level
!				call logger(logmsg)
			endif
		enddo
	endif
!	do k = 1,lastcogID
!		kcell = cognate_list(k)
!		write(nfout,'(i4,i6,2x,3i4)') k,kcell,cellist(kcell)%site
!	enddo
	call mover(ok)
	return
endif

if (mod(istep,240) == 0) then
    if (log_traffic) then
        write(nftraffic,'(5i8,3f8.3)') istep, NBcells, Nexits, total_in, total_out, &
                InflowTotal, Vascularity
    endif
    total_in = 0
    total_out = 0
    call scanner
	call cpu_time(t1)
    call UpdateFields
	call cpu_time(t2)
	tfields = t2 - t1
!    write(*,'(a,f8.1,a,f8.1)') 'Times: mover: ',tmover,'  fields: ',tfields
    tmover = 0
!    call CheckBdryList
endif

if (dbug) write(nflog,*) 'call mover'
call cpu_time(t1)
call mover(ok)
call cpu_time(t2)
tmover = tmover + t2 - t1
if (dbug) write(nflog,*) 'did mover'
if (.not.ok) then
	call logger("mover returned error")
	res = 1
	return
endif

if (.not.evaluate_residence_time) then
	if (dbug) write(nflog,*) 'call updater'
    call updater(ok)
	if (dbug) write(nflog,*) 'did updater'
    if (.not.ok) then
		call logger('updater returned error')
		res = 1
		return
	endif
endif

if (use_traffic) then
    if (vary_vascularity) then
        call vascular
    endif
!    if (use_portal_egress) then
!		if (dbug) write(nflog,*) 'call portal_traffic'
!        call portal_traffic(ok)
!	    if (.not.ok) then
!			call logger('portal_traffic returned error')
!			res = 1
!			return
!		endif
!		if (dbug) write(nflog,*) 'did portal_traffic'
!    else
		if (dbug) write(nflog,*) 'call traffic'
        call traffic(ok)
!        if (istep > 4*240) then
!		    call CheckBdryList
!		endif
		if (dbug) write(nflog,*) 'did traffic'
	    if (.not.ok) then
			call logger('traffic returned error')
			res = 1
			return
		endif
!    endif
endif
if (dbug) call check_xyz(3)

if (dbug) write(nflog,*) 'call balancer'
call balancer(ok)
if (dbug) write(nflog,*) 'did balancer' 
if (.not.ok) then
	call logger('balancer returned error')
	res = 1
	return
endif
call set_globalvar

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif
!write(nftemp,*) 'did winsock_init'

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine


!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run.
! ncpu = the number of processors to use
! infile = file with the input data
! outfile = file to hold the output
! runfile = file to pass info to the master program (e.g. Python) as the program executes.
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: error

ok = .true.
initialized = .false.
par_zig_init = .false.
Mnodes = ncpu
inputfile = infile
outputfile = outfile
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call logger("read_Bcell_params")
call read_Bcell_params(ok)
if (.not.ok) return
call logger("did read_Bcell_params")

Fcognate = BC_COGNATE_FRACTION

ndivisions = 0
!call CD69_setparameters(K1_S1P1,K2_S1P1,K1_CD69,K2_CD69)
call array_initialisation(ok)
if (.not.ok) return
call logger('did array_initialisation')

if (calibrate_motility) then
	call motility_calibration
	stop
endif

if (log_traffic) then
    open(nftraffic,file='traffic.out',status='replace')
endif

call PlaceCells(ok)
call logger('did placeCells')
if (.not.ok) return

call CreateBdryList

NXcells = 0
!if (use_portal_egress) then
!    call placeExits
!	call adjustExits
!else
    Nexits = 0
    Lastexit = 0
!endif
check_inflow = 0
check_egress = 0
last_portal_update_time = -999

chemo_N = 8
call ChemoSetup

if (TAGGED_EXIT_CHEMOTAXIS .and. .not.use_chemotaxis) then
	call logger('ERROR: TAGGED_EXIT_CHEMOTAXIS requires USE_CHEMOTAXIS')
	ok = .false.
	return
endif
if (TAGGED_EXIT_CHEMOTAXIS .and. .not.evaluate_residence_time) then
	call logger('ERROR: TAGGED_EXIT_CHEMOTAXIS requires EVALUATE_RESIDENCE_TIME')
	ok = .false.
	return
endif
if (vary_vascularity) then
	call initialise_vascularity
endif
!write(*,*) 'NBcells, NDCalive, NTsites: ',NBcells, NDCalive, Nsites
!write(*,*) 'nlist: ',nlist

call set_globalvar
call make_split
!call scanner
call init_counters
if (save_input) then
    call save_inputfile(inputfile)
    call save_parameters
	call save_inputfile(fixedfile)
endif

call AllocateConcArrays

call ChemoSteadystate

firstSummary = .true.
initialized = .true.

write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',NBcells0
call logger(logmsg)

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine wrapup
integer :: ierr, ichemo
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
if (allocated(xoffset)) deallocate(xoffset)
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(xdomain)) deallocate(xdomain)
if (allocated(zdomain)) deallocate(zdomain)
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(Tres_dist)) deallocate(Tres_dist)
if (allocated(cellist)) deallocate(cellist,stat=ierr)
if (allocated(DClist)) deallocate(DClist,stat=ierr)
if (allocated(FDClist)) deallocate(FDClist,stat=ierr)
if (ierr /= 0) then
    write(*,*) 'cellist deallocate error: ',ierr
    stop
endif
ierr = 0
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
if (allocated(nz_sites)) deallocate(nz_sites)
if (allocated(nz_totsites)) deallocate(nz_totsites)
if (allocated(nz_cells)) deallocate(nz_cells)
if (allocated(nz_excess)) deallocate(nz_excess)
if (allocated(cognate_list)) deallocate(cognate_list)
if (allocated(life_dist)) deallocate(life_dist)
if (allocated(divide_dist)) deallocate(divide_dist)
if (allocated(exitlist)) deallocate(exitlist)
if (allocated(chemo_r)) deallocate(chemo_r)
if (allocated(chemo_p)) deallocate(chemo_p)
do ichemo = 1,MAX_CHEMO
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
enddo
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)

!if (allocated(cytp)) deallocate(cytp)
!if (allocated(xminmax)) deallocate(xminmax)
!if (allocated(inblob)) deallocate(inblob)
!if (allocated(sitelist)) deallocate(sitelist)
!if (allocated(neighbours)) deallocate(neighbours)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
inquire(nftraffic,OPENED=isopen)
if (isopen) close(nftraffic)
inquire(nfdcbind,OPENED=isopen)
if (isopen) close(nfdcbind)
inquire(nfchemo,OPENED=isopen)
if (isopen) close(nfchemo)

if (par_zig_init) then
	call par_zigfree
endif

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call SaveGenDist
if (evaluate_residence_time) then
	call write_Tres_dist
endif
call wrapup

if (res == 0) then
	call logger(' Execution successful!')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
!			call logger("closed PORT_1")
		endif
	endif
endif

end subroutine

!-----------------------------------------------------------------------------------------
! This is the DLL procedure that can be called from an external non-Fortran program to
! make a simulation run.
! Called from Python with:
!     mydll.EXECUTE(byref(ncpu),infile,n1,outfile,n2,resfile,n3,runfile,n4)
! Note that the arguments n1,n2,n3,n4, the lengths of the filename strings, are
! hidden arguments (they are not explicitly in the Fortran subroutine argument list).
! Every Python string needs to be followed by a hidden length parameter, and unless
! the declared length of a Fortran string matches the actual length of that passed from
! Python, the form character*(*) must be used.
!-----------------------------------------------------------------------------------------
subroutine execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"EXECUTE" :: execute
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success
integer :: i, res

use_CPORT1 = .false.	! DIRECT CALLING FROM C++
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

open(nflog,file='bcell-abm.log',status='replace')
awp_0%is_open = .false.
awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
!write(logmsg,*) 'resultfile: ', resfile 
!call logger(logmsg)
if (use_tcp) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
	res = 1
endif
if (test_vascular) then
	write(*,*) 'vascular_test'
	call vascular_test
	stop
endif
return
!call terminate_run(res)

end subroutine

end module

