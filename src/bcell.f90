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
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nbc0, inflow
integer :: cog_size
real :: d, rr(3), aRadius
type(Bcog_type) :: cog

ok = .false.
call logger("call rng_initialisation")
call rng_initialisation
call logger("did rng_initialisation")
!call check_rng

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the chemokine concentration gradient fields.
if (allocated(occupancy)) deallocate(occupancy)
do ichemo = 1,MAX_CHEMO
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
enddo

cog_size = (sizeof(cog) + 2)/4
nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
ngaps = 0
max_ngaps = NY*NZ
nBlist = 0
MAX_COG = 0.5*NX*NY*NZ

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
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
Radius%x = aRadius
Radius%y = aRadius/ELLIPSE_RATIO
Radius%z = aRadius

!max_nBlist = 1.5*NX*NY*NZ	! This should be computed from aRadius, Tres and ndays
! First guess NBcells0
nbc0 = (4./3.)*PI*Radius%x*Radius%y*Radius%z
! Use NBcells0 + 10*(initial influx/hour)*24*ntdays
inflow = nbc0/RESIDENCE_TIME		! cells per hour
max_nBlist = nbc0 + 10*inflow*24*days
!write(nflog,*) 'max_nBlist: ',max_nBlist
!max_nTlist = 1000	! guess
allocate(DC_list(MAX_DC))
allocate(FDC_list(MAX_FDC))
allocate(MRC_list(MAX_MRC))
allocate(occupancy(NX,NY,NZ))
allocate(Bcell_list(max_nBlist))
!allocate(Tcell_list(max_nTlist))
allocate(gaplist(max_ngaps))
allocate(cognate_list(MAX_COG))

do k = 1,max_nBlist
	nullify(Bcell_list(k)%cptr)
enddo

call make_reldir

Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
ncogseed = 0
lastcogID = 0
lastBCID = 0
!lastTCID = 0
k_nonrandom = 0
lastNBcells = 0
nadd_sites = 0
lastbalancetime = 0
localres%dN_EffCogBC = 0
localres%dN_EffCogBCGen = 0
localres%N_EffCogBC = 0
localres%N_EffCogBCGen = 0
totalres%dN_EffCogBC = 0
totalres%dN_EffCogBCGen = 0
totalres%N_EffCogBC = 0
totalres%N_EffCogBCGen = 0
totalres%N_dead = 0
totalres%dN_dead = 0

if (evaluate_residence_time) then
    allocate(Tres_dist(int(days*24)))
    Tres_dist = 0
endif

if (use_crowding_correction) then
	allocate(pressure_grad(3,NX,NY,NZ))
endif
ok = .true.

end subroutine

!--------------------------------------------------------------------------------
! Generates the arrays zoffset() and zdomain().
! The domains (slices) are numbered 0,...,2*Mnodes-1
! wz(k) = width of the slice for kth domain
! zoffset(k) = offset of kth domain occupancy array in the occupancy array.
! zdomain(x) = domain that global z lies in.
! The kth domain (slice) extends from z = zoffset(k)+1 to z = zoffset(k+1)
! The idea is to set the domain boundaries such that each domain has roughly the
! same number of available sites.
! This is the initial split, which will continue to be OK if:
! not using a blob, or Mnodes <= 2
! blobrange(:,:) holds the info about the ranges of x, y and z that the blob occupies.
! blobrange(1,1) <= x <= blobrange(1,2)
! blobrange(2,1) <= y <= blobrange(2,2)
! blobrange(3,1) <= z <= blobrange(3,2)
!--------------------------------------------------------------------------------
subroutine make_split(force)
logical :: force
integer :: k, wsum, kdomain, nsum, Ntot, N, last, x, y, z
integer, allocatable :: scount(:)
integer, allocatable :: wz(:), ztotal(:)
integer :: Mslices
real :: dNT, diff1, diff2
logical :: show = .false.

!write(*,*) 'make_split: istep,Mnodes: ',istep,Mnodes
if (Mnodes == 1) then
    Mslices = 1
    zdomain = 0
else
	Mslices = 2*Mnodes
endif
dNT = abs(NBcells - lastNBcells)/real(lastNBcells+1)
if (.not.force .and. dNT < 0.03) then
    return
endif
lastNBcells = NBcells
if (Mslices > 1) then
	allocate(wz(0:Mslices))
	allocate(ztotal(0:Mslices))
	allocate(scount(NX))
endif
blobrange(:,1) = 99999
blobrange(:,2) = 0
nsum = 0
do z = 1,NZ
    k = 0
    do y = 1,NY
        do x = 1,NX
            if (occupancy(x,y,z)%indx(1) /= OUTSIDE_TAG) then
                k = k + 1
                blobrange(1,1) = min(blobrange(1,1),x)
                blobrange(1,2) = max(blobrange(1,2),x)
                blobrange(2,1) = min(blobrange(2,1),y)
                blobrange(2,2) = max(blobrange(2,2),y)
                blobrange(3,1) = min(blobrange(3,1),z)
                blobrange(3,2) = max(blobrange(3,2),z)
            endif
        enddo
    enddo
    if (Mslices > 1) then
	    scount(z) = k
	    nsum = nsum + scount(z)
	endif
enddo
if (Mslices == 1) return

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
! The blob centre is at (x0,y0,z0), and the blob radius is R = Radius%x
! NOT USED
!--------------------------------------------------------------------------------
integer function slice_count(x)
integer :: x
real :: r2

r2 = Radius%x**2 - (x-x0)**2
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
type(Bcell_type) :: cell
logical :: ok

write(*,*) 'motility_calibration'

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
        call make_split(.true.)
        if (nBlist > 0) then
	        write(*,*) 'make tag list: NBcells,nBlist,ntagged: ',NBcells,nBlist,ntagged

            allocate(tagseq(NBcells))
            allocate(tagid(ntagged))
            allocate(tagsite(3,ntagged,0:nvar))
            tagseq = 0
            k = 0
	        kpath = 0
            do ic = 1,nBlist
                if (Bcell_list(ic)%ctype == TAGGED_CELL) then
                    id = Bcell_list(ic)%ID
                    k = k+1
                    tagid(k) = id
                    tagseq(id) = k
                    tagsite(:,k,0) = Bcell_list(ic)%site
					if (motility_save_paths) then
						if (kpath < npaths) then
							kpath = kpath + 1
							pathcell(kpath) = ic
						endif
					endif
                endif
            enddo
        endif

        ns = min(ns,nBlist)
        allocate(prevsite(3,ns))
        do ic = 1,ns
            prevsite(:,ic) = Bcell_list(ic)%site
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
!                write(*,*) 'Cell 1: ',Bcell_list(1)%site
                do ic = 1,ns
                    ds = Bcell_list(ic)%site - prevsite(:,ic)
                    prevsite(:,ic) = Bcell_list(ic)%site
                    d = sqrt(real(ds(1)*ds(1) + ds(2)*ds(2) + ds(3)*ds(3)))
                    ssum = ssum + d*DELTA_X/DELTA_T
                enddo
                if (motility_save_paths) then
                    k = (imin-1)*nsteps_per_min + isub
                    if (k >= nvar0*nsteps_per_min .and. k < nvar0*nsteps_per_min + npos) then
                        write(nfpath,'(160i4)') (Bcell_list(pathcell(kpath))%site(1:2),kpath=1,npaths)
                    endif
                endif
            enddo
            write(*,*) 'speed: ',ssum/(ns*nsteps_per_min*imin)

            do ic = 1,nBlist
                cell = Bcell_list(ic)
                if (cell%ctype == TAGGED_CELL) then
                    id = cell%ID
                    k = tagseq(id)
                    tagsite(:,k,imin) = cell%site
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
integer :: x, y, z, istep
integer, parameter :: Ntimes = 50   !500
real(DP) :: t1, t2, Cnumeric(Ntimes),Canalytic(Ntimes)
logical :: ok

write(*,*) 'diffusion_calibration: '

x = NX/4
!call analytical_soln(x,Canalytic,Ntimes)
call placeCells(ok)
if (.not.ok) stop
call make_split(.true.)

call cpu_time(t1)
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            occupancy(x,y,z)%indx = 1
        enddo
    enddo
enddo

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
! To test addBTcell()
!-----------------------------------------------------------------------------------------
subroutine add_random_cells(n, ctype, gen, stage, status, region)
integer :: n, ctype, gen, stage, status, region
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
        call addBTcell(site,ctype,gen,stage,status,region,kcell,ok)
        if (dbug) write(*,'(a,7i6)') 'after add_random_cells: ',site,occupancy(x,y,z)%indx,slots
        call checkslots('add_random_cells: ',site)
        k = k+1
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
integer :: k, x, y, z, indx(2), site(3), kcell, ctype, gen, status, region, kpar=0
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
            ctype = select_Bcell_type(kpar)
        endif
        if (ctype /= NONCOG_TYPE_TAG) then
            ncogseed = ncogseed + 1
        endif
        status = BCL6_LO
        call addBTcell(site,ctype,gen,NAIVE,status,region,kcell,ok)
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
integer :: x, y, z, k, kcell, indx(2), ctype, gen, stage, region, site(3), n, slot
integer :: node_inflow, node_outflow, add(3), net_inflow, ihr
real(DP) :: R, df, prob
real :: tnow, noncog_exit_prob
logical :: can_leave, left
integer :: kpar=0
integer :: nb, nc, nx
type(boundary_type), pointer :: bdry
type(Bcog_type), pointer :: p

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
if (.not.ok) then
	write(*,*) 'cellInflux error'
	return
endif
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
				p => Bcell_list(kcell)%cptr
				if (associated(p)) then
					stage = get_stage(p)
					if (stage == PLASMA) then
						can_leave = .true.
					endif
				elseif (Bcell_list(kcell)%ctype /= COG_CD4_CELL) then	! do not allow CD4 T cells to leave (for now)
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
						if (Bcell_list(kcell)%ctype == RES_TAGGED_CELL) then
							noutflow_tag = noutflow_tag + 1
							restime_tot = restime_tot + tnow - Bcell_list(kcell)%entrytime
							ihr = (tnow - Bcell_list(kcell)%entrytime)/60. + 1
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
if (dbug) then
	write(nfout,'(a,5i6,f8.4)') 'traffic: ',istep,nbcells,node_inflow,nb,nx,noncog_exit_prob
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
        if (Bcell_list(kcell)%ctype == RES_TAGGED_CELL) then
            noutflow_tag = noutflow_tag + 1
            restime_tot = restime_tot + tnow - Bcell_list(kcell)%entrytime
            ihr = (tnow - Bcell_list(kcell)%entrytime)/60. + 1
            Tres_dist(ihr) = Tres_dist(ihr) + 1
        endif
    endif
enddo
nadd_sites = nadd_sites - node_outflow
!write(*,'(a,4i6)') 'istep,in,out: ',istep,node_inflow,node_outflow,NBcells
end subroutine

!-----------------------------------------------------------------------------------------
! Determine whether a cognate T cell is licensed to exit the DCU.
!-----------------------------------------------------------------------------------------
logical function exitOK(p,ctype)
type(Bcog_type), pointer :: p
integer :: ctype

exitOK = .true.
if (ctype == COG_TYPE_TAG) then
    exitOK = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! The cell kcell in slot at esite(:) is a candidate to exit.  If it meets the criteria,
! left = .true.
! Currently, when a T cell exits the LN its ID is set to 0.
! If use_gaplist is true, a gap is recorded in gaplist(:).
!-----------------------------------------------------------------------------------------
subroutine CellExit(kcell,slot,esite)
integer :: kcell, slot, esite(3)
integer :: x, y, z, ctype, gen, stage, region, status
real :: tnow
logical :: cognate, activated
type(Bcog_type), pointer :: p

tnow = istep*DELTA_T
if (evaluate_residence_time) then
    cognate = .false.
elseif (associated(Bcell_list(kcell)%cptr)) then
    cognate = .true.
    p => Bcell_list(kcell)%cptr
	stage = get_stage(p)
	status = get_status(p)
    gen = get_generation(p)
    ctype = Bcell_list(kcell)%ctype
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
    call set_stage(p,LEFT)
	call set_region(p,GONE)
	write(nflog,*) 'cell left: ',kcell,Bcell_list(kcell)%ID,stage
endif
Bcell_list(kcell)%exists = .false.
!write(nflog,*) 'cell left: ',kcell,Bcell_list(kcell)%ID
if (.not.cognate) then
	Bcell_list(kcell)%ID = 0
endif
if (use_gaplist) then
	ngaps = ngaps + 1
	if (ngaps > max_ngaps) then
		call logger('Error: gaplist dimension exceeded')
		stop
	endif
	gaplist(ngaps) = kcell
endif
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
	VEGF_baserate = VEGF_beta*NBcells0
    VEGFmass = VEGF_baserate/VEGF_decayrate    ! steady-state VEGF level M_G0
    Cvegf0 = VEGFmass/NBcells0	! taking K_V = 1.0
    vasc_decayrate = vasc_maxrate*hill(Cvegf0,vasc_beta*Cvegf0,vasc_n)	! delta_V
endif
Vascularity = 1.00
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
integer :: k, idc, naddDC, naddex, nremex, dexit
real :: tnow
integer :: kpar = 0
logical :: blob_changed

ok = .true.
if (dbug) write(*,*) 'balancer: ',istep

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
    if (blob_changed) then
		if (dbug) write(nflog,*) 'call make_split'
        call make_split(.false.)
    endif
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
type(Bcog_type), pointer :: p
integer :: ctype, gen, i

gen = get_generation(p)
if (ctype > NCTYPES) then
    write(*,*) 'efferent: bad cell type:', ctype
    stop
endif
totalres%dN_EffCogBC(ctype)  = totalres%dN_EffCogBC(ctype) + 1
totalres%dN_EffCogBCGen(gen) = totalres%dN_EffCogBCGen(gen) + 1
totalres%N_EffCogBC(ctype)   = totalres%N_EffCogBC(ctype) + 1
totalres%N_EffCogBCGen(gen)  = totalres%N_EffCogBCGen(gen) + 1
end subroutine

!-----------------------------------------------------------------------------------------
! Using the complete list of cells, Bcell_list(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine show_snapshot(ok)
logical :: ok
integer :: kcell, ctype, ncog, noncog, ncd4, ntot, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead
real :: stim(FINISHED), tgen, tnow, fac, act, cyt_conc, mols_pM
type(Bcog_type), pointer :: p
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
ncd4 = 0
ntot = 0
nst = 0
stim = 0
gendist = 0
div_gendist = 0
do kcell = 1,nBlist
    if (.not.Bcell_list(kcell)%exists) cycle
    p => Bcell_list(kcell)%cptr
    ntot = ntot + 1
    ctype = Bcell_list(kcell)%ctype
    if (ctype == COG_TYPE_TAG) then
        ncog = ncog + 1
        stage = get_stage(p)
        nst(stage) = nst(stage) + 1
        gen = get_generation(p)
        if (gen < 0 .or. gen > BC_MAX_GEN) then
            write(logmsg,'(a,2i6)') 'show_snapshot: bad gen: ',kcell,gen
            call logger(logmsg)
            ok = .false.
            return
        endif
        gendist(gen) = gendist(gen) + 1
    elseif (ctype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
	elseif (ctype == COG_CD4_CELL) then
		ncd4 = ncd4 + 1
    else
        write(*,*) 'ERROR: show_snapshot: bad ctype: ',ctype
        stop
    endif
enddo
do i = 1,FINISHED
    if (nst(i) > 0) then
        stim(i) = stim(i)/nst(i)
    else
        stim(i) = 0
    endif
enddo
tgen = sum(gendist)
do i = BC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogBCGen(1:BC_MAX_GEN))
do i = BC_MAX_GEN,1,-1
    if (totalres%N_EffCogBCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead

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
write(*,'(a,2i8,4x,i8)') 'Recent efferent: ',totalres%dN_EffCogBC(2:3),sum(totalres%dN_EffCogBC)
write(*,'(a,2i8,4x,i8)') 'Total efferent:  ',totalres%N_EffCogBC(2:3),teffgen
write(*,'(a,10i6)')   'gen dist: ',(i,i=1,10)
write(*,'(a)')        'In node:  '
write(*,'(10x,10f6.3)') gendist(1:ngens)/tgen
write(*,'(a)')        'Efferent: '
write(*,'(10x,10f6.3)') fac*totalres%N_EffCogBCGen(1:neffgens)
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

! To plot outflow variation with time
deallocate(gendist)
totalres%dN_EffCogBC = 0
totalres%dN_EffCogBCGen = 0
totalres%dN_Dead = 0

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_stim_dist
integer :: kcell, ctype, ncog, noncog, ncd4
real :: s
type(Bcog_type), pointer :: p

ncog = 0
noncog = 0
ncd4 = 0
do kcell = 1,nBlist
    if (.not.Bcell_list(kcell)%exists) cycle
    p => Bcell_list(kcell)%cptr
    ctype = Bcell_list(kcell)%ctype
    if (ctype == COG_TYPE_TAG) then
        ncog = ncog + 1
    elseif (ctype == NONCOG_TYPE_TAG) then
        noncog = noncog + 1
    elseif (ctype == COG_CD4_CELL) then
		ncd4 = ncd4 + 1
    else
        write(*,*) 'ERROR: compute_stim_dist: bad ctype: ',ctype
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
do kcell = 1,nBlist
    if (.not.Bcell_list(kcell)%exists) cycle
    z = Bcell_list(kcell)%site(3)
    cognate = (associated(Bcell_list(kcell)%cptr))
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
    j = select_Bcell_type(kpar)
    if (j == COG_CD4_CELL) ncog = ncog + 1
    if (mod(k,1000) == 0) then
        write(*,*) ncog,real(ncog)/k
    endif
enddo
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
write(nfout,'(a,f6.3)') 'chemo_K_exit:            ',chemo_K_exit
write(nfout,'(a,f6.3)') '  K1_S1P1:               ',K1_S1P1
write(nfout,'(a,3i8)') 'Results: noutflow_tag, ninflow_tag,kmax: ',noutflow_tag,ninflow_tag,kmax
write(nfout,'(a,f6.1)') 'Residence time: ',Tres
write(nfout,'(a)') ' Hour   Probability'
do k = 1,kmax
    write(nfout,'(f6.1,e12.4)') k-0.5,Tres_dist(k)/noutflow_tag
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The distributions gendist and tcrdist are all for the current DCU population.
!-----------------------------------------------------------------------------------------
subroutine write_results
integer :: kcell, ctype, ncog, ntot, i
integer :: gen
real :: hour
type(Bcog_type), pointer :: p
integer :: gendist(BC_MAX_GEN)
character*(60) :: fmtstr = '(f6.2,2i8,4x,15f7.4,4x,10f7.4,4x,10f7.4,4x,10i7)'

write(fmtstr(14:15),'(i2)') BC_MAX_GEN
hour = istep*DELTA_T/60
gendist = 0
ntot = 0
ncog = 0
do kcell = 1,nBlist
    if (.not.Bcell_list(kcell)%exists) cycle
    p => Bcell_list(kcell)%cptr
    ntot = ntot + 1
    ctype = Bcell_list(kcell)%ctype
    if (ctype == COG_TYPE_TAG) then
        ncog = ncog + 1
        gen = get_generation(p)
        gendist(gen) = gendist(gen) + 1
    endif
enddo
write(nfres,fmtstr) hour,ntot,ncog,real(gendist)/ncog
end subroutine

!-----------------------------------------------------------------------------------------
! A different approach to tagging cells for residence time computation
!-----------------------------------------------------------------------------------------
subroutine tag_cells
integer :: kcell
type(Bcell_type),pointer :: cell

do kcell = 1,nBlist
    cell => Bcell_list(kcell)
	if (.not.cell%exists) cycle
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
		gen = get_generation(Bcell_list(kcell)%cptr)
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
! Using the complete list of cells, Bcell_list(), extract info about the current state of the
! paracortex.  This info must be supplemented by counts of cells that have died and cells that
! have returned to the circulation.
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)
logical :: ok
integer :: kcell, ctype, ncog, noncog, ncd4, ntot, nbnd, stage, region, i, iseq, error
integer :: gen, ngens, neffgens, teffgen, dNdead, Ndead, nvasc
real :: stim(2*STAGELIMIT), tgen, tnow, fac, act, cyt_conc, mols_pM
type(Bcog_type), pointer :: p
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
ncd4 = 0
ntot = 0
nbnd = 0
nst = 0
stim = 0
gendist = 0
div_gendist = 0
do kcell = 1,nBlist
    if (.not.Bcell_list(kcell)%exists) cycle
    p => Bcell_list(kcell)%cptr
    if (associated(p)) then
		stage = get_stage(p)
		region = get_region(p)
	else
		stage = 0
		region = FOLLICLE
	endif
	if (region == FOLLICLE) then
	    ntot = ntot + 1
	else
		write(logmsg,*) 'Error: get_summary: kcell, region: ',kcell,region
		call logger(logmsg)
		stop
	endif
    ctype = Bcell_list(kcell)%ctype
    if (ctype == COG_TYPE_TAG) then
        ncog = ncog + 1
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
    elseif (ctype == COG_CD4_CELL) then
		ncd4 = ncd4 + 1
    else
        write(logmsg,*) 'ERROR: get_summary: bad ctype: ',ctype
        call logger(logmsg)
        stop
    endif
enddo

tgen = sum(gendist)
do i = BC_MAX_GEN,1,-1
    if (gendist(i) /= 0) exit
enddo
ngens = i

teffgen = sum(totalres%N_EffCogBCGen(1:BC_MAX_GEN))
do i = BC_MAX_GEN,1,-1
    if (totalres%N_EffCogBCGen(i) /= 0) exit
enddo
neffgens = i

dNdead = totalres%dN_Dead
Ndead = totalres%N_Dead

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
endif
deallocate(gendist)
deallocate(div_gendist)

totalres%dN_EffCogBC = 0
totalres%dN_EffCogBCGen = 0
totalres%dN_Dead = 0

summaryData(1:9) = (/int(tnow/60),istep,ntot,ncogseed,ncog,Ndead,int(InflowTotal*60/DELTA_T), int(100*vascularity), teffgen/)
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
! Improving this:
! blobrange(:,:) holds the info about the ranges of x, y and z that the blob occupies.
! blobrange(1,1) <= x <= blobrange(1,2)
! blobrange(2,1) <= y <= blobrange(2,2)
! blobrange(3,1) <= z <= blobrange(3,2)
!--------------------------------------------------------------------------------
subroutine get_scene(nBC_list,BC_list,nFDCMRC_list,FDCMRC_list,nbond_list,bond_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nFDCMRC_list, nBC_list, nbond_list, FDCMRC_list(*), BC_list(*), bond_list(*)
integer :: k, kc, kcell, site(3), j, jb, idc, fdcsite(3)
integer :: col(3)
integer :: x, y, z
real :: dcstate
integer :: ifdcstate, ibcstate, ctype, stage, region
real :: bcell_diam = 0.9
real :: FDC_diam = 1.8
integer :: gen, bnd(2), noncnt, last_id1, last_id2
logical :: show_Tcells = .true.
logical :: ok
integer, parameter :: axis_centre = -2	! identifies the ellipsoid centre
integer, parameter :: axis_end    = -3	! identifies the ellipsoid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the ellipsoid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 5			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the follicle

nFDCMRC_list = 0
! Need some markers to delineate the follicle extent.  These nax "cells" are used to convey (the follicle centre
! and) the approximate ellipsoidal blob limits in the 3 axis directions.
do k = 1,nax
	select case (k)
!	case (1)
!		x = Centre(1) + 0.5
!		y = Centre(2) + 0.5
!		z = Centre(3) + 0.5
!		site = (/x, y, z/)
!		ibcstate = axis_centre
	case (1)
!		x = Centre(1) - Radius%x - 2
		x = blobrange(1,1) - 1
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (2)
!		x = Centre(1) + Radius%x + 2
		x = blobrange(1,2) + 1
		y = Centre(2) + 0.5
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (3)
		x = Centre(1) + 0.5
!		y = Centre(2) - Radius%y - 2
		y = blobrange(2,1) - 1
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_bottom
	case (4)
		x = Centre(1) + 0.5
!		y = Centre(2) + Radius%y + 2
		y = blobrange(2,2) + 1
		z = Centre(3) + 0.5
		site = (/x, y, z/)
		ibcstate = axis_end
	case (5)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
!		z = Centre(3) - Radius%z - 2
		z = blobrange(3,1) - 1
		site = (/x, y, z/)
		ibcstate = axis_end
	case (6)
		x = Centre(1) + 0.5
		y = Centre(2) + 0.5
!		z = Centre(3) + Radius%z + 2
		z = blobrange(3,2) + 1
		site = (/x, y, z/)
		ibcstate = axis_end
	end select

	j = ninfo*(k-1)
	BC_list(j+1) = k-1
	BC_list(j+2:j+4) = site
	BC_list(j+5) = ibcstate
	last_id1 = k-1
enddo
k = last_id1 + 1

! B cell section

! Cognate cells
do kc = 1,lastcogID
	kcell = cognate_list(kc)
	if (Bcell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = Bcell_list(kcell)%site
		call BcellColour(kcell,col)
		BC_list(j+1) = kc + last_id1
		BC_list(j+2:j+4) = site
		BC_list(j+5) = rgb(col)
		last_id2 = kc + last_id1
!		write(nflog,'(a,8i7)') 'cog: ',k,BC_list(j+1),kcell,get_stage(Bcell_list(kcell)%cptr),kc,last_id1
!		if (istep > 14000) then
!			write(logmsg,'(a,7i8)') 'cog: ',k,kcell,col,get_stage(Bcell_list(kcell)%cptr)
!			call logger(logmsg)
!		endif
	endif
enddo
noncnt = 0
if (use_Tcells .and. show_Tcells) then
    ! CD4 T cells are displayed
    do kcell = 1,nBlist
        if (associated(Bcell_list(kcell)%cptr)) cycle
        if (Bcell_list(kcell)%ctype /= COG_CD4_CELL) cycle
	    if (Bcell_list(kcell)%exists) then
	        noncnt = noncnt + 1
		    k = k+1
		    j = ninfo*(k-1)
		    site = Bcell_list(kcell)%site
		    call BcellColour(kcell,col)
		    BC_list(j+1) = noncnt + last_id2
		    BC_list(j+2:j+4) = site
		    BC_list(j+5) = rgb(col)
!    		write(nflog,*) 'CD4: ',k,kcell,BC_list(j+1)
!		    write(logmsg,'(a,7i8)') 'CD4: ',noncnt,k,kcell,col,rgb(col)
!		    call logger(logmsg)
        else
            write(nflog,*) 'cell nonexistent: ',kcell
	    endif
    enddo
endif
noncnt = 0
if (show_noncognate) then
    ! Non-cognate cells, a specified fraction noncogfraction are displayed
    do kcell = 1,noncog_display_fraction*nBlist
        if (associated(Bcell_list(kcell)%cptr)) then
            write(nflog,*) 'cell associated: ',kcell
            cycle
        endif
        if (Bcell_list(kcell)%ctype == COG_CD4_CELL) cycle
	    if (Bcell_list(kcell)%exists) then
	        noncnt = noncnt + 1
		    k = k+1
		    j = ninfo*(k-1)
		    site = Bcell_list(kcell)%site
		    call BcellColour(kcell,col)
!		    BC_list(j+1) = kcell + last_id2
		    BC_list(j+1) = noncnt + last_id2
		    BC_list(j+2:j+4) = site
		    BC_list(j+5) = rgb(col)
!    		write(nflog,*) 'non: ',k,kcell,BC_list(j+1)
!		    write(logmsg,'(a,7i8)') 'noncog: ',noncnt,k,kcell,col,rgb(col)
!		    call logger(logmsg)
        else
            write(nflog,*) 'cell nonexistent: ',kcell
	    endif
    enddo
endif
nBC_list = k
!write(logmsg,*) 'nBC_list: ',nBC_list
!call logger(logmsg)
! FDC section
k = 0
if (NFDC > 0) then
    do kcell = 1,NFDC
        if (FDC_list(kcell)%alive) then
			k = k+1
			j = ninfo*(k-1)
            site = FDC_list(kcell)%site
!            write(logmsg,*) 'FDC site: ',kcell,site
!            call logger(logmsg)
			FDCMRC_list(j+1) = kcell-1
			FDCMRC_list(j+2:j+4) = site
			FDCMRC_list(j+5) = 100
        endif
    enddo
endif
if (NMRC > 0) then
    do kcell = 1,NMRC
        if (MRC_list(kcell)%alive) then
			k = k+1
			j = ninfo*(k-1)
            site = MRC_list(kcell)%site
!            write(logmsg,*) 'MRC site: ',kcell,site
!            call logger(logmsg)
			FDCMRC_list(j+1) = NFDC + kcell-1
			FDCMRC_list(j+2:j+4) = site
			FDCMRC_list(j+5) = 200
        endif
    enddo
endif
nFDCMRC_list = k
!write(logmsg,*) 'istep, nFDC_list: ',istep,nFDCMRC_list
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cognate B cell colour depends on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine BcellColour(kcell,col)
integer :: kcell, col(3)
integer :: stage, status
type(Bcog_type), pointer :: p
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

p => Bcell_list(kcell)%cptr
if (associated(p)) then
    stage = get_stage(p)
    status = get_status(p)
    select case(stage)
    case (NAIVE)
	    col = DEEPBLUE			! 1
    case (ANTIGEN_MET, CCR7_UP)
	    col = LIGHTBLUE			! 2
    case (TCELL_MET, EBI2_UP)
	    col = LIGHTGREEN		! 3
    case (DIVIDING) 
	    if (status == BCL6_HI) then
		    col = YELLOW		! 4
	    else
		    col = DEEPGREEN		! 5
	    endif
    case (GCC_COMMIT, BCL6_UP)
	    col = YELLOW			! 4
    case (PLASMA)
	    col = RED		! 6
    case (FINISHED)
	    col = GRAY				! 7
    case default
	    col = WHITE				! 8
    end select
elseif (Bcell_list(kcell)%ctype == COG_CD4_CELL) then
	col = LIGHTORANGE
else
	col = WHITE				    ! non-cognate
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine squeezer_test
integer :: i, j, kcell, site(3), indx(2), slot, region, nb0, n, kpar=0
type(Bcog_type), pointer :: p
logical :: ok, flag

write(*,*) 'doing squeezer_test'
flag = .false.
nb0 = nbcells
do j = 1,5
write(*,*) 'nBlist: ',nBlist, NBcells
do i = 1,200
	kcell = random_int(1,nBlist,kpar)
	if (.not.Bcell_list(kcell)%exists) cycle
	site = Bcell_list(kcell)%site
	indx = occupancy(site(1),site(2),site(3))%indx
	if (indx(1) == kcell) then
		slot = 1
	else
		slot = 2
	endif
!	call CellExit(kcell,slot,site)
	call BcellDeath(kcell)
	if (Bcell_list(kcell)%ID == 12391) then
		write(*,*) 'Dead cell ID=12391: ',kcell
	endif
enddo
write(*,*)
call squeezer(.true.)
write(*,*) 'did squeezer'
write(*,*)
n = nb0-nbcells
call cellInflux(n,ok)
write(*,*) 'did cellInflux'
do kcell = 1,nBlist
    p => Bcell_list(kcell)%cptr
	if (associated(p)) then
        region = get_region(p)
        if (Bcell_list(kcell)%exists) then
			if (region /= FOLLICLE) then
				write(*,*) 'Cell exists, bad region: ',kcell,Bcell_list(kcell)%ID,region
				stop
			endif
		else
			if (region == FOLLICLE) then
				write(*,*) 'Cell does not exist, bad region: ',kcell,Bcell_list(kcell)%ID,region
				stop
			endif
		endif
	endif
enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine disableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from bcell_main()	
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one time step (DELTA_T)
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
real(DP) :: t1, t2, tmover=0, tfields=0
integer :: k, kcell, stage, kpar = 0
logical :: ok, sim_dbug
logical, save :: first = .true.

if (first) then
	first = .false.
	tmover = 0
	tfields = 0
endif
res = 0
ok = .true.
istep = istep + 1
sim_dbug = .false.
dbug = .false.
if (istep < 0) then
	dbug = .true.
	sim_dbug = .true.
endif
if (sim_dbug) then
	write(logmsg,*) 'simulate_step: ',istep
	call logger(logmsg)
	call checkcellcount(ok)
	if (.not.ok) stop
endif

if (sim_dbug) then
	kcell = cognate_list(1)
!	write(*,*) Bcell_list(kcell)%site
	call check_cognate_list
endif

if (test_chemotaxis) then
	if (istep == 1) then
		do k = 1,lastcogID
			kcell = cognate_list(k)
			if (kcell > 0) then
				call set_stage(Bcell_list(kcell)%cptr,PLASMA)
				call set_status(Bcell_list(kcell)%cptr,BCL6_LO)
				call ReceptorLevel(kcell,NAIVE_TAG,ACTIVATED_TAG,1.,0.,1.,Bcell_list(kcell)%receptor_level)
			endif
		enddo
	endif
!	do k = 1,lastcogID
!		kcell = cognate_list(k)
!		stage = get_stage(Bcell_list(kcell)%cptr)
!		write(nfout,'(i4,i6,2x,3i4)') k,kcell,stage
!	enddo
	call mover(ok)
	return
endif

if (test_case(1)) then    ! cognate cells are made insensitive to all chemokines (note: other NAIVE cells sense all except CXCL13)
    if (istep == 1) then
		do k = 1,lastcogID
			kcell = cognate_list(k)
			if (kcell > 0) then
!				Bcell_list(kcell)%ctype = TESTCELL2
				Bcell_list(kcell)%receptor_level = 0
!				if (k <= 5) then
!    				call set_stage(Bcell_list(kcell)%cptr,PLASMA)
!    		    endif
            endif
        enddo
    endif
	call mover(ok)
	return
endif
if (test_case(2)) then    
	! Testing the effects of the crowding correction
	! Use only S1P, make non-cognate cells highly attracted, cognate cells are PLASMA, no chemotaxis
	! (Note: chemokine parameters must be set appropriately).
    if (istep == 1) then
		do k = 1,lastcogID
			kcell = cognate_list(k)
			if (kcell > 0) then
				Bcell_list(kcell)%receptor_level = 0
   				call set_stage(Bcell_list(kcell)%cptr,PLASMA)
            endif
        enddo
    endif
!    if (mod(istep,10) == 1) then
		call pressure_gradient
!	endif
	call mover(ok)
	return
endif

if (mod(istep,240) == 0) then
    if (log_traffic) then
        write(nftraffic,'(4i8,3f8.3)') istep, NBcells, total_in, total_out, &
                InflowTotal, Vascularity
    endif
    total_in = 0
    total_out = 0
!	call cpu_time(t1)
	if (use_SS_fields) then
		call make_split(.true.)		! just to be sure
		if (sim_dbug) write(nflog,*) 'call UpdateSSFields'
		call UpdateSSFields
		if (sim_dbug) write(nflog,*) 'did UpdateSSFields'
	endif
!	call cpu_time(t2)
!	tfields = t2 - t1
!    write(*,'(a,f8.1,a,f8.1)') 'Times: mover: ',tmover,'  fields: ',tfields
    tmover = 0
!    call CheckBdryList
!	call show_lineage(logID)
!	call checkreceptor()
	call checkcellcount(ok)
	if (.not.ok) then
		res = 1
		return
	endif
endif

if (.not.use_SS_fields) then
	if (sim_dbug) write(nflog,*) 'call UpdateFields'
	call UpdateFields(DELTA_T)
	if (sim_dbug) write(nflog,*) 'did UpdateFields'
endif
if (sim_dbug) then
	write(nflog,*) 'call mover'
	call checkcellcount(ok)
endif
call cpu_time(t1)
call mover(ok)
call cpu_time(t2)
tmover = tmover + t2 - t1
if (sim_dbug) then
	write(nflog,*) 'did mover'
	call checkcellcount(ok)
endif
if (.not.ok) then
	call logger("mover returned error")
	res = 1
	return
endif

if (.not.evaluate_residence_time) then
	if (sim_dbug) write(nflog,*) 'call updater'
    call updater(ok)
	if (sim_dbug) write(nflog,*) 'did updater'
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
	if (sim_dbug) write(nflog,*) 'call traffic'
    call traffic(ok)
	if (sim_dbug) write(nflog,*) 'did traffic'
    if (.not.ok) then
		call logger('traffic returned error')
		res = 1
		return
	endif
endif
if (sim_dbug) call check_xyz(3)

call balancer(ok)
if (sim_dbug) then
	write(nflog,*) 'did balancer' 
	if (.not.ok) then
		res = 1
		return
	endif
endif
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

chemo_N = 8
call ChemoSetup

if (vary_vascularity) then
	call initialise_vascularity
endif

call set_globalvar
call make_split(.true.)
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

call checkcellcount(ok)
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
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(zdomain)) deallocate(zdomain)
!if (allocated(occupancy)) deallocate(occupancy)
if (allocated(Tres_dist)) deallocate(Tres_dist)
if (allocated(Bcell_list)) deallocate(Bcell_list,stat=ierr)
!if (allocated(Tcell_list)) deallocate(Tcell_list,stat=ierr)
if (allocated(DC_list)) deallocate(DC_list,stat=ierr)
if (allocated(FDC_list)) deallocate(FDC_list,stat=ierr)
if (allocated(MRC_list)) deallocate(MRC_list,stat=ierr)
if (ierr /= 0) then
    write(logmsg,*) 'deallocate error: ',ierr
    call logger(logmsg)
endif
ierr = 0
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)
if (allocated(cognate_list)) deallocate(cognate_list)
if (allocated(life_dist)) deallocate(life_dist)
if (allocated(divide_dist)) deallocate(divide_dist)
if (allocated(chemo_r)) deallocate(chemo_r)
if (allocated(chemo_p)) deallocate(chemo_p)
!do ichemo = 1,MAX_CHEMO
!	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
!	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
!	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
!enddo
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)

if(allocated(pressure_grad)) deallocate(pressure_grad)
call logger('deallocated all arrays')

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
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
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
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
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
if (test_squeezer) then
	write(*,*) 'squeezer_test'
	call squeezer_test
	stop
endif
end subroutine

end module

