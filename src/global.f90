!==========================================================================================
! Modification notes
! ------------------

!==========================================================================================

module global

use omp_lib
use par_zig_mod
use winsock

implicit none

INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
integer, parameter :: REAL_KIND = 4

! General parameters
integer, parameter :: BIG_INT = 2**30
integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MODEL = MOORE26_MODEL
integer, parameter :: MAXRELDIR = 26

! B cell activation stage
integer, parameter :: NAIVE		  = 1
integer, parameter :: ANTIGEN_MET = 2
integer, parameter :: CCR7_UP     = 3
integer, parameter :: TCELL_MET   = 4
integer, parameter :: EBI2_UP     = 5
integer, parameter :: DIVIDING    = 6
integer, parameter :: GCC_COMMIT  = 7
integer, parameter :: PLASMA      = 8
integer, parameter :: BCL6_UP     = 9
integer, parameter :: DEAD        = 10
integer, parameter :: LEFT        = 11
integer, parameter :: FINISHED    = 12
integer, parameter :: STAGELIMIT  = 12

integer, parameter :: BCL6_LO = 1
integer, parameter :: BCL6_HI = GCC_COMMIT

real, parameter :: T_CCR7_UP   = 2*60
real, parameter :: T_EBI2_UP = 2*60
real, parameter :: T_BCL6_UP   = 2*60
real, parameter :: T_FIRST_DIVISION = 12*60
real, parameter :: T_DIVISION       = 10*60

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging) 
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)

integer, parameter :: SLOT_NUM1 = 1
integer, parameter :: SLOT_NUM2 = 2
integer, parameter :: BOTH = 3

integer, parameter :: MAX_CHEMO = 4
integer, parameter :: S1P    = 1
integer, parameter :: CCL21  = 2
integer, parameter :: OXY    = 3
integer, parameter :: CXCL13 = 4
integer, parameter :: MAX_RECEPTOR = 5
integer, parameter :: S1PR1    = 1		! S1P
integer, parameter :: CCR7     = 2		! CCL21
integer, parameter :: EBI2     = 3		! OXY
integer, parameter :: CXCR5    = 4		! CXCL13
integer, parameter :: S1PR2    = 5		! S1P (negative)

real, parameter :: receptor_level(5,4) = reshape((/ 1.,1.,1.,0.,0., .2,2.,2.,0.,0., .2,.5,2.,1.,0., 0.,0.,.2,2.,2. /), (/5,4/))	! Taka

integer, parameter :: NCTYPES = 4
integer, parameter :: NONCOG_TYPE_TAG  = 1
integer, parameter :: COG_TYPE_TAG  = 2
integer, parameter :: COG_CD4_TAG  = 2
integer, parameter :: COG_CD8_TAG  = 3
integer, parameter :: TAGGED_CELL = 100
integer, parameter :: RES_TAGGED_CELL = 101
integer, parameter :: OUTSIDE_TAG = -9999

integer, parameter :: NAIVE_TAG = 1
integer, parameter :: ANTIGEN_TAG = 2
integer, parameter :: ACTIVATED_TAG = 3
integer, parameter :: GCC_TAG = 4

integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))
integer, parameter :: jumpvec2D(3,8) = reshape((/ 1,0,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, -1,-1,0, 0,-1,0, 1,-1,0 /), (/3,8/))
integer, parameter :: nfcell = 10, nfout = 11, nfvec = 12, nfpath = 13, nfres = 14, nfdcbind = 15, &
						nftraffic = 16, nfrun = 17, nftravel = 18, nfcmgui = 19, nfpos = 20, nflog=21, nfchemo=22
real, parameter :: BIG_TIME = 100000
real, parameter :: BALANCER_INTERVAL = 10
real, parameter :: DELTA_T = 0.25       ! minutes
logical, parameter :: use_gaplist = .true.
logical, parameter :: use_add_count = .true.    ! keep count of sites to add/remove, do the adjustment at regular intervals 
logical, parameter :: save_input = .true.
integer, parameter :: MAX_DC = 1000
integer, parameter :: MAX_FDC = 1000
integer, parameter :: DCDIM = 4         ! MUST be an even number
real, parameter :: DCRadius = 2		! (grids) This is just the approx size in the lattice, NOT the SOI
real, parameter :: FDCRadius = 2	! (grids) This is just the approx size in the lattice, NOT the SOI

! Diffusion parameters
logical, parameter :: use_ode_diffusion = .false.	! otherwise use the original method in fields.f90 
integer, parameter :: NDIFFSTEPS = 6    ! divisions of DELTA_T for diffusion computation

! B cell parameters
logical, parameter :: random_cognate = .false.          ! number of cognate seed cells is random or determined
integer, parameter :: MMAX_GEN = 20     ! max number of generations (for array dimension only)

! B cell region
integer, parameter :: FOLLICLE = 1
integer, parameter :: GONE = 2
real, parameter :: ELLIPSE_RATIO = 2.0
real, parameter :: ENTRY_ALPHA = 0.5
real, parameter :: EXIT_ALPHA = 0.5
!integer, parameter :: BASE_NFDC = 50
logical, parameter :: use_FDCs = .true.

! Differentiation probabilities
real, parameter :: PLASMA_PROB = 0.4

! Data above this line almost never change
!==============================================================================================================

! Run parameters

! Parameters and switches for calibration
logical, parameter :: calibrate_motility = .false.
logical, parameter :: motility_param_range = .false.
logical, parameter :: motility_save_paths = .false.
logical, parameter :: calibrate_diffusion = .false.
logical, parameter :: compute_travel_time = .false.
integer, parameter :: n_multiple_runs = 1

! Parameters and switches for testing
logical, parameter :: test_vascular = .false.
logical, parameter :: test_squeezer = .false.
logical, parameter :: turn_off_chemotaxis = .false.		! to test the chemotaxis model when cells are not attracted to exits
logical, parameter :: test_chemotaxis = .false.			! to test multiple chemotactic effects on one or a few cells

! Debugging parameters
integer :: idbug = 0

logical, parameter :: evaluate_residence_time = .false.
integer, parameter :: istep_res1 = 4*60*24*3			! 3 days (was 5000)
integer, parameter :: istep_res2 = istep_res1 + 4*60*24	! 1 day of tagging

! Parameters for controlling data capture for graphical purposes
logical, parameter :: save_pos_cmgui = .false.          ! To make movies
integer, parameter :: ntres = 60    ! 15 min
logical, parameter :: log_traffic = .true.

type vector3_type
	real :: x, y, z
end type

! Ex-globalvar definitions
integer :: NBcells0
integer :: NBcells
integer :: Nsites
integer :: NDC
integer :: NDCalive
integer :: NFDC
type(vector3_type) :: Radius
real :: InflowTotal
real :: OutflowTotal
real :: VEGFmass
real :: Vascularity
real :: dVdt
real :: Cvegf

!-------------------------------------------------------------
! Result type
type result_type
    integer :: dN_EffCogBC(NCTYPES)
    integer :: dN_EffCogBCGen(MMAX_GEN)
    integer :: N_EffCogBC(NCTYPES)
    integer :: N_EffCogBCGen(MMAX_GEN)
    integer :: dN_Dead
    integer :: N_Dead
end type

type cog_type
    sequence
	integer :: ID			! ID of the originating naive cell
	real :: dietime			! time that the cell dies
	real :: dividetime		! time that the cell divides
	real :: stagetime		! time that a cell can pass to next stage
	integer(2) :: generation
	integer(2) :: stage
	integer(2) :: region
	integer(2) :: status	! BCL6_LO, BCL6_HI, PLASMA
	integer :: cogID		! index in the list of cognate cells
end type

type cell_type
    sequence
    integer :: ID
    logical :: exists
    integer :: site(3)
    integer :: step
    integer(2) :: ctype
	integer(2) :: lastdir
    real :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
    real :: receptor_level(MAX_RECEPTOR)  ! level of receptor (susceptibility to the chemokine signal)
    type(cog_type),pointer :: cptr    ! because NULL is used by winsock (from ifwinty).  NULLIFY() instead.
end type

type DC_type
    sequence
    integer :: ID               ! unique ID number
    integer :: site(3)          ! DC location
    integer :: nsites
    integer :: nbound           ! current number of bound T cells
    integer :: ncogbound        ! current number of bound cognate T cells
    real :: density             ! current antigen density
    real :: dietime             ! time DC will die
    real :: stimulation         ! amount of stimulation provided by a DC
    real :: secretion
    logical :: capable          ! can DC deliver TCR stimulation?
    logical :: alive            ! is DC alive?
end type

type FDC_type
    sequence
    integer :: ID               ! unique ID number
    integer :: site(3)          ! FDC location
    integer :: nsites
!    real :: density             ! current antigen density
!    real :: dietime             ! time DC will die
!    real :: stimulation         ! amount of stimulation provided by a DC
	real :: secretion
    logical :: alive            ! is FDC alive?
end type

type boundary_type
    integer :: site(3)
    logical :: chemo_influx(MAX_CHEMO)
    logical :: entry_ok
    logical :: exit_ok
    type (boundary_type), pointer :: next
end type

type occupancy_type
    integer(2) :: DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
    integer :: indx(2)
    integer :: FDC_nbdry				! number of FDCs that the site is adjacent to
    type (boundary_type), pointer :: bdry
end type

type exit_type
    integer :: ID               ! unique ID number
    integer :: site(3)          ! location
end type

type dist_type
	integer :: class
	real :: p1, p2, p3
end type

type counter_type
    logical :: logscale
    integer :: period
    integer :: nbins
    real :: binmin, binstep
    integer, allocatable :: ndist(:)
    real :: total
end type

!---------------------------------------------------
! Parameters to read from cell parameters input file
!---------------------------------------------------
real :: BC_AVIDITY_MEDIAN = 1.0         ! median B cell avidity
real :: BC_AVIDITY_SHAPE = 1.2			! shape -> 1 gives normal dist with small variance
real :: BC_COGNATE_FRACTION = 0.0001	! fraction of T cells that are cognate initially
real :: BC_STIM_RATE_CONSTANT = 0.83	! rate const for BCR stimulation (-> molecules/min)
real :: BC_STIM_HALFLIFE = 24			! hours
integer :: BC_MAX_GEN = 10              ! maximum number of TC generations

type(dist_type) :: divide_dist1
type(dist_type) :: divide_dist2

real :: BC_FRACTION
real :: BC_RADIUS
real :: BLOB_RADIUS
real :: FLUID_FRACTION
integer :: BASE_NFDC

real(DP) :: GAMMA                           ! controls crowding
real(DP) :: BETA                            ! speed: 0 < beta < 1
real(DP) :: RHO                             ! persistence: 0 < rho < 1

logical :: use_traffic = .true.
logical :: computed_outflow

real :: RESIDENCE_TIME                  ! T cell residence time in hours -> inflow rate
! Vascularity parameters
real :: Inflammation_days1 = 4          ! Days of plateau level - parameters for VEGF_MODEL = 1
real :: Inflammation_days2 = 5          ! End of inflammation
real :: Inflammation_level = 1.0		! This is the level of inflammation (scaled later by NBcells0)
integer :: VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity
real :: chemo_K_exit                    ! level of chemotactic influence towards exits

real :: days                            ! number of days to simulate
integer :: seed(2)                      ! seed vector for the RNGs
integer :: NT_GUI_OUT					! interval between GUI outputs (timesteps)
integer :: SPECIES						! animal species source of T cells
character*(128) :: fixedfile

!---------------------------------------------------
! end of parameters to read from input file
!---------------------------------------------------

!---------------------------------------------------
! More input parameters
!---------------------------------------------------
integer :: NX = 100

integer :: NLN_RESPONSE 				! number of LNs in the response

! B cell parameters
real :: BC_life_median1					! median lifetime of naive T cells
real :: BC_life_median2					! median lifetime of activated T cells
real :: BC_life_shape					! shape parameter for lifetime of T cells
integer :: NBC_LN = 3.0e07				! number of B cells in a LN
integer :: NBC_BODY = 1.6e09			! number of circulating B cells in the whole body
integer :: BC_TO_DC = 1000				! ratio of B cells to DCs
integer :: BC_TO_FDC = 600				! ratio of B cells to FDCs

! T cell parameters
logical :: TCR_splitting = .false.      ! enable sharing of integrated TCR signal between progeny cells
real :: transient_stagetime				! stagetime for TRANSIENT (Stage 1)
real :: clusters_stagetime				! stagetime for CLUSTERS (Stage 2)
real :: transient_bindtime				! bindtime for TRANSIENT (Stage 1)
real :: clusters_bindtime				! bindtime for CLUSTERS (Stage 2)
real :: swarms_bindtime					! bindtime for SWARMS (Stage 3)
real :: TC_life_median1					! median lifetime of naive T cells
real :: TC_life_median2					! median lifetime of activated T cells
real :: TC_life_shape					! shape parameter for lifetime of T cells
integer :: NTC_LN = 3.0e07				! number of T cells in a LN
integer :: NTC_BODY = 1.6e09			! number of circulating T cells in the whole body
real :: K1_S1P1 = 0.005					! S1P1/CD69 system parameters
real :: K2_S1P1 = 0.05
real :: K1_CD69 = 0.04
real :: K2_CD69 = 0.01

! DC parameters
logical, parameter :: DC_motion = .false.
logical, parameter :: RANDOM_DCFLUX = .false.
integer, parameter :: NDCsites = 7		! Number of lattice sites occupied by the DC core (soma). In fact DC vol = 1400 = 5.6*250
real, parameter :: DC_DCprox = 2.0      ! closest placement of DCs, units DC_RADIUS (WAS 1.0 for ICB DCU paper)
real, parameter :: bdry_DCprox = 2.0	! closest placement of DC to bdry, units DC_RADIUS
real, parameter :: bdry_FDCprox = 4.0	! closest placement of DC to bdry, units DC_RADIUS

! Egress parameters
real :: exit_fraction = 1.0/1000.       ! number of exits as a fraction of T cell population
real :: Ksurfaceportal = 40			! calibration factor for number of surface portals

!---------------------------------------------------
! end of more parameters to be read from input file
!---------------------------------------------------

! Geometry data
integer :: NY, NZ
integer, allocatable :: zdomain(:),zoffset(:)
integer :: blobrange(3,2)
real :: DELTA_X, PI
real :: TagRadius
real :: x0,y0,z0   ! centre in global coordinates (units = grids)
real :: Centre(3)
real :: Vc, Ve

! Motility data
integer :: nreldir, njumpdirs
integer :: jumpvec(3,27)    ! 14 is no-jump case (0,0,0)
integer :: reldir(6,MAXRELDIR)
real(DP) :: dirprob(0:MAXRELDIR)
integer :: DCoffset(3,NDCsites)

integer :: nreldir2D, njumpdirs2D
integer :: reldir2D(8,8)
real(DP) :: dirprob2D(0:8)
logical :: diagonal_jumps

! Chemotaxis data
integer :: chemo_N
real, allocatable :: chemo_r(:,:,:)
real, allocatable :: chemo_p(:,:,:,:)

! Cell data
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cellist(:)
integer, allocatable :: cognate_list(:)
integer, allocatable :: gaplist(:)
type(DC_type), allocatable :: DClist(:)
type(FDC_type), allocatable :: FDClist(:)
integer :: lastID, MAX_COG, lastcogID, nlist, n2Dsites, ngaps, ntagged=0, ID_offset, ncogseed, nbdry
integer :: lastNBcells, k_nonrandom
integer :: max_nlist, max_ngaps
integer :: nadd_sites, ndivisions
real :: Kdeath_lo(MMAX_GEN), Kdeath_hi(MMAX_GEN)
real :: lastbalancetime
real :: scale_factor	! scaling from model to one (or more) whole LNs
real :: Fcognate		! fraction of T cells in circulation that are cognate
logical :: use_cognate

! Result data
type(result_type) :: localres, totalres
integer :: noutflow_tag, ninflow_tag
real :: restime_tot
real, allocatable :: Tres_dist(:)
logical :: firstSummary
integer :: logID

! Travel time computation
integer :: ntravel
integer :: N_TRAVEL_COG, N_TRAVEL_DIST
integer :: k_travel_cog
real, allocatable :: travel_cog(:), travel_dist(:,:,:)

! Miscellaneous data
logical :: initialized
integer ::  Nsteps, nsteps_per_min, istep
integer :: Mnodes
integer :: IDtest
integer :: total_in = 0, total_out = 0
integer :: ndivided(MMAX_GEN) = 0   ! to store times between divisions
real :: tdivided(MMAX_GEN) = 0
real :: BCRdecayrate                ! rate of decay of integrated BCR stimulation (/min)
logical :: vary_vascularity = .true.     ! to allow inflammation to change vascularity (false if VEGF_MODEL = 0)

! Vascularity parameters
real :: VEGF_alpha = 4.0e-7         ! rate constant for dependence on inflammation (/min) (alpha_G in hev.m) (was 5.0e-7)
real :: VEGF_beta = 5.0e-8			! rate constant for basal VEGF production (beta_G in hev.m) (was 4.0e-8)
real :: VEGF_decayrate = 0.002      ! VEGF decay rate (/min)	(was 0.002)
real :: vasc_maxrate = 0.001        ! max rate constant for vascularity growth (/min)  (was 0.003)
real :: vasc_beta = 2.0				! Hill function parameter
integer :: vasc_n = 2               ! Hill function exponent

real :: vasc_decayrate				! vascularity decay rate (/min) (deduced)
real :: VEGF_baserate				! base rate of production of VEGF (> 0 for VEGF-vascularity model VEGF_MODEL = 1)
real :: Cvegf0					    ! steady-state VEGF concentration (VEGF_MODEL = 1)

character*(128) :: inputfile
character*(128) :: outputfile
character*(2048) :: logmsg
TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send, simulation_start, par_zig_init
logical :: dbug = .false.
real :: base_exit_prob 					! testing different chemo_K

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps

contains

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))

end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real :: dummy
real(DP) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm(r)
real :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm2(r)
real :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real function cdistance(site)
integer :: site(3)
real :: r(3)

r = site - Centre
cdistance = norm(r)
end function

!--------------------------------------------------------------------------------
! Crude check for a site inside, just looking at the allowable ranges of x, y and z.
!--------------------------------------------------------------------------------
logical function inside_xyz(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY .or. site(3) < 1 .or. site(3) > NZ) then
    inside_xyz = .false.
else
    inside_xyz = .true.
endif
end function

!-----------------------------------------------------------------------------------------
! To determine if a site falls within the ellipsoid with major axis radius = Radius%x
! Currently the ellipsoid is shaped like a flattened football, i.e. the ellipsoid equation is:
! (x^2+ z^2)/a^2 + y^2/b^2 <= 1
! Axis orientation: 
! The y-axis is on the line passing through the centres of the paracortex and the follicle.
! The x-axis is parallel to one long axis of the follicle.
! The z-axis is perpendicular to the x- and y-axes. 
!
!-----------------------------------------------------------------------------------------
logical function InsideEllipsoid(site)
integer :: site(3)
real :: r(3)

r = site - Centre
if (r(1)*r(1)/(Radius%x*Radius%x) + r(2)*r(2)/(Radius%y*Radius%y) + r(3)*r(3)/(Radius%z*Radius%z) <= 1) then
    InsideEllipsoid = .true.
else
    InsideEllipsoid = .false.
endif
end function

!--------------------------------------------------------------------------------
! A site is taggable if it is less than a specified distance TagRadius from
! the centre.
!--------------------------------------------------------------------------------
logical function taggable(site)
integer :: site(3)
real :: d, r(3)

r = site - Centre
d = norm(r)
if (d <= TagRadius) then
    taggable = .true.
else
    taggable = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! Returns the user-defined structure type to be used with a cell of type ctype.
!-----------------------------------------------------------------------------------------
integer function struct_type(ctype)
integer :: ctype

select case(ctype)
case(NONCOG_TYPE_TAG,TAGGED_CELL,RES_TAGGED_CELL)
    struct_type = NONCOG_TYPE_TAG
case(COG_CD4_TAG,COG_CD8_TAG)
    struct_type = COG_TYPE_TAG
case default
    write(logmsg,*) 'ERROR: struct_type: unrecognised cell type: ',ctype
    call logger(logmsg)
    stop
end select
end function

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
! A gap is indicated by cellist(k)%ID = 0.  Note that a gap should never correspond
! to a cognate cell, since we keep info for cognate cells after they leave or die.
! If use_gaplist = false, ngaps is always = 0, and there is no squeezing.
!-----------------------------------------------------------------------------------------
subroutine squeezer(force)
logical :: force
integer :: last, k, site(3), indx(2), i, n, region, res
logical :: ok

!write(*,*) 'squeezer'
if (ngaps == 0) return
if (.not.force .and. (ngaps < max_ngaps/10)) return
write(nflog,*) 'squeezer: ',ngaps,max_ngaps,nlist

n = 0
do k = 1,nlist
    if (cellist(k)%ID == 0) then    ! a gap
		if (associated(cellist(k)%cptr)) then
			write(*,*) 'Error: squeezer: gap cell is cognate: ',k,cellist(k)%cptr%cogID
			stop
		endif
        n = n+1
    endif
enddo
if (n /= ngaps) then
	write(*,*) 'Error: squeezer: inconsistent gap counts: ',n,ngaps
	stop
endif
last = nlist
k = 0
n = 0
do
    k = k+1
    if (cellist(k)%ID == 0) then    ! a gap
        if (k == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: k: ',k
                stop
            endif
            if (cellist(last)%ID == 0) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        call copycell2cell(cellist(last),cellist(k),k)
		if (associated(cellist(last)%cptr)) then
			region = get_region(cellist(last)%cptr)
			deallocate(cellist(last)%cptr)
			if (.not.associated(cellist(k)%cptr)) then
				write(*,*) 'Error: squeezer: cell copy target cptr npt associated: ',last,k
				stop
			endif
		else
			region = FOLLICLE
		endif
		
		cellist(last)%ID = 0
		cellist(last)%exists = .false.
		
		if (region == FOLLICLE) then
	        site = cellist(last)%site
	        indx = occupancy(site(1),site(2),site(3))%indx
	        do i = 1,2
	            if (indx(i) == last) indx(i) = k
	        enddo
	        occupancy(site(1),site(2),site(3))%indx = indx
	    endif
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
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
    if (.not.associated(cellist(kcell)%cptr)) then
		write(*,*) 'check_cognate_list: cptr not associated: istep: ',istep,lastcogID,k,kcell
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Copy the contents of cellist(kfrom) to the entry cellist(kto)
! Need to fix up the corresponding cognate_list() entry as well.
! Note that:
!   cell = cellist(kcell)
!   kcog = cell%cptr%cogID
!   k = cognate_list(kcog)
! =>
!   k = kcell
! The new version is safer since it protects against failing to change this subroutine
! when the components of cell_type are changed (as long as there is only one pointer, %cptr)
!-----------------------------------------------------------------------------------------
subroutine copycell2cell(cell_from,cell_to,kcell)
integer :: kcell
type(cell_type) :: cell_from, cell_to
integer :: ctype, stype, kcog
integer :: region_to, region_from
logical :: new_version = .true.
logical :: check = .false.

ctype = cell_from%ctype
stype = struct_type(ctype)

if (new_version) then
	cell_to = cell_from
	nullify(cell_to%cptr)
endif
if (stype == NONCOG_TYPE_TAG .and. associated(cell_to%cptr)) then
    deallocate(cell_to%cptr)
endif
if (stype == COG_TYPE_TAG) then
    if (.not.associated(cell_to%cptr)) then
        allocate(cell_to%cptr)
    endif
    cell_to%cptr = cell_from%cptr
    kcog = cell_to%cptr%cogID
    cognate_list(kcog) = kcell
elseif (stype /= NONCOG_TYPE_TAG) then
    write(logmsg,*) 'ERROR: copycell2cell: istep, ID, ctype, stype: ',istep,cell_from%ID,ctype,stype
    call logger(logmsg)
    stop
endif
if (.not.new_version) then
	cell_to%ID = cell_from%ID
	cell_to%exists = cell_from%exists
	cell_to%site = cell_from%site
	cell_to%step = cell_from%step
	cell_to%ctype = cell_from%ctype
	cell_to%lastdir = cell_from%lastdir
	cell_to%entrytime = cell_from%entrytime
	cell_to%receptor_level = cell_from%receptor_level
endif
if (check) then
	if (cell_to%ID /= cell_from%ID) then
	    write(logmsg,*) 'copycell2cell: check: bad ID: ',cell_to%ID,cell_from%ID
	    call logger(logmsg)
	    stop
	endif
	if (cell_to%exists /= cell_from%exists) then
	    write(logmsg,*) 'copycell2cell: check: bad exists: ',cell_to%exists,cell_from%exists
	    call logger(logmsg)
	    stop
    endif
    if (stype == COG_TYPE_TAG) then
        region_to = get_region(cell_to%cptr)
        region_from = get_region(cell_from%cptr)
        if (region_to /= region_from) then
    	    write(logmsg,*) 'copycell2cell: check: bad region: ',region_to,region_from
	        call logger(logmsg)
	        stop
	    endif
	endif       
endif
if (cell_from%ctype == 0) then
    write(logmsg,*) 'ERROR: copycell2cell: ctype = 0'
    call logger(logmsg)
    stop
endif
end subroutine

!--------------------------------------------------------------------------------
! Returns:
! 0 if no slots are occupied
! 1 if slot 1 is occupied
! 2 if slot 2 is occupied
! 3 if both slots are occupied
!--------------------------------------------------------------------------------
integer function getslots(site)
integer :: site(3)
integer :: k

getslots = 0
do k = 1,2
    if (occupancy(site(1),site(2),site(3))%indx(k) > 0) then
        getslots = getslots + k
    endif
enddo
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
integer function blobNeighbours(x,y,z)
integer :: x,y,z
integer :: nb, k, xx, yy, zz

nb = 0
blobNeighbours = 0
if (x<1 .or. x>NX) return
if (y<1 .or. y>NY) return
if (z<1 .or. z>NZ) return
if (occupancy(x,y,z)%indx(1) < 0) return	! outside or DC
do k = 1,27
	if (k == 14) cycle
	xx = x + jumpvec(1,k)
	yy = y + jumpvec(2,k)
	zz = z + jumpvec(3,k)
	if (occupancy(xx,yy,zz)%indx(1) >= 0) nb = nb + 1
enddo
blobNeighbours = nb
end function

!-----------------------------------------------------------------------------------------
! Is site near a FDC?
! The criterion for a FDC site might be different from an exit site.
! prox = DC_DCprox*FDC_RADIUS for DC - DC
!-----------------------------------------------------------------------------------------
logical function tooNearFDC(site,kdc,prox)
integer :: site(3), kdc
real :: prox
integer :: idc
real :: r(3), d

if (kdc == 0) then
    tooNearFDC = .false.
    return
endif
do idc = 1,kdc
    if (.not.FDClist(idc)%alive) cycle
    r = site - FDClist(idc)%site
    d = norm(r)     ! units sites
    if (d < prox) then
        tooNearFDC = .true.
        return
    endif
enddo
tooNearFDC = .false.
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetRadius(n)
integer :: n
real :: aRadius

aRadius = (ELLIPSE_RATIO*n*3/(4*PI))**0.33333
Radius%x = aRadius
Radius%y = aRadius/ELLIPSE_RATIO
Radius%z = aRadius
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the centre to the ellipsoid surface (approx) in the direction given by r(:)
! v(:) = (x,y,z) is a point on the surface if x^2/a^2 + (y^2 + z^2)/b^2 = 1
! Let v(:) = beta*r(:), then
! r(1)^2/a^2 + (r(2)^2 + r(3)^2)/b^2 = 1/beta^2
! The "radius" in the specified direction is then beta*|r| = beta*norm(r)
!-----------------------------------------------------------------------------------------
real function EllipsoidRadius(r)
real :: r(3)
real :: beta

beta = 1./(r(1)*r(1)/(Radius%x*Radius%x) + r(2)*r(2)/(Radius%y*Radius%y) + r(3)*r(3)/(Radius%z*Radius%z))
EllipsoidRadius = beta*norm(r)
end function

!-----------------------------------------------------------------------------------------
! prox is the minimum distance from the ellipsoid boundary (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearBdry(site,prox)
integer :: site(3)
real :: prox
real :: d, r(3), eradius

r = site - Centre
eradius = EllipsoidRadius(r)
d = cdistance(site)
if (eradius - d < prox) then
    tooNearBdry = .true.
    return
endif
tooNearBdry = .false.
end function

!-----------------------------------------------------------------------------------------
! Locate a free slot in a site adjacent to site1: site2 (freeslot)
! Returns freeslot = 0 if there is no free space in an adjacent site.
!-----------------------------------------------------------------------------------------
subroutine get_free_slot(xlim,site1,site2,freeslot)
integer :: xlim, site1(3), site2(3), freeslot
logical :: Lwall, Rwall
integer :: i, indx2(2)

if (site1(1) == 1) then
    Lwall = .true.
else
    Lwall = .false.
endif
if (site1(1) == xlim) then
    Rwall = .true.
else
    Rwall = .false.
endif

if (site1(1) < 1) then
    write(logmsg,*) 'get_free_slot: bad site1: ',site1
	call logger(logmsg)
    stop
endif
do i = 1,27
    if (i == 14) cycle       ! i = 14 corresponds to the no-jump case
	site2 = site1 + jumpvec(:,i)
    if (Lwall .and. site2(1) < 1) then
        cycle
    endif
    if (Rwall .and. site2(1) > xlim) then
        cycle
    endif
    if (site2(2) < 1 .or. site2(2) > NY .or. site2(3) < 1 .or. site2(3) > NZ) cycle
	indx2 = occupancy(site2(1),site2(2),site2(3))%indx
	if (indx2(1) >= 0) then         ! not OUTSIDE_TAG or DC
	    if (indx2(1) == 0) then     ! slot 1 is free
	        freeslot = 1
            return
	    elseif (indx2(2) == 0) then ! slot 2 is free
	        freeslot = 2
            return
        endif
    endif
enddo
freeslot = 0
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
integer function neighbourhoodCount(site)
integer :: site(3)
integer :: site2(3), k, count

count = 0
do k = 1,27
	site2 = site + jumpvec(:,k)
	if (site2(1) < 1 .or. site2(1) > NX) cycle
	if (site2(2) < 1 .or. site2(2) > NY) cycle
	if (site2(3) < 1 .or. site2(3) > NZ) cycle
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) >= 0) then
		count = count + 1
	endif
enddo
neighbourhoodCount = count
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function cell_count()
integer :: kcell, ntot
ntot = 0
do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    ntot = ntot + 1
enddo
cell_count = ntot
end function

!-----------------------------------------------------------------------------------------
! The type (cognate or non-cognate) of a T cell can be random or strictly determined.
!-----------------------------------------------------------------------------------------
integer function select_cell_type(kpar)
integer :: kpar
integer :: nratio

if (random_cognate) then
    select_cell_type = random_cell_type(kpar)
else
	if (mod(istep,6*4*60) == 1) then	! update Fcognate every 6 hours - otherwise mod(k_nonrandom,nratio) messed up
		Fcognate = BC_COGNATE_FRACTION - (scale_factor*Ncogseed)/NBC_BODY
	endif
    nratio = 1./Fcognate
    k_nonrandom = k_nonrandom + 1
    if (mod(k_nonrandom,nratio) == 0) then
        select_cell_type = COG_TYPE_TAG
    else
        select_cell_type = NONCOG_TYPE_TAG
    endif
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function random_cell_type(kpar)
integer :: kpar
real(DP) :: R

R = par_uni(kpar)
if (R > Fcognate) then
    random_cell_type = NONCOG_TYPE_TAG
else
    random_cell_type = COG_TYPE_TAG
endif
end function

!-----------------------------------------------------------------------------------------
! Note: The number of B cells entering a follicle is a fraction of the number of
! T cells entering the LN, e.g. approx. 0.5/15, if there are 15 follicles.
! The influx rate is determined by Tres, as for T cells.
!-----------------------------------------------------------------------------------------
subroutine set_globalvar
real :: inflow0

if (use_traffic) then
    inflow0 = NBcells0*DELTA_T/(residence_time*60)
else
    inflow0 = 0
endif

call generate_traffic(inflow0)
end subroutine

!-----------------------------------------------------------------------------------------
! Total T cell inflow and outflow are generated from the vascularity and baseline
! inflow, inflow0.
!-----------------------------------------------------------------------------------------
subroutine generate_traffic(inflow0)
real :: inflow0
real :: act, expansion, actfactor, tnow
real :: inflow, outflow

! Note: if inflammation signal = 0 the vascularity (and inflow) should be constant
tnow = istep*DELTA_T
inflow = inflow0*Vascularity   ! level of vascularity (1 = steady-state)
outflow = NBcells*DELTA_T/(RESIDENCE_TIME*60)
InflowTotal = inflow
OutflowTotal = outflow
end subroutine

!-----------------------------------------------------------------------------------------
! Pre-programmed inflammation signal is flat at Inflammation_100k (for NBcells0 = 100k)
! for Inflammation_days1 then ends at Inflammation_days2
!-----------------------------------------------------------------------------------------
real function get_inflammation()
real :: tnow, plateau

tnow = istep*DELTA_T    ! mins
tnow = tnow/(24*60)     ! days
plateau = Inflammation_level*NBcells0
if (tnow < Inflammation_days1) then
    get_inflammation = plateau
elseif (tnow > Inflammation_days2) then
    get_inflammation = 0
else
    get_inflammation = plateau*(Inflammation_days2 - tnow)/(Inflammation_days2 - Inflammation_days1)
endif
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine set_stage(p,stage)
type(cog_type), pointer :: p
integer :: stage

p%stage = stage
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
integer function get_stage(p)
type(cog_type), pointer :: p

get_stage = p%stage
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine set_status(p,status)
type(cog_type), pointer :: p
integer :: status

p%status = status
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
integer function get_status(p)
type(cog_type), pointer :: p

get_status = p%status
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine set_region(p,region)
type(cog_type), pointer :: p
integer :: region

p%region = region
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
integer function get_region(p)
type(cog_type), pointer :: p

get_region = p%region
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine set_generation(p,gen)
type(cog_type), pointer :: p
integer :: gen

p%generation = gen
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
integer function get_generation(p)
type(cog_type), pointer :: p

get_generation = p%generation
end function

!-----------------------------------------------------------------------------------------
! Log the stage and state of all cells in the lineage of the naive cell with %ID = id
! A cognate cell should never be removed from the list (i.e. ID > 0 always):
!   when it dies, then %exists = .false. %stage -> DEAD
! When a cognate cell leaves the blob it stays in the list (ID > 0) but %exists = .false.
!-----------------------------------------------------------------------------------------
subroutine show_lineage(id)
integer :: id
integer :: k, kcell, gen, stage, status
type(cog_type), pointer :: p
logical, save :: first = .true.

if (first) then
	write(nfout,*) 'Progeny of cognate cell: ID: ',id
	first = .false.
endif
write(nfout,*) 'lineage: istep: ',istep
do k = 1,lastcogID
	kcell = cognate_list(k)
	if (kcell < 1) then
		write(logmsg,*) 'Error: show_lineage: kcell < 1: ',kcell,k,logID
		call logger(logmsg)
		stop
	endif
	p => cellist(kcell)%cptr
	if (.not.associated(p)) then
		write(logmsg,*) 'Error: show_lineage: p not associated: ',kcell,k
		call logger(logmsg)
		stop
	endif
	if (p%ID /= id) cycle
	stage = get_stage(p)
	status = get_status(p)
	gen = get_generation(p)	
	write(nfout,'(i6,i6,3i4)') k, kcell, gen, stage, status
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Display the properties of a cognate cell.  The calling program must ascertain that kcell
! is cognate.
!-----------------------------------------------------------------------------------------
subroutine show_cognate_cell(kcell)
integer :: kcell
type(cog_type), pointer :: p
type(cell_type) :: bcell
integer :: gen, stage, region

bcell = cellist(kcell)
if (.not.associated(bcell%cptr)) then
    write(logmsg,*) 'ERROR: show_cognate_cell: cptr not associated: ',kcell
    call logger(logmsg)
    stop
endif
p => bcell%cptr
write(logmsg,*) 'Cognate cell: ',p%cogID,kcell,cellist(kcell)%ID
call logger(logmsg)
write(logmsg,'(a,i10,a,3i4,a,i2)') '  ID: ',bcell%ID,' site: ',bcell%site,' ctype: ',bcell%ctype
call logger(logmsg)
gen = get_generation(p)
stage = get_stage(p)
write(logmsg,'(a,i8,a,i2,a,i2)') '   cogID: ',p%cogID,' gen: ',gen,' stage: ', stage
call logger(logmsg)
write(logmsg,'(a,4f10.2)') '   times: entry,die,div,stage: ',bcell%entrytime,p%dietime,p%dividetime,p%stagetime
call logger(logmsg)

end subroutine

!-----------------------------------------------------------------------------------------
! Called whenever balancer carries out add_sites or removeSites.
! cognate_list(k) = 0 when the kth cognate cell has gone (left or died).
! If Mnodes = 1 this is called only once, after placeCells.  After that the list
! is maintained directly when cognate cells arrive or leave.
! This is a list of all the cognate cells in the follicle.
!-----------------------------------------------------------------------------------------
subroutine make_cognate_list(ok)
logical :: ok
integer :: kcell, ctype, stype, cogID
type(cog_type), pointer :: p

ok = .true.
cognate_list(1:lastcogID) = 0
do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    ctype = cellist(kcell)%ctype
    stype = struct_type(ctype)
    if (stype == COG_TYPE_TAG) then
        p => cellist(kcell)%cptr
        cogID = p%cogID
        if (cogID == 0) then
            lastcogID = lastcogID + 1
            if (lastcogID > MAX_COG) then
                write(logmsg,'(a,i6)') 'Error: make_cognate_list: cognate_list dimension exceeded: ',MAX_COG
                call logger(logmsg)
                ok = .false.
                return
            endif
            cogID = lastcogID
            p%cogID = cogID
        endif
        cognate_list(cogID) = kcell
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_pos
character*(8) :: msg
integer :: error

if (save_pos_cmgui) then
	call save_exnode
elseif (use_TCP) then
	call save_cell_positions
	msg = 'VTK'
	clear_to_send = .false.
    call winsock_send(awp_1,msg,len_trim(msg),error)
endif
end subroutine

!--------------------------------------------------------------------------------
! We display only B cells that are still in the FOLLICLE
!--------------------------------------------------------------------------------
subroutine save_cell_positions
integer :: k, kcell, site(3), j
integer :: ibcstate, stype, ctype, stage, region
real :: bcell_diam = 0.9
integer :: gen, bnd(2)
logical :: ex
character*(12) :: fname = 'cell_pos.dat'
character*(9) :: removefile = 'TO_REMOVE'

if (simulation_start) then
	inquire(file=fname,exist=ex)
	if (ex) then
		call unlink(fname)
	endif
	inquire(file=removefile,exist=ex)
	if (ex) then
		call unlink(removefile)
	endif
endif
simulation_start = .false.

if (.not.clear_to_send) then
	! wait until the file called removefile exists, then remove it
	inquire(file=removefile,exist=ex)
	if (.not.ex) then
		do
			inquire(file=removefile,exist=ex)
			if (.not.ex) then
!				call millisleep(10) ! no good at all
			else
				exit
			endif
		enddo
	endif
	call unlink(removefile)
	clear_to_send = .true.
endif

! B cell section
do k = 1,lastcogID
	kcell = cognate_list(k)
	if (kcell > 0) then
		region = get_region(cellist(kcell)%cptr)
		if (region /= FOLLICLE) cycle
		site = cellist(kcell)%site
		gen = cellist(kcell)%cptr%generation
		ibcstate = gen
		write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',k-1, site, bcell_diam, ibcstate
	endif
enddo

write(nfpos,'(a2,i6)') 'E ',istep
close(nfpos)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_exnode
integer :: k, kcell, site(3), nd, j
real :: bcstate
real :: bcell_diam = 0.9
real :: spectrum_max = 10, spectrum_freefraction = 0.6
integer :: gen, stage, region, bnd(2)
character*(64) :: fname = '\CMGUI\DCU\dcu_00000.exnode'

write(fname(16:20),'(i5.5)') istep
write(*,*) fname
open(nfcmgui,file=fname,status='replace')

nd = 0

! B cell section
write(nfcmgui,*) 'Group name : bcell'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
write(nfcmgui,*) '    value.'
write(nfcmgui,*) '  3) bcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
		stage = get_stage(cellist(kcell)%cptr)
		region = get_region(cellist(kcell)%cptr)
		if (region /= FOLLICLE) cycle
        nd = nd+1
        site = cellist(kcell)%site
		gen = cellist(kcell)%cptr%generation
        bcstate = (gen-1.0)/(BC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
        write(nfcmgui,'(a,i8,3i4,f4.1,f6.2)') 'Node: ',nd, site, bcell_diam, bcstate
    endif
enddo

! Bond section
write(nfcmgui,*) 'Group name : Bond'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) vector, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  3) dcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'

close(nfcmgui)
end subroutine

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(DP) :: p(:)
integer :: k
real(DP) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!-----------------------------------------------------------------------------------------
! Given c (0 < c < 0.5), the fraction of the sphere volume that is cut off, compute
! the distance of the cut from the centre, d, as x = d/R, where R = sphere radius.
!-----------------------------------------------------------------------------------------
subroutine get_slice(c,x)
real :: c,x,xp
real, parameter :: epsilon = 0.0001
integer :: k

x = 0.5
xp = x
k = 0
do
    k = k+1
    if (k > 100) then
        write(logmsg,*) 'ERROR: get_slice: failed to converge'
		call logger(logmsg)
        stop
    endif
    x = 1 - 4*c/(2-(1+x)*x**3)
    if (abs(x-xp) < epsilon) exit
    xp = x
enddo

end subroutine

!--------------------------------------------------------------------------------
! Interpretation in terms of chemokine gradient:
! The weighting is on jump directions that are closest to the direction of the
! net chemokine gradient vector (cells move up the chemokine concentration gradient).
! The vector r provides the gradient vector direction.  In practice the net gradient
! vector vsum is first normalised, then scaled by chemo_N to yield (x,y,z) that
! can be used to access the chemo_p lookup table.  Note that only a small fraction
! of the table entries will ever be accessed - these are those (x,y,z) near the
! sphere boundary.
!--------------------------------------------------------------------------------
subroutine ChemoSetup
integer :: x, y, z, k, r(3), s(3)
real :: r2, s2, rmod,smod,cosa
real, allocatable :: w(:)

write(logmsg,*) 'ChemoSetup'
call logger(logmsg)
allocate(chemo_r(0:chemo_N,0:chemo_N,0:chemo_N))
allocate(chemo_p(-chemo_N:chemo_N,-chemo_N:chemo_N,-chemo_N:chemo_N,njumpdirs))
allocate(w(njumpdirs))

do x = -chemo_N,chemo_N
    do y = -chemo_N,chemo_N
        do z = -chemo_N,chemo_N
            r = (/x,y,z/)
            r2 = dot_product(r,r)
            rmod = sqrt(r2)
            w = 0
            do k = 1,njumpdirs
                if (k == 14) cycle
                s = jumpvec(:,k)
                s2 = dot_product(s,s)
                smod = sqrt(s2)
                cosa = dot_product(r,s)/(rmod*smod)
                if (cosa > 0) then
                    w(k) = cosa*cosa/smod
                endif
            enddo
            w = w/sum(w)
            chemo_p(x,y,z,:) = w	! Note that these probs sum to 1
            if (x >= 0 .and. y >= 0 .and. z >= 0) then
                chemo_r(x,y,z) = sqrt(r2)
            endif
        enddo
    enddo
enddo
deallocate(w)
!call test_chemo
!stop
end subroutine

!--------------------------------------------------------------------------------
! Computes the jump probabilities (absolute directions) accounting for chemotaxis
! On input:
!   p(:) holds the jump probabilities not accounting for  chemotaxis
!   v(:) was the site offset relative to the exit, used only to approximate direction
!        by allowing a discrete number of directions (chemo_N determines these, in fact
!        only a small subset of the array positions are used - those roughly falling
!        on a sphere of radius chemo_N)
!   f is the amount of chemotactic influence
! Note that f incorporates both the magnitude of chemokine gradients and the cell's
! susceptibility to chemotaxis.
! On return p(:) holds the modified jump probabilities.
! Note: when njumpdirs = 27, jump 14 corresponds to (0,0,0) - unused.
! Note: code modifications now base v and f on net chemotactic attraction of 
! multiple exits and DCs.
!--------------------------------------------------------------------------------
subroutine chemo_probs_pre(p,v,f)
real(DP) :: p(:)
integer :: v(:)
real :: f
integer :: k
real(DP) :: pc(MAXRELDIR+1)

if (f == 0) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
pc(1:njumpdirs) = chemo_p(v(1),v(2),v(3),:)
do k = 1,njumpdirs
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*p(k) + f*pc(k)
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! B cell chemotaxis
!--------------------------------------------------------------------------------
subroutine chemo_probs(p,v,f)
real(DP) :: p(:)
real :: v(3), f
integer :: k
real(DP) :: pc(MAXRELDIR+1)
real :: r(3), r2, rmod, s(3), s2, smod, cosa
logical :: dbug

if (f == 0) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs: f > 0: ',f
	call logger(logmsg)
	return
endif
dbug = .false.
p = p/sum(p)

r = v
r2 = dot_product(r,r)
rmod = sqrt(r2)
pc = 0
do k = 1,njumpdirs
    if (k == 14) cycle
    s = jumpvec(:,k)
    s2 = dot_product(s,s)
    smod = sqrt(s2)
    cosa = dot_product(r,s)/(rmod*smod)
    if (cosa > 0) then
        pc(k) = cosa*cosa/smod
    endif
enddo
pc = pc/sum(pc)
if (dbug) then
	write(*,*) 'pc:'
	write(*,'(10f7.3)') pc
endif
    
do k = 1,njumpdirs
    if (k == 14) cycle
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*p(k) + f*pc(k)
    endif
enddo
if (dbug) then
	write(*,*) 'p: f: ',f
	write(*,'(10f7.3)') p
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    endif
else
	write(*,*) trim(msg)
endif
inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL, INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL                   :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_xyz(k)
integer :: k
integer :: kcell, xyzsum(3)

xyzsum = 0
do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    xyzsum = xyzsum + cellist(kcell)%site
enddo
write(nfres,'(2i6,3i12)') istep,k,xyzsum
end subroutine


!--------------------------------------------------------------------------------
! Returns the number of cells on the site offset by (x,y,z) from site(:)
!--------------------------------------------------------------------------------
integer function sitecells(site,x,y,z)
integer :: site(3),x,y,z
integer :: k,indx(2)

indx = occupancy(site(1)+x,site(2)+y,site(3)+z)%indx
sitecells = 0
do k = 1,2
    if (indx(k) > 0) sitecells = sitecells + 1
enddo
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_rng
integer :: kpar = 0
integer :: i
real(DP):: R

do i = 1,10
    R = par_uni(kpar)
    write(*,*) i,R,par_shr3(kpar)
enddo

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_zdistribution
real :: tot(NZ)
integer :: kcell, site(3)

tot = 0
do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    site = cellist(kcell)%site
    tot(site(3)) = tot(site(3)) + 1
enddo
tot = tot/nlist
end subroutine

!-----------------------------------------------------------------------------------------
! For a location xyz, check that the occupancy() info is consistent with the info in
! the cellist() entries corresponding to occupancy()%indx (if any).  There could be
! 0, 1, or 2 cellist() entries.
!-----------------------------------------------------------------------------------------
subroutine checkslots(msg,xyz)
character*(*) :: msg
integer :: xyz(3)
integer :: indx(2), k, cells, slots, site(3)
logical :: occ(2)

occ = .false.
cells = 0
slots = getslots(xyz)
indx = occupancy(xyz(1),xyz(2),xyz(3))%indx
do k = 1,2
    if (indx(k) > 0) then
        cells = cells + k
        occ(k) = .true.
        site = cellist(indx(k))%site
        if (xyz(1) /= site(1) .or. xyz(2) /= site(2) .or. xyz(3) /= site(3)) then
            write(*,'(a,a,8i6)') msg,' checkslots: site error: ',k,xyz,site
            stop
        endif
    elseif (indx(1) < 0) then	! outside or DC
        write(*,*) msg,' checkslots: indx: ',xyz,indx
        stop
    endif
enddo
if (slots /= cells) then
    write(*,'(a,a,6i4,2L2)') msg,' checkslots: mismatch: ',xyz,slots,cells,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Checks global cell locations for tag validity.
!-----------------------------------------------------------------------------------------
subroutine check_tagged
integer :: k,d2,n,site(3)
type(cell_type) :: cell

do k = 1,nlist
    cell = cellist(k)
    if (cell%ctype == TAGGED_CELL) then
        n = n+1
        site = cell%site
        if (.not.taggable(site)) then
            write(*,*) 'Bad tagged cell: ',site,d2,Radius%x*Radius%x
        endif
    endif
enddo
write(*,*) 'did check_tagged'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check(msg,x,y,z)
character*(*) :: msg
integer :: x,y,z,slots,k,indx(2),site(3)
logical :: occ(2)

site = (/x,y,z/)
slots = getslots(site)
occ = .false.
indx = occupancy(x,y,z)%indx
do k = 1,2
    if (cellist(indx(k))%ctype > 0) occ(k) = .true.
enddo
if (slots == 1 .and. (.not.occ(1) .or. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 2 .and. (.not.occ(2) .or. occ(1))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 3 .and. .not.(occ(1) .and. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checker
integer :: kcell,x,y,z,slots,slot,cells,k,indx(2),site(3)
integer, allocatable :: tot(:)
logical :: occ(2)

allocate(tot(NX))
do x = 1,NX
    tot(x) = 0
    do y = 1,NY
        do z = 1,NZ
            site = (/x,y,z/)
            slots = getslots(site)
            cells = 0
            occ = .false.
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                if (indx(k) > 0) then
                    tot(x) = tot(x) + 1
                    cells = cells + k
                    occ(k) = .true.
                    site = cellist(indx(k))%site
                    if (x /= site(1) .or. y /= site(2) .or. z /= site(3)) then
                        write(*,'(a,2i2,2i7,6i4)') 'checker: site error: ',k,indx(k),cellist(indx(k))%ID,x,y,z,site
                        stop
                    endif
                endif
            enddo
            if (slots /= cells) then
                write(*,'(a,6i4,2L2)') 'checker: mismatch: ',x,y,z,slots,cells,occ
                stop
            endif
        enddo
    enddo
enddo

do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    site = cellist(kcell)%site
    indx = occupancy(site(1),site(2),site(3))%indx
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: checker: bad indx: ',kcell,site,indx
        stop
    endif
enddo
deallocate(tot)
write(*,*) 'checked OK: ',' nlist: ',nlist
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcellcount(ok)
logical :: ok
integer :: kcell, ntot

ok = .true.
ntot = 0
do kcell = 1,nlist
	if (cellist(kcell)%exists /= (cellist(kcell)%ID /= 0)) then
		if (.not.associated(cellist(kcell)%cptr)) then
			write(logmsg,*) 'Error: checkcellcount: inconsistent %exists, %ID: ',kcell,cellist(kcell)%exists,cellist(kcell)%ID
			call logger(logmsg)
			ok = .false.
			return
		endif
	endif
    if (.not.cellist(kcell)%exists) cycle
	ntot = ntot + 1
enddo
if (ntot /= NBcells) then
	write(logmsg,*) 'Error: inconsistent cell counts: NBcells, ntot: ',NBcells,ntot
	call logger(logmsg)
	ok = .false.
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcellsite(kcell)
integer :: kcell
integer :: id,site(3)

site = cellist(kcell)%site
id = cellist(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'cell: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
site = cellist(kcell)%site
id = cellist(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'big_: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_big_occ
integer :: x,y,z,k,kcell,indx(2)
logical :: OK

OK = .true.
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                kcell = indx(k)
                if (kcell < 0) then
                    if (kcell /= OUTSIDE_TAG) then
                        OK = .false.
                        write(*,*) x,y,z,k,kcell
                    endif
                endif
            enddo
        enddo
    enddo
enddo
if (.not.OK) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcell(str)
character*(*) :: str
integer :: ictest = 12

write(*,*) 'checkcell: ',str,'  ',ictest
!write(*,*) cellist(ictest)%ID,cellist(ictest)%PTR1%entrytime
end subroutine

!-----------------------------------------------------------------------------------------
! Checks the two sites (before and after jump) to see if the one free slot is always #2
! HAPPENS ALL THE TIME
!-----------------------------------------------------------------------------------------
subroutine check_site_indx(site1,site2)
integer :: site1(3),site2(3)
integer :: indx(2)

indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
indx = occupancy(site2(1),site2(2),site2(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_cells
type(cog_type), pointer :: p
integer :: kcell, ctype, n1, n2

n1 = 0
n2 = 0
do kcell = 1,nlist
    if (.not.cellist(kcell)%exists) cycle
    p => cellist(kcell)%cptr
    if (associated(p)) then
		n1 = n1 + 1
	endif
    ctype = cellist(kcell)%ctype
    if (ctype == COG_TYPE_TAG) then
		n2 = n2 + 1
	endif
enddo
write(*,*) 'n1, n2: ',n1,n2
end subroutine

end module
