!==========================================================================================
! Modification notes
! ------------------

!==========================================================================================

module global

use omp_lib
!use omp_CD69
!use omp_IL2
!use omp_IL7
!use omp_IL_dummy
use par_zig_mod
use winsock

implicit none

INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )

! General parameters
integer, parameter :: BIG_INT = 2**30
integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MODEL = MOORE26_MODEL
integer, parameter :: MAXRELDIR = 26
integer, parameter :: TRAFFIC_MODE_1 = 1
integer, parameter :: TRAFFIC_MODE_2 = 2
integer, parameter :: EXIT_EVERYWHERE = 1
integer, parameter :: EXIT_LOWER_SURFACE = 2
integer, parameter :: EXIT_BLOB_PORTALS = 3

!integer, parameter :: TASK_TAG = 1
!integer, parameter :: CELL_TAG = 2
!integer, parameter :: LOC_TAG  = 3
!integer, parameter :: OCC_TAG  = 4
!integer, parameter :: OCNT_TAG = 5
!integer, parameter :: CCNT_TAG = 6
!integer, parameter :: ADD_TAG   = 8
!integer, parameter :: CYT_L2R_TAG   = 10
!integer, parameter :: CYT_R2L_TAG   = 11
!integer, parameter :: GLOBAL1_TAG   = 12
!integer, parameter :: GLOBAL2_TAG   = 13
!integer, parameter :: WX_TAG   = 19
!integer, parameter :: CYT_INFO_TAG   = 20
!integer, parameter :: CYT_DATA_TAG   = 21

!integer, parameter :: RES1_TAG   = 22
!integer, parameter :: RES2_TAG   = 23
!integer, parameter :: RES3_TAG   = 24
!integer, parameter :: RES4_TAG   = 25
!integer, parameter :: RES5_TAG   = 26
!integer, parameter :: VIS1_TAG   = 27
!integer, parameter :: VIS2_TAG   = 28
!integer, parameter :: MOL1_TAG   = 29
!integer, parameter :: MOL2_TAG   = 30

! B cell activation stage
integer, parameter :: NAIVE		  = 1
integer, parameter :: ANTIGEN_MET = 2
integer, parameter :: CCR7_UP     = 3
integer, parameter :: TCELL_MET   = 4
integer, parameter :: DIVIDING    = 5
integer, parameter :: GCC_COMMIT  = 6
integer, parameter :: PLASMA      = 7
integer, parameter :: BCL6_UP     = 8
integer, parameter :: FINISHED    = 9
integer, parameter :: STAGELIMIT  = 9

integer, parameter :: BCL6_LO = 1
integer, parameter :: BCL6_HI = GCC_COMMIT

real, parameter :: T_CCR7_UP   = 2*60
real, parameter :: T_EBI2_UP = 2*60
real, parameter :: T_BCL6_UP   = 2*60
real, parameter :: T_FIRST_DIVISION = 12*60
real, parameter :: T_DIVISION       = 10*60

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging) 
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: TCP_PORT_2 = 5002
integer, parameter :: TCP_PORT_3 = 5003

integer, parameter :: STAGE_BYTE = 1
integer, parameter :: GENERATION_BYTE = 2
integer, parameter :: SLOT_NUM1 = 1
integer, parameter :: SLOT_NUM2 = 2
integer, parameter :: BOTH = 3

integer, parameter :: MAX_CYT = 6
!integer, parameter :: IL2_TAG  = 1
!integer, parameter :: IL4_TAG  = 2
!integer, parameter :: IL7_TAG  = 3
!integer, parameter :: IL9_TAG  = 4
!integer, parameter :: IL15_TAG = 5
!integer, parameter :: IL21_TAG = 6
!character*(5), parameter :: cyt_name(MAX_CYT) = (/ 'IL-2 ','IL-4 ','IL-7 ','IL-9 ','IL-15','IL-21' /)

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

!real, parameter :: receptor_level(5,4) = reshape((/ 1.,1.,1.,1.,0., .2,3.,1.,1.,0., .2,0.,3.,1.,0., 0.,0.,0.,1.,2. /), (/5,4/)) ! my guess
real, parameter :: receptor_level(5,4) = reshape((/ 1.,1.,1.,1.,0., .2,2.,2.,1.,0., .2,.5,2.,1.,0., 0.,0.,.2,1.,2. /), (/5,4/))	! Taka

integer, parameter :: NCTYPES = 4
integer, parameter :: NONCOG_TYPE_TAG  = 1
integer, parameter :: COG_TYPE_TAG  = 2
integer, parameter :: COG_CD4_TAG  = 2
integer, parameter :: COG_CD8_TAG  = 3
integer, parameter :: TAGGED_CELL = 100
integer, parameter :: RES_TAGGED_CELL = 101
!integer, parameter :: OUTSIDE_TAG = -BIG_INT + 1
integer, parameter :: OUTSIDE_TAG = -9999

integer, parameter :: NAIVE_TAG = 1
integer, parameter :: ANTIGEN_TAG = 2
integer, parameter :: ACTIVATED_TAG = 3
integer, parameter :: GCC_TAG = 3

integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))
integer, parameter :: jumpvec2D(3,8) = reshape((/ 1,0,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, -1,-1,0, 0,-1,0, 1,-1,0 /), (/3,8/))
integer, parameter :: nfcell = 10, nfout = 11, nfvec = 12, nfpath = 13, nfres = 14, nfdcbind = 15, &
						nftraffic = 16, nfrun = 17, nftravel = 18, nfcmgui = 19, nfpos = 20, nflog=21, nfchemo=22
real, parameter :: BIG_TIME = 100000
real, parameter :: BALANCER_INTERVAL = 10
integer, parameter :: SCANNER_INTERVAL = 100
real, parameter :: DELTA_T = 0.25       ! minutes
logical, parameter :: use_add_count = .true.    ! keep count of sites to add/remove, do the adjustment at regular intervals 
logical, parameter :: save_input = .true.
integer, parameter :: MAX_DC = 1000
integer, parameter :: MAX_FDC = 1000
integer, parameter :: DCDIM = 4         ! MUST be an even number
real, parameter :: DCRadius = 2		! (grids) This is just the approx size in the lattice, NOT the ROI
real, parameter :: FDCRadius = 2	! (grids) This is just the approx size in the lattice, NOT the ROI

! Diffusion parameters
logical, parameter :: use_cytokines = .false.
logical, parameter :: use_diffusion = .false.
integer, parameter :: NDIFFSTEPS = 6    ! divisions of DELTA_T for diffusion computation

! B cell parameters
integer, parameter :: traffic_mode = TRAFFIC_MODE_2
logical, parameter :: use_blob = .true.
logical, parameter :: random_cognate = .false.          ! number of cognate seed cells is random or determined
integer, parameter :: MMAX_GEN = 20     ! max number of generations (for array dimension only)
integer, parameter :: NGEN_EXIT = 8     ! minimum non-NAIVE T cell generation permitted to exit (exit_rule = 1)
real, parameter :: CHEMO_MIN = 0.05		! minimum level of chemotactic influence (at r = chemo_radius)
integer :: exit_rule = 1                ! 1 = use NGEN_EXIT, 2 = use EXIT_THRESHOLD, 3 = use S1P1
logical :: COMPUTE_OUTFLOW = .false.
!real, parameter :: exit_prox = 4

! B cell region
integer, parameter :: FOLLICLE = 1
real, parameter :: ELLIPSE_RATIO = 2.0
real, parameter :: ENTRY_ALPHA = 0.5
real, parameter :: EXIT_ALPHA = 0.5

! Differentiation probabilities
real, parameter :: PLASMA_PROB = 0.4

! GUI parameters
character*(12), parameter :: stopfile = 'stop_dll'
character*(13), parameter :: pausefile = 'pause_dll'

! Data above this line almost never change
!==============================================================================================================

! Run parameters

!logical, parameter :: vary_vascularity = .true.
!logical, parameter :: use_chemotaxis = .false. ! now based on exit_region == EXIT_CHEMOTAXIS
!logical, parameter :: fix_avidity = .true.
!integer, parameter :: avidity_nlevels = 8
!!logical, parameter :: avidity_logscale = .true.
!!real, parameter :: avidity_min = -1    ! if avidity_logscale then avidity_min is log10(actual min)
!!real, parameter :: avidity_step = 0.185861   ! and log10(actual max) = log10(actual min) + (nlevels-1)*avidity_step
!logical, parameter :: avidity_logscale = .false.
!real, parameter :: avidity_min = 0.3
!real, parameter :: avidity_step = 0.2

! Parameters and switches for calibration
logical, parameter :: calibrate_motility = .false.
logical, parameter :: motility_param_range = .false.
logical, parameter :: motility_save_paths = .false.
logical, parameter :: calibrate_diffusion = .false.
logical, parameter :: compute_travel_time = .false.
integer, parameter :: n_multiple_runs = 1

! Parameters and switches for testing
logical, parameter :: test_vascular = .false.
logical, parameter :: turn_off_chemotaxis = .false.		! to test the chemotaxis model when cells are not attracted to exits

! Debugging parameters
!logical, parameter :: dbug = .false.
logical, parameter :: dbug_cog = .false.
logical, parameter :: avid_debug = .false.
integer, parameter :: CHECKING = 0
integer :: idbug = 0

real, parameter :: TAGGED_CHEMO_FRACTION = 0.02		! was 0.1 for exit chemotaxis
real, parameter :: TAGGED_CHEMO_ACTIVITY = 1.0

! To investigate the effect of chemotaxis on residence time.
! The situation to be simulated is one in which most cells are not subject to exit chemotaxis,
! but a small fraction of tagged cells are.  The question is: what is the effect on the residence time
! of the tagged cells?
logical, parameter :: TAGGED_EXIT_CHEMOTAXIS = .false.	! ==> evaluate_residence_time = .true.
logical, parameter :: evaluate_residence_time = .false.
integer, parameter :: istep_res1 = 4*60*24*3			! 3 days (was 5000)
integer, parameter :: istep_res2 = istep_res1 + 4*60*24	! 1 day of tagging

! Parameters for controlling data capture for graphical purposes
logical, parameter :: save_pos_cmgui = .false.          ! To make movies
integer, parameter :: save_interval_hours = 48
real, parameter :: save_length_hours = 0.5      ! 30 minutes
logical, parameter :: generate_exnode = .false.
logical, parameter :: evaluate_stim_dist = .false.
integer, parameter :: ntaglimit = 100000
integer, parameter :: ntres = 60    ! 15 min
logical, parameter :: log_results = .false.
logical, parameter :: log_traffic = .true.
integer, parameter :: TCR_nlevels = 10
real, parameter :: TCR_limit = 2000
integer, parameter :: MAX_AVID_LEVELS = 30

!-------------------------------------------------------------
! Cytokine section
!integer, parameter :: N_CYT = 2
!integer, parameter :: CYT_TAG(1:N_CYT) = (/IL2_TAG, IL7_TAG/)
!integer, parameter :: CYT_NP = IL2_NP + IL7_NP
! diffusion calibration
!integer, parameter :: N_CYT = 1
!integer, parameter :: CYT_TAG(1:N_CYT) = (/IL2_TAG/)
!integer, parameter :: CYT_NP = 0    !IL2_NP
!logical, parameter :: CD25_SWITCH = .true.

! Ex-globalvar definitions
integer :: NBcells0
integer :: NBcells
integer :: Nsites
integer :: Nexits
integer :: Lastexit
integer :: NDC
integer :: NDCalive
integer :: NFDC
integer :: NFDCalive
real :: aRadius
real :: bRadius
real :: InflowTotal
real :: OutflowTotal
real :: VEGFmass
real :: Vascularity
real :: dVdt
real :: Cvegf

!-------------------------------------------------------------
! Result type
type result_type
    integer :: dN_EffCogTC(NCTYPES)
    integer :: dN_EffCogTCGen(MMAX_GEN)
    integer :: N_EffCogTC(NCTYPES)
    integer :: N_EffCogTCGen(MMAX_GEN)
    integer :: dN_Dead
    integer :: N_Dead
end type

type cog_type
    sequence
	real :: avidity			! level of TCR avidity with DC
	real :: stimulation		! TCR stimulation level
!    real :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
	real :: dietime			! time that the cell dies
	real :: dividetime		! time that the cell divides
	real :: stagetime		! time that a cell can pass to next stage
!	real :: stimrate        ! rate of TCR stimulation
!	real :: CD69            ! level of CD69 expression
!	real :: S1P1            ! level of S1P1 expression
!	real :: CCR7            ! level of CCR7 expression
!    real :: IL_state(CYT_NP)    ! receptor model state variable values
!    real :: IL_statep(CYT_NP)   ! receptor model state variable time derivative values
	integer(2) :: generation
	integer(2) :: stage
	integer(2) :: region
	integer(2) :: status	! BCL6_LO, BCL6_HI, PLASMA
!    integer :: status       ! holds data in bytes: 1=stage, 2=generation
	integer :: cogID		! index in the list of cognate cells
end type

type cell_type
    sequence
    integer :: ID
    integer :: site(3)
    integer :: step
    integer(2) :: ctype
	integer(2) :: lastdir
    real :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
    real :: receptor_level(MAX_RECEPTOR)  ! level of receptor (susceptibility to the chemokine signal)
!    type(cog_type),    pointer :: cptr => NULL()    ! pointer to cognate cell data
    type(cog_type),    pointer :: cptr    ! because NULL is used by winsock (from ifwinty).  NULLIFY() instead.
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
    logical :: capable          ! can DC deliver TCR stimulation?
    logical :: alive            ! is DC alive?
end type

type FDC_type
    sequence
    integer :: ID               ! unique ID number
    integer :: site(3)          ! FDC location
    integer :: nsites
!    integer :: nbound           ! current number of bound T cells
!    integer :: ncogbound        ! current number of bound cognate T cells
!    real :: density             ! current antigen density
!    real :: dietime             ! time DC will die
!    real :: stimulation         ! amount of stimulation provided by a DC
!    logical :: capable          ! can DC deliver TCR stimulation?
    logical :: alive            ! is FDC alive?
end type

type boundary_type
    integer :: site(3)
    logical :: chemo_influx(MAX_CHEMO)
    real :: chemo_rate(MAX_CHEMO)
    logical :: S1P
    real :: S1Prate
    logical :: CXCL13
    real :: CXCL13rate
    logical :: CCL21
    real :: CCL21rate
    logical :: OXY
    real :: OXYrate
    logical :: entry_ok
    logical :: exit_ok
!    type (boundary_type), pointer :: previous
    type (boundary_type), pointer :: next
end type

type occupancy_type
    integer(2) :: DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
!    integer(2) :: cDC(0:DCDIM-1)   ! cDC(0) = number of DCs within chemo_radius, cDC(k) = ID of kth DC
    integer :: indx(2)
    integer :: exitnum          ! > 0 => index of closest exit, = 0 => no close exit, < 0 => an exit
!    logical :: bdry
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

!TYPE winsockport
!    LOGICAL           :: is_open
!    INTEGER           :: handle
!    TYPE(T_IN_ADDR)   :: ip_addr    ! server IP address u_long: wp%ip_addr.S_addr = inet_addr(ip_address)
!    INTEGER   :: ip_addr    ! server IP address u_long: wp%ip_addr.S_addr = inet_addr(ip_address)
!    INTEGER           :: ip_port    ! IP port on server
!	INTEGER           :: protocol   ! TCP or UDP for ethernet
!END TYPE winsockport


!---------------------------------------------------
! Parameters to read from cell parameters input file
!---------------------------------------------------
real :: BC_AVIDITY_MEDIAN = 1.0         ! median B cell avidity
real :: BC_AVIDITY_SHAPE = 1.2			! shape -> 1 gives normal dist with small variance
real :: BC_COGNATE_FRACTION = 0.0001	! fraction of T cells that are cognate initially
real :: BC_STIM_RATE_CONSTANT = 0.83	! rate const for BCR stimulation (-> molecules/min)
real :: BC_STIM_HALFLIFE = 24			! hours
integer :: BC_MAX_GEN = 10              ! maximum number of TC generations

real :: TC_AVIDITY_MEDIAN = 1.0			! median T cell avidity
real :: TC_AVIDITY_SHAPE = 1.2			! shape -> 1 gives normal dist with small variance
real :: TC_CD8_FRACTION = 0.33			! fraction of all T cells that are CD8
real :: TC_COGNATE_FRACTION = 0.003	! fraction of T cells that are cognate initially
real :: TC_CONVERSION_TIME = 48			! time (in hours) over which T cells become cognate
real :: TC_STIM_RATE_CONSTANT = 0.83	! rate const for TCR stimulation (-> molecules/min)
real :: TC_STIM_WEIGHT = 0.1			! contribution of stimulation to act level
real :: TC_STIM_HALFLIFE = 24			! hours
integer :: TC_MAX_GEN = 20              ! maximum number of TC generations

real :: DC_ANTIGEN_MEAN = 10			! mean DC antigen density
real :: DC_ANTIGEN_SHAPE = 1.2			! DC antigen density shape param
real :: DC_ANTIGEN_MEDIAN				! median DC antigen density
!real :: DC_LIFETIME_MEAN = 3.5			! days
!real :: DC_LIFETIME_SHAPE  = 1.2		! days
!real :: DC_ACTIV_TAPER = 12				! time (hours) over which DC activity decays to zero
!real :: DC_BIND_DELAY = 2.5				! delay after unbinding before next binding (mins)
!real :: DC_BIND_ALFA = 0.95				! binding prob parameter
!real :: DC_MULTIBIND_PROB = 0.0			! reducing factor to bind prob for each current DC binding
!real :: DC_DENS_BY_STIM = 0.0002        ! rate of reduction of density by TCR stimulation
!real :: DC_DENS_HALFLIFE                ! half-life of DC activity (hours)

integer :: optionA
integer :: optionB
integer :: optionC
real :: IL2_PRODUCTION_TIME			    ! duration of IL2/CD25 production (hrs)
real :: IL2_THRESHOLD			        ! TCR stimulation needed to initiate IL-2/CD25 production
real :: ACTIVATION_THRESHOLD    		! combined stimulation needed for activation
real :: FIRST_DIVISION_THRESHOLD(2)		! activation level needed for first division
real :: DIVISION_THRESHOLD(2)			! activation level needed for subsequent division
real :: EXIT_THRESHOLD(2)               ! Activation/CD69 level S below which exit is permitted
real :: STIMULATION_LIMIT				! maximum activation level
real :: CD25_DIVISION_THRESHOLD         ! CD25 store level needed for division of activated cell
real :: CD25_SURVIVAL_THRESHOLD         ! CD25 store level needed for survival of activated cell

type(dist_type) :: divide_dist1
type(dist_type) :: divide_dist2
real :: CD8_DIVTIME_FACTOR = 1.5		! CD8 divide time as multiple of CD4 divide time

real :: BC_FRACTION
real :: BC_RADIUS
real :: TC_FRACTION
real :: TC_RADIUS
real :: BLOB_RADIUS
real :: FLUID_FRACTION

real(DP) :: GAMMA                           ! controls crowding
real(DP) :: BETA                            ! speed: 0 < beta < 1
real(DP) :: RHO                             ! persistence: 0 < rho < 1

logical :: use_traffic = .true.
logical :: use_chemotaxis
logical :: computed_outflow

real :: RESIDENCE_TIME                  ! T cell residence time in hours -> inflow rate
! Vascularity parameters
real :: Inflammation_days1 = 4          ! Days of plateau level - parameters for VEGF_MODEL = 1
real :: Inflammation_days2 = 5          ! End of inflammation
real :: Inflammation_level = 1.0		! This is the level of inflammation (scaled later by NBcells0)
integer :: exit_region                   ! determines blob region for cell exits
real :: efactor                         ! If constant_efactor = true, this is the factor for the p correction
integer :: VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity
real :: chemo_radius					! radius of chemotactic influence (um)
real :: chemo_K_exit                    ! level of chemotactic influence towards exits

logical :: fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
logical :: avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
integer :: avidity_nlevels              ! If fix_avidity, number of discrete avidity values
real :: avidity_min                     ! minimum value
real :: avidity_step                    ! step between equi-spaced values
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
!integer, parameter :: cDCDIM = 4        ! MUST be an even number
integer, parameter :: NDCsites = 7		! Number of lattice sites occupied by the DC core (soma). In fact DC vol = 1400 = 5.6*250
!integer, parameter :: NDCcore = 7      ! In fact DC vol = 1400 = 5.6*250
real, parameter :: DC_DCprox = 2.0      ! closest placement of DCs, units DC_RADIUS (WAS 1.0 for ICB DCU paper)
real, parameter :: bdry_DCprox = 2.0	! closest placement of DC to bdry, units DC_RADIUS
real, parameter :: bdry_FDCprox = 4.0	! closest placement of DC to bdry, units DC_RADIUS
!real, parameter :: exit_DCprox = 4.0    ! closest placement of DC to exit, units sites
!real, parameter :: exit_prox = 4.0      ! closest placement of exit to exit, units chemo_radius

! Egress parameters
real :: exit_fraction = 1.0/1000.       ! number of exits as a fraction of T cell population
real :: Ksurfaceportal = 40			! calibration factor for number of surface portals

!---------------------------------------------------
! end of more parameters to be read from input file
!---------------------------------------------------

! Cytokine data
!real, allocatable :: cyt(:,:,:,:)
!real, allocatable :: cyt_constit(:), cyt_mols(:), dcyt_mols(:)
!real :: cyt0(MAX_CYT) = (/3.0,1.0,1.0,1.0,1.0,1.0/)
!real :: K_diff(MAX_CYT), delta_diff(MAX_CYT)     ! = D.dt/(dx*dx)
!integer :: Ncytokines, cytokines(MAX_CYT), cyt_seq(MAX_CYT), NP_offset(MAX_CYT+1)
!real :: cyt_init(MAX_CYT), cyt_mean(MAX_CYT)

! Geometry data
integer :: NY, NZ
integer, allocatable :: xdomain(:),xoffset(:),zdomain(:),zoffset(:)
integer, allocatable :: nz_sites(:), nz_totsites(:), nz_cells(:), nz_excess(:)
real :: DELTA_X, PI
real :: TagRadius
real :: x0,y0,z0   ! centre in global coordinates (units = grids)
real :: Centre(3)
real :: Vc, Ve

! Motility data
integer :: nreldir, njumpdirs
integer :: jumpvec(3,27)    ! 14 is no-jump case (0,0,0)
!integer :: jumpvec2D(3,8)
integer :: reldir(6,MAXRELDIR)
real(DP) :: dirprob(0:MAXRELDIR)
integer :: DCoffset(3,NDCsites)

integer :: nreldir2D, njumpdirs2D
integer :: reldir2D(8,8)
real(DP) :: dirprob2D(0:8)
logical :: diagonal_jumps
real :: ep_factor      ! (25k) 2.4 when K1_S1P1 = 0.01, 2.8 when K1_S1P1 = 0.001  ! based on no-DC case!
                       ! (50k) 2.3 when K1_S1P1 = 0.01, chemo_K_exit = 0.3

! Chemotaxis data
integer :: chemo_N
real :: chemo_exp
real, allocatable :: chemo_r(:,:,:)
real, allocatable :: chemo_p(:,:,:,:)

! Cell data
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cellist(:)
integer, allocatable :: cognate_list(:)
integer, allocatable :: gaplist(:)
type(DC_type), allocatable :: DClist(:)
type(FDC_type), allocatable :: FDClist(:)
integer :: lastID, MAX_COG, lastcogID, nlist, n2Dsites, ngaps, ntagged=0, ID_offset, ncogseed, nbdry, NXcells
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
type(counter_type) :: avid_count,avid_count_total
logical :: firstSummary

! Travel time computation
integer :: ntravel
integer :: N_TRAVEL_COG, N_TRAVEL_DIST
integer :: k_travel_cog
real, allocatable :: travel_cog(:), travel_dist(:,:,:)

! Miscellaneous data
type(exit_type), allocatable :: exitlist(:)
integer :: max_exits        ! size of the exitlist(:) array
real :: min_transit_time = 60		! minimum time a T cell spends in the DCU (min)
real :: CD69_threshold				! level of CD69 below which egress can occur
real :: last_portal_update_time		! time that the number of exit portals was last updated
logical :: initialized, steadystate
integer :: navid = 0
integer ::  Nsteps, nsteps_per_min, istep
integer :: Mnodes
integer :: IDtest
integer :: total_in = 0, total_out = 0
integer :: nIL2thresh = 0           ! to store times to IL2 threshold
real :: tIL2thresh = 0
integer :: ndivided(MMAX_GEN) = 0   ! to store times between divisions
real :: tdivided(MMAX_GEN) = 0
real :: TCRdecayrate                ! rate of decay of integrated TCR stimulation (/min)
real :: BCRdecayrate                ! rate of decay of integrated BCR stimulation (/min)
real :: max_TCR = 0
real :: avidity_level(MAX_AVID_LEVELS)  ! discrete avidity levels for use with fix_avidity
logical :: vary_vascularity = .true.     ! to allow inflammation to change vascularity (false if VEGF_MODEL = 0)

integer :: check_egress(1000)		! for checking traffic
integer :: check_inflow

! Vascularity parameters
real :: VEGF_alpha = 4.0e-7         ! rate constant for dependence on inflammation (/min) (alpha_G in hev.m) (was 5.0e-7)
real :: VEGF_beta = 5.0e-8			! rate constant for basal VEGF production (beta_G in hev.m) (was 4.0e-8)
real :: VEGF_decayrate = 0.002      ! VEGF decay rate (/min)	(was 0.002)
!real :: vasc_maxrate = 0.0006       ! max rate constant for vascularity growth (/min)
real :: vasc_maxrate = 0.001       ! max rate constant for vascularity growth (/min)  (was 0.003)
!real :: vasc_beta = 1.5				! Hill function parameter
real :: vasc_beta = 2.0				! Hill function parameter
integer :: vasc_n = 2               ! Hill function exponent
! NOTE: all the parameter changes I've been trying do not get the peak N/N0 above 4, with 
! inflam = 1, days 3,4.

real :: vasc_decayrate				! vascularity decay rate (/min) (deduced)
real :: VEGF_baserate				! base rate of production of VEGF (> 0 for VEGF-vascularity model VEGF_MODEL = 1)
real :: Cvegf0					! steady-state VEGF concentration (VEGF_MODEL = 1)
real :: Kflow1                          ! Inflow dependence on expansion
real :: Kflow2                          ! Outflow dependence on expansion
real :: Kflow3                          ! Maximum amplification of inflow from DC activity
real :: Bflow                           ! Hill function parameter for amplification
integer :: Nflow                        ! Hill function power coefficient

character*(128) :: inputfile
character*(128) :: outputfile
!character*(128) :: resultfile
character*(2048) :: logmsg
TYPE(winsockport) :: awp_0, awp_1, awp_2, awp_3
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send, simulation_start, par_zig_init
logical :: dbug = .false.
logical :: use_portal_egress			! use fixed exit portals rather than random exit points
logical :: FIXED_NEXITS = .false.		! the number of exit portals is held fixed
!real :: base_exit_prob = 0.00097		! prob of exit of cell at bdry site (Tres = 24)
!real :: base_exit_prob = 0.00195		! prob of exit of cell at bdry site (Tres = 12) OK for NO chemotaxis
!real :: base_exit_prob = 0.0017		! prob of exit of cell at bdry site (Tres = 12) with chemotaxis (all suscept = 0.1)
!real :: base_exit_prob = 0.0014		! prob of exit of cell at bdry site (Tres = 12) with chemotaxis (all suscept > 0.1)
real :: base_exit_prob = 0.0014		! testing different chemo_K
real :: INLET_R_FRACTION = 0.7			! fraction of blob radius within which ingress occurs
real :: INLET_EXIT_LIMIT = 5			! if RELAX_INLET_EXIT_PROXIMITY, this determines how close an inlet point can be to an exit portal.

!DEC$ ATTRIBUTES DLLEXPORT :: ntravel, N_TRAVEL_COG, N_TRAVEL_DC, N_TRAVEL_DIST, k_travel_cog, k_travel_dc
!DEC$ ATTRIBUTES DLLEXPORT :: travel_dc, travel_cog, travel_dist
!DEC$ ATTRIBUTES DLLEXPORT :: nsteps	!istep
contains

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
!use IFPORT
integer :: n1,n2,kpar
integer :: k,R

!k = irand()     ! intrinsic
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
!    call random_number(R)
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
! To determine if a site falls within the ellipsoid with major axis radius = aRadius
! The ellipsoid equation is:
! x^2/a^2 + (y^2 + z^2)/b^2 <= 1
!-----------------------------------------------------------------------------------------
logical function insideEllipsoid(site)
integer :: site(3)
real :: r(3)

r = site - Centre
if (r(1)*r(1)/(aRadius*aRadius) + (r(2)*r(2) + r(3)*r(3))/(bRadius*bRadius) <= 1) then
    insideEllipsoid = .true.
else
    insideEllipsoid = .false.
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
!-----------------------------------------------------------------------------------------
subroutine squeezer(force)
logical :: force
integer :: last, k, site(3), indx(2), i, n, region

!write(*,*) 'squeezer'
if (ngaps == 0) return
if (.not.force .and. (ngaps < max_ngaps/2)) return
if (dbug) write(nflog,*) 'squeezer: ',ngaps,max_ngaps,nlist

n = 0
do k = 1,nlist
    if (cellist(k)%ID == 0) then    ! a gap
        n = n+1
!        write(*,*) 'gap at : ',k
    endif
enddo

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
!        cellist(k) = cellist(last)
		if (associated(cellist(last)%cptr)) then
!			call get_region(cellist(last)%cptr,region)
			region = get_region(cellist(last)%cptr)
		else
			region = FOLLICLE
		endif
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
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

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
!-----------------------------------------------------------------------------------------
subroutine copycell2cell(cell_from,cell_to,kcell)
integer :: kcell
type(cell_type) :: cell_from, cell_to
integer :: ctype, stype, kcog

ctype = cell_from%ctype
stype = struct_type(ctype)

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
cell_to%ID = cell_from%ID
cell_to%site = cell_from%site
cell_to%ctype = cell_from%ctype
cell_to%lastdir = cell_from%lastdir
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
! Distance from the centre to the ellipsoid surface (approx) in the direction given by r(:)
! v(:) = (x,y,z) is a point on the surface if x^2/a^2 + (y^2 + z^2)/b^2 = 1
! Let v(:) = beta*r(:), then
! r(1)^2/a^2 + (r(2)^2 + r(3)^2)/b^2 = 1/beta^2
! The "radius" in the specified direction is then beta*|r| = beta*norm(r)
!-----------------------------------------------------------------------------------------
real function EllipsoidRadius(r)
real :: r(3)
real :: beta

beta = 1./sqrt(r(1)*r(1)/(aRadius*aRadius) + (r(2)*r(2) + r(3)*r(3))/(bRadius*bRadius))
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
!	write(*,'(a,3i4,3f8.1)') 'tooNearBdry: ',site,prox,eradius,d
    return
endif
tooNearBdry = .false.
end function

!-----------------------------------------------------------------------------------------
! Is site near an exit?  prox is the minimum separation (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearExit(site,prox)
integer :: site(3)
real :: prox
integer :: iexit
real :: r(3), d

if (Nexits == 0) then
    tooNearExit = .false.
    return
endif
do iexit = 1,Lastexit
    if (exitlist(iexit)%ID == 0) cycle  ! exit no longer exists
    r = site - exitlist(iexit)%site
    d = norm(r)
    if (d < prox) then
        tooNearExit = .true.
        return
    endif
enddo
tooNearExit = .false.
end function

!-----------------------------------------------------------------------------------------
! Locate a free slot in a site adjacent to site1: site2 (freeslot)
! Returns freeslot = 0 if there is no free space in an adjacent site.
! The occupancy array occ() can be either occupancy() or big_occupancy(), hence
! the need for xlim (= NXX or NX)
!-----------------------------------------------------------------------------------------
subroutine get_free_slot(xlim,site1,site2,freeslot)
!type(occupancy_type) :: occ(:,:,:)
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
!subroutine showsite(node,site)
!integer :: node, site(3)
!
!write(*,*) 'showsite: ',site
!write(*,*) 'indx: ',occupancy(site(1),site(2),site(3))%indx
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function cell_count()
integer :: kcell, ntot
ntot = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle             ! skip gaps in the list
    ntot = ntot + 1
enddo
cell_count = ntot
end function

!-----------------------------------------------------------------------------------------
! The type (cognate or non-cognate) of a T cell can be random or strictly determined.
!-----------------------------------------------------------------------------------------
integer function select_cell_type(kpar)
integer :: kpar
!integer, save :: k = 0
integer :: nratio

if (random_cognate) then
    select_cell_type = random_cell_type(kpar)
else
!    nratio = 1./BC_COGNATE_FRACTION
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

if (.not.steadystate) then     ! surrogate for modeling an immune response
    call generate_traffic(inflow0)
!    if (NBcells < NBcells0) then
!		InflowTotal = InflowTotal*real(NBcells0)/NBcells
!	endif
else
    InflowTotal = inflow0
    OutflowTotal = inflow0
endif
if (istep == 1 .and. .not.use_TCP) then
	write(logmsg,'(a,i8,12f8.2)') 'NBcells,Inflow,Outflow: ',NBcells0,InflowTotal,OutflowTotal
    call logger(logmsg)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Total T cell inflow and outflow are generated from the vascularity and baseline
! inflow, inflow0.
!-----------------------------------------------------------------------------------------
subroutine generate_traffic(inflow0)
real :: inflow0
real :: act, expansion, actfactor, tnow
real :: inflow, outflow
!real, parameter :: T1 = 8*60, T2 = 16*60, T3 = 24*60


!traffic_mode == TRAFFIC_MODE_2 or TRAFFIC_MODE_3
! Note: if inflammation signal = 0 the vascularity (and inflow) should be constant
tnow = istep*DELTA_T
inflow = inflow0*Vascularity   ! level of vascularity (1 = steady-state)
outflow = NBcells*DELTA_T/(RESIDENCE_TIME*60)
InflowTotal = inflow
OutflowTotal = outflow
!if (mod(istep,240) == 0) then
!	write(logmsg,*) 'generate_traffic: inflow: ',inflow0,Vascularity,InflowTotal
!	call logger(logmsg)
!endif
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
! Display the properties of a cognate cell.  The calling program must ascertain that kcell
! is cognate.
!-----------------------------------------------------------------------------------------
subroutine show_cognate_cell(kcell)
integer :: kcell
type (cog_type), pointer :: p
!integer :: cogID
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
!gen = get_generation(p)
gen = p%generation
!call get_stage(p,stage,region)
stage = get_stage(p)
write(logmsg,'(a,i8,a,i2,a,i2)') '   cogID: ',p%cogID,' gen: ',gen,' stage: ', stage
call logger(logmsg)
write(logmsg,'(a,4f10.2)') '   times: entry,die,div,stage: ',bcell%entrytime,p%dietime,p%dividetime,p%stagetime
call logger(logmsg)
write(logmsg,'(a,3f8.2)') 'avidity, stimulation: ', p%avidity,p%stimulation
call logger(logmsg)

end subroutine

!-----------------------------------------------------------------------------------------
! Called whenever balancer carries out add_sites or removeSites.
! cognate_list(k) = 0 when the kth cognate cell has gone (left or died).
! If Mnodes = 1 this is called only once, after placeCells.  After that the list
! is maintained directly when cognate cells arrive or leave.
!-----------------------------------------------------------------------------------------
subroutine make_cognate_list(ok)
logical :: ok
integer :: kcell, ctype, stype, cogID
type (cog_type), pointer :: p

!write(*,*) 'make_cognate_list: ', lastcogID
ok = .true.
cognate_list(1:lastcogID) = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
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
!integer :: k, kcell, site(3)
character*(8) :: msg
integer :: error

!write(nfout,*) 'STEP: ',istep,lastcogID
!do k = 1,lastcogID
!    kcell = cognate_list(k)
!    if (kcell > 0) then
!        site = cellist(kcell)%site
!        gen = get_generation(cellist(kcell)%cptr)
!    else
!        site = 0
!        gen = 0
!    endif
!    write(nfout,'(i8,3i4)') kcell,site,gen
!enddo
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
! We display only T cells that are still in the FOLLICLE
!--------------------------------------------------------------------------------
subroutine save_cell_positions
!!!use ifport
integer :: k, kcell, site(3), j
integer :: itcstate, stype, ctype, stage, region
real :: bcell_diam = 0.9
!real :: spectrum_max = 10, spectrum_freefraction = 0.9
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
!		call logger('wait')
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

	! T cell section
	do k = 1,lastcogID
		kcell = cognate_list(k)
		if (kcell > 0) then
!			call get_stage(cellist(kcell)%cptr,stage,region)
			region = get_region(cellist(kcell)%cptr)
			if (region /= FOLLICLE) cycle
			site = cellist(kcell)%site
!			gen = get_generation(cellist(kcell)%cptr)
			gen = cellist(kcell)%cptr%generation
			itcstate = gen
			! Need tcstate to convey non-activated status, i.e. 0 = non-activated
			write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',k-1, site, bcell_diam, itcstate
		endif
	enddo

write(nfpos,'(a2,i6)') 'E ',istep
close(nfpos)

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_exnode
integer :: k, kcell, site(3), nd, j
real :: tcstate
real :: bcell_diam = 0.9
real :: spectrum_max = 10, spectrum_freefraction = 0.6
integer :: gen, stage, region, bnd(2)
character*(64) :: fname = '\CMGUI\DCU\dcu_00000.exnode'

write(fname(16:20),'(i5.5)') istep
write(*,*) fname
open(nfcmgui,file=fname,status='replace')

nd = 0

! T cell section
write(nfcmgui,*) 'Group name : bcell'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
write(nfcmgui,*) '    value.'
write(nfcmgui,*) '  3) tcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  3) tcgen, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  4) tcbound, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
!		call get_stage(cellist(kcell)%cptr,stage,region)
		stage = get_stage(cellist(kcell)%cptr)
		region = get_region(cellist(kcell)%cptr)
		if (region /= FOLLICLE) cycle
        nd = nd+1
        site = cellist(kcell)%site
!        gen = get_generation(cellist(kcell)%cptr)
		gen = cellist(kcell)%cptr%generation
        tcstate = (gen-1.0)/(BC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
!        write(nfcmgui,'(a,i8,3i4,f4.1,i3,i2)') 'Node: ',nd, site, bcell_diam, gen, tcbound
        write(nfcmgui,'(a,i8,3i4,f4.1,f6.2)') 'Node: ',nd, site, bcell_diam, tcstate
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

!call random_number(R)
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
! Need to compare chemo_probs() with chemo_probs_pre().
! Good news!  They are in agreement.
! Note: originally vsum was the offset vector of the cell from the attracting site.  
! This has been changed.  Now vsum is the normalised net chemokine gradient vector,
! and the probabilities of jumps in directions close to the net gradient vector are
! greatest.  If the angle between a jump direction and the gradient vector exceeds
! pi/2 (i.e. the dot product (cosine) is negative) the jump probability is set to 0.
!--------------------------------------------------------------------------------
subroutine test_chemo
integer :: k,rv(3)
real(DP) :: p0(MAXRELDIR+1), p(MAXRELDIR+1)
real :: v(3), vsum(3), f

f = 1
p0 = 1

vsum = (/1., .0, .0/)
v = vsum/norm(vsum)
rv = chemo_N*v
p = p0
write(*,*) 'chemo_probs_pre'
call chemo_probs_pre(p,rv,f)     ! this is the precomputed version
write(*,'(7f8.4)') p
p = p0
write(*,*) 'chemo_probs'
call chemo_probs(p,v,f)
!write(*,'(7f8.4)') p
do k = 1,njumpdirs
	write(*,'(4i4,f8.4)') k,jumpvec(:,k),p(k)
enddo
end subroutine

!--------------------------------------------------------------------------------
! Original interpretation, in terms of exit chemotaxis:
! The chemotactic step probabilities for all possible sites within an exit's
! SOI are precomputed.  The array is indexed by the offset of each site from
! the attractant site, given by (x,y,z).
! Note that the weight given to a jump direction is found from the cosine^2 of
! the angle between the jump direction and the vector from the attractant site to the
! cell site.  The offset vector is used only to provide direction information.
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
!                if (cosa < 0) then
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
!dbug = (sum(p) < 0.001 .and. v(1) == 0)
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
!    if (cosa < 0) then
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


!--------------------------------------------------------------------------------
! Computes the functional dependence of chemotactic effect on distance r from
! the exit.  r is in units of site spacing.
!--------------------------------------------------------------------------------
real function chemo_g(r)
real :: r

chemo_g = min(1.0,(1.0/r)**chemo_exp)
end function

!--------------------------------------------------------------------------------
! Determines the degree to which a cell is subject to chemotaxis.
! (Determined by CD69 level, or S1P1 level.)
! If use_chemotaxis is true, i.e. chemotaxis is used to control cell exit, any cell
! that gets close enough to the exit will leave the paracortex.  Is this
! acceptable?
! For a noncognate cell, should rise from 0 to 1 in an hour or so.
!--------------------------------------------------------------------------------
real function chemo_active_exit(cell)
type(cell_type), pointer :: cell
real :: tnow, t

if (turn_off_chemotaxis) then
    chemo_active_exit = 0
    return
endif

if (TAGGED_EXIT_CHEMOTAXIS) then
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
    if (cell%ctype == RES_TAGGED_CELL) then     ! testing effect of S1P1
        chemo_active_exit = (1 - exp(-K1_S1P1*t))*TAGGED_CHEMO_ACTIVITY
    else
		chemo_active_exit = 0
	endif
	return
endif

!if (associated(cell%cptr)) then     ! cognate cell
!    chemo_active_exit = cell%cptr%S1P1
!else
!    tnow = istep*DELTA_T
!    t = tnow - cell%entrytime
!    chemo_active_exit = 1 - exp(-K1_S1P1*t)
!endif
end function

!--------------------------------------------------------------------------------
! Returns the level of CCR7 ligand (i.e. CCL19/21) at a distance r sites from an
! exit site, within the exit SOI.  The value must be in the range (0,1).
! The parameters CCR7_R1 and CCR7_R2 specify the range of r over which the
! ligand level ranges linearly from 0 to 1.  This is a simple way to program a
! decreased level of CCR7 near an exit, following Cyster's ideas.
! Note: currently this always returns 1
!--------------------------------------------------------------------------------
real function CCR7_ligand(r)
real :: r
real, parameter :: CCR7_R1 = 0, CCR7_R2 = 0

if (r < CCR7_R1) then
    CCR7_ligand = 0
elseif (r < CCR7_R2) then
    CCR7_ligand = (r - CCR7_R1)/(CCR7_R2 - CCR7_R1)
else
    CCR7_ligand = 1
endif
end function

!--------------------------------------------------------------------------------
! Determines e(:), the location of the nearest exit, if there is a close one.
! Otherwise returns e(:) = 0.
! Note: currently only one exit number is stored in occupancy()
!--------------------------------------------------------------------------------
subroutine nearest_exit(site,in_exit_SOI,e)
integer :: site(3), e(3)
logical :: in_exit_SOI
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
    in_exit_SOI = .false.
    e = (/0,0,0/)
    return
endif
in_exit_SOI = .true.
e = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
end subroutine

!--------------------------------------------------------------------------------
! The criterion for a near exit is based on chemo_radius
!--------------------------------------------------------------------------------
subroutine near_exits(site,ne,ee)
integer :: site(3), ne, ee(3,*)
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
	ne = 0
    return
endif
ne = 1
ee(:,1) = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
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
    if (cellist(kcell)%ID == 0) cycle
    xyzsum = xyzsum + cellist(kcell)%site
enddo
write(nfres,'(2i6,3i12)') istep,k,xyzsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkExits(msg)
character*(*) :: msg
integer :: iexit,site(3)
logical :: ok = .true.

write(logmsg,'(a,i8,a,a)') 'checkExits: ',istep,'  ',msg
call logger(logmsg)
do iexit = 1,Lastexit
    if (exitlist(iexit)%ID == 0) cycle
	site = exitlist(iexit)%site
	if (occupancy(site(1),site(2),site(3))%exitnum /= -exitlist(iexit)%ID) then
		write(logmsg,*) 'checkExits: ',iexit,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
		ok = .false.
	endif
enddo
if (.not.ok) stop
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_exit(iexit)
integer :: iexit
integer :: x,y,z,a,nc,nctot(0:10),site(3)
integer :: n = 5

site = exitlist(iexit)%site
nctot(0) = sitecells(site,0,0,0)
do a = 1,n
    nc = 0
    do z = -a,a,2*a
        do x = -a,a
            do y = -a,a
                nc = nc + sitecells(site,x,y,z)
            enddo
        enddo
    enddo
    do z = -a+1,a-1
        y = -a
        do x = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = a
        do y = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        y = a
        do x = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = -a
        do y = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
    enddo
    nctot(a) = nc
enddo
write(nfout,'(10i6)') nctot(0:n),sum(nctot(0:n))
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
    if (cellist(kcell)%ID == 0) cycle
    site = cellist(kcell)%site
    tot(site(3)) = tot(site(3)) + 1
enddo
tot = tot/nlist
!write(*,'(a,i8,2f6.3)') 'z distribution: ',nlist,tot(NZ/2+10),tot(NZ/2-10)
!write(*,'(10f6.3)') tot
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
            write(*,*) 'Bad tagged cell: ',site,d2,aRadius*aRadius
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
integer :: ic,x,y,z,slots,slot,cells,k,indx(2),site(3)
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

do ic = 1,nlist
    if (cellist(ic)%ID == 0) cycle  ! gap
    site = cellist(ic)%site
    indx = occupancy(site(1),site(2),site(3))%indx
    if (ic == indx(1)) then
        slot = 1
    elseif (ic == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: checker: bad indx: ',ic,site,indx
        stop
    endif
enddo
deallocate(tot)
write(*,*) 'checked OK: ',' nlist: ',nlist
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
type (cog_type), pointer :: p
integer :: kcell, ctype, n1, n2

n1 = 0
n2 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
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
