! B cell behaviour
module behaviour

use global
use fields
use motility

implicit none

!INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )


!integer, parameter :: POST_DIVISION = SWARMS    ! on division, cells enter SWARMS

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

!integer, parameter :: Nbindtime_nc = 10
!real, parameter :: CYTOKINESIS_TIME = 20.0

type(dist_type), allocatable :: stage_dist(:,:,:), life_dist(:), divide_dist(:)

contains

!--------------------------------------------------------------------------------------
! Compute parameters for probability distributions of lifetime and time-to-divide
! In the case of NORMAL distributions, p1 = mean, p2 = std dev
! For X LOGNORMAL distributed, p1 and p2 are mean and std dev of Y, and X = exp(Y)
! Lifetime and dividetime parameters are given in hours, then converted to minutes
! when used.  Program time is in minutes.
!--------------------------------------------------------------------------------------
subroutine setup_dists
!real :: divide_median1, divide_median2
!real :: divide_shape1, divide_shape2
!real :: life_tc1 = 15, life_tc2 = 18
!real :: life_median1 = 48, life_median2 = 24
!real :: life_median1 = 196, life_median2 = 196      !<------------- Note: hard-coded values
!real :: life_shape1 = 1.5, life_shape2 = 1.4
integer :: i

!life_dist(1)%class = EXPONENTIAL_DIST
!life_dist(i)%p1 = life_tc1
!do i = 1,TC_MAX_GEN
!	life_dist(i)%class = EXPONENTIAL_DIST
!	life_dist(i)%p1 = life_tc2
!enddo

allocate(life_dist(BC_MAX_GEN))        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(divide_dist(BC_MAX_GEN))      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

life_dist(1)%class = LOGNORMAL_DIST
life_dist(1)%p1 = log(60*BC_life_median1)
life_dist(1)%p2 = log(BC_life_shape)
do i = 2,BC_MAX_GEN
	life_dist(i)%class = LOGNORMAL_DIST
	life_dist(i)%p1 = log(60*BC_life_median2)
	life_dist(i)%p2 = log(BC_life_shape)
enddo


divide_dist(1)%class = divide_dist1%class
divide_dist(1)%p1 = divide_dist1%p1
divide_dist(1)%p2 = divide_dist1%p2
divide_dist(2)%class = divide_dist2%class
divide_dist(2)%p1 = divide_dist2%p1
divide_dist(2)%p2 = divide_dist2%p2
do i = 3,BC_MAX_GEN
	divide_dist(i)%class = divide_dist(2)%class
	divide_dist(i)%p1 = divide_dist(2)%p1
	divide_dist(i)%p2 = divide_dist(2)%p2
enddo

end subroutine

!--------------------------------------------------------------------------------------
! If a cognate cell is NAIVE it is long-lived - effectively it does not die.
! When a cognate cell has contact with a cognate DC its lifetime is immediately limited.
!--------------------------------------------------------------------------------------
real function BClifetime(ptr)
type (cog_type), pointer :: ptr
integer :: gen, stage, region
real :: p1, p2
integer :: kpar = 0

BClifetime = BIG_TIME
return

!stage = get_stage(ptr)
call get_stage(ptr,stage,region)
if (stage == NAIVE) then
    BClifetime = BIG_TIME
    return
endif
gen = get_generation(ptr)
if (gen < 1 .or. gen > BC_MAX_GEN) then
    write(logmsg,*) 'BClifetime: bad gen: ',gen
    call logger(logmsg)
!    stop
endif
p1 = life_dist(gen)%p1
p2 = life_dist(gen)%p2
select case (life_dist(gen)%class)
case (NORMAL_DIST)
	BClifetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	BClifetime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	BClifetime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function dividetime(gen,ictype)
integer :: gen,ictype
real :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(gen)%p1
p2 = divide_dist(gen)%p2
!write(*,*) 'dividetime: ',istep,gen,divide_dist(gen)%p1,divide_dist(gen)%p2
select case (divide_dist(gen)%class)
case (NORMAL_DIST)
	dividetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
!    write(*,*) 'dividetime: ',istep,gen,p1,p2
	dividetime = rv_lognormal(p1,p2,kpar)
!	write(*,*) 'dividetime: ',dividetime
case (CONSTANT_DIST)
	dividetime = p1
end select
!if (ictype == CD8_CELL) then
!	dividetime = dividetime*CD8_DIVTIME_FACTOR
!endif
end function


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(DP) function rv_normal(p1,p2,kpar)
integer :: kpar
real :: p1,p2
real(DP) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(DP) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real :: p1,p2
real(DP) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(DP) function my_rnor()
real(DP) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
!    call random_number(R)
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(DP) function rv_exponential(p1)
real :: p1
real(DP) :: r
integer :: kpar = 0

r = par_rexp(kpar)
rv_exponential = p1*r
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(DP) function cum_prob_lognormal(a,p1,p2)
real :: a, p1, p2
real(DP) :: b, prob

!p1 = log(m)
!p2 = log(s)
b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!--------------------------------------------------------------------------------
! The probability histogram for a (supposedly) lognormal variable X is given
! by the probability values p(i), i=1,..,n and the associated interval
! delimiting values x(i), i=1,..,n+1
! where p(i) = Pr{x(i) <= X < x(i+1)}
! The cumulative probability for a lognormal variate with median m, shape s is
! Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------
subroutine lognfit(n,x,p,mrange,srange,mean,shape)
integer :: n
real :: x(*), p(*), mrange(2), srange(2)
real :: mean, shape
integer :: i
real :: alo, ahi, p1, p2, c, err

p1 = log(mrange(1))
p2 = log(srange(1))
err = 0
do i = 1,n
	alo = x(i)
	ahi = x(i+1)
	c = cum_prob_lognormal(ahi,p1,p2) - cum_prob_lognormal(alo,p1,p2)
	err = err + (c - p(i))**2
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_fixed_BCparams(ok)
logical :: ok
logical :: ext

inquire(file=fixedfile,exist=ext)
if (.not.ext) then
	call logger("Fixed parameter input file not found")
	ok = .false.
	return
endif
open(nfcell,file=fixedfile,status='old',err=99)
ok = .true.
read(nfcell,*) BC_life_median1		!= 196				! median lifetime of naive B cells
read(nfcell,*) BC_life_median2		!= 196				! median lifetime of activated B cells
read(nfcell,*) BC_life_shape		!= 1.4				! shape parameter for lifetime of B cells
read(nfcell,*) NBC_LN				!= 3.0e07			! number of B cells in a LN
read(nfcell,*) NBC_BODY				!= 1.6e09			! number of circulating B cells in the whole body
read(nfcell,*) NLN_RESPONSE								! number of LNs in the response
read(nfcell,*) BC_RADIUS					            ! radius of B cell (um)
read(nfcell,*) BC_MAX_GEN				                ! maximum B cell generation
read(nfcell,*) divide_dist1%class
read(nfcell,*) divide_dist2%class
read(nfcell,*) GAMMA						            ! controls crowding

close(nfcell)
return
99	continue
ok = .false.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_Bcell_params(ok)
logical :: ok
real :: sigma, divide_mean1, divide_shape1, divide_mean2, divide_shape2
integer :: i, shownoncog, ncpu_dummy, iuse(MAX_CHEMO)
integer :: usetraffic, usechemo, computedoutflow
character(4) :: logstr

ok = .true.
!write(*,*) 'Read cell parameter file: ',inputfile
!call logger('Reading cell parameter file')
!call logger(inputfile)
open(nfcell,file=inputfile,status='old')
read(nfcell,*) BC_AVIDITY_MEDIAN            ! median of avidity distribution (only if fix_avidity = false)
read(nfcell,*) BC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
read(nfcell,*) BC_COGNATE_FRACTION			
read(nfcell,*) BC_STIM_RATE_CONSTANT		! rate const for TCR stimulation (-> molecules/min)
read(nfcell,*) BC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
read(nfcell,*) divide_mean1
read(nfcell,*) divide_shape1
read(nfcell,*) divide_mean2
read(nfcell,*) divide_shape2
read(nfcell,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfcell,*) RHO							! persistence: 0 < rho < 1	(0.95)
read(nfcell,*) NX							! rule of thumb: about 4*BLOB_RADIUS
read(nfcell,*) BLOB_RADIUS					! initial T cell blob size (sites)
read(nfcell,*) BC_FRACTION					! B cell fraction of follicle volume
read(nfcell,*) FLUID_FRACTION				! fraction of paracortex that is fluid
read(nfcell,*) usetraffic					! use T cell trafficking
read(nfcell,*) usechemo                     ! use chemotaxis
read(nfcell,*) computedoutflow				! compute outflow (with inflow)
read(nfcell,*) RESIDENCE_TIME               ! T cell residence time in hours -> inflow rate
! Vascularity parameters
read(nfcell,*) Inflammation_days1	        ! Days of plateau level - parameters for VEGF_MODEL = 1
read(nfcell,*) Inflammation_days2	        ! End of inflammation
read(nfcell,*) Inflammation_level	        ! This is the level of inflammation
read(nfcell,*) iuse(S1P)
read(nfcell,*) chemo(S1P)%bdry_conc
read(nfcell,*) chemo(S1P)%diff_coeff
read(nfcell,*) chemo(S1P)%halflife
read(nfcell,*) chemo(S1P)%strength
read(nfcell,*) iuse(CCL21)
read(nfcell,*) chemo(CCL21)%bdry_conc
read(nfcell,*) chemo(CCL21)%diff_coeff
read(nfcell,*) chemo(CCL21)%halflife
read(nfcell,*) chemo(CCL21)%strength
read(nfcell,*) iuse(OXY)
read(nfcell,*) chemo(OXY)%bdry_conc
read(nfcell,*) chemo(OXY)%diff_coeff
read(nfcell,*) chemo(OXY)%halflife
read(nfcell,*) chemo(OXY)%strength
read(nfcell,*) iuse(CXCL13)
read(nfcell,*) chemo(CXCL13)%bdry_conc
read(nfcell,*) chemo(CXCL13)%diff_coeff
read(nfcell,*) chemo(CXCL13)%halflife
read(nfcell,*) chemo(CXCL13)%strength
!read(nfcell,*) VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity
read(nfcell,*) chemo_radius			        ! radius of chemotactic influence (um)
read(nfcell,*) base_exit_prob               ! base probability of exit at a boundary site
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) fixedfile					! file with "fixed" parameter values
close(nfcell)

! Chemokines
use_chemotaxis = .true.
do i = 1,MAX_CHEMO
	chemo(i)%used = (iuse(i) == 1)
	chemo(i)%decay_rate = DecayRate(chemo(i)%halflife)
enddo
chemo(S1P)%name    = 'S1P'
chemo(CCL21)%name  = 'CCL21'
chemo(OXY)%name    = 'OXY'
chemo(CXCL13)%name = 'CXCL13'
use_S1P    = chemo(S1P)%used
use_CCL21  = chemo(CCL21)%used
use_OXY    = chemo(OXY)%used
use_CXCL13 = chemo(CXCL13)%used

VEGF_MODEL = 1
if (usetraffic == 1) then
	USE_TRAFFIC = .true.
else
	USE_TRAFFIC = .false.
endif
!if (usechemo == 1) then
!	use_chemotaxis = .true.
!	! Interim measure:
!	chemo_K = 1.0
!else
!	use_chemotaxis = .false.
!	chemo_K_exit = 0
!endif
if (computedoutflow == 1) then
	computed_outflow = .true.
else
	computed_outflow = .false.
endif
!exit_region = EXIT_LOWER_SURFACE
use_portal_egress = .false.

call read_fixed_BCparams(ok)
if (.not.ok) then
	write(logmsg,'(a,a)') 'Error reading fixed input data file: ',fixedfile
	call logger(logmsg)
	return
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
if (BC_MAX_GEN > MMAX_GEN) then
    write(logmsg,*) 'BC_MAX_GEN > MMAX_GEN: ',BC_MAX_GEN,MMAX_GEN
    call logger(logmsg)
    ok = .false.
    return
endif
PI = 4*atan(1.0)
DELTA_X = (4*PI/(3*BC_FRACTION))**0.33333*BC_RADIUS
open(nfout,file=outputfile,status='replace')
write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

chemo_radius = chemo_radius/DELTA_X				! convert from um to lattice grids
chemo_N = max(3,int(chemo_radius + 0.5))	    ! convert from um to lattice grids
chemo_exp = log(1/CHEMO_MIN)/log(chemo_radius)

sigma = log(divide_shape1)
divide_dist1%p1 = log(60*divide_mean1/exp(sigma*sigma/2))
divide_dist1%p2 = sigma
sigma = log(divide_shape2)
divide_dist2%p1 = log(60*divide_mean2/exp(sigma*sigma/2))
divide_dist2%p2 = sigma

if (BC_COGNATE_FRACTION == 0) then
    use_cognate = .false.
else
    use_cognate = .true.
endif

if (BC_STIM_HALFLIFE > 0) then
!    BCRdecayrate = log(2.0)/(BC_STIM_HALFLIFE*60)    ! rate/min
    BCRdecayrate = DecayRate(BC_STIM_HALFLIFE)    ! rate/min
else
    BCRdecayrate = 0
endif

if (VEGF_MODEL == 0) then
    vary_vascularity = .false.
else
    vary_vascularity = .true.
endif

Ve = DELTA_X*DELTA_X*DELTA_X
Vc = FLUID_FRACTION*Ve

Nsteps = days*60*24/DELTA_T

call setup_dists

ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
! Convert halflife in hours to a decay rate /min
!----------------------------------------------------------------------------------------
real function DecayRate(halflife)
real :: halflife

DecayRate = log(2.0)/(halflife*60)
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine save_inputfile(cellfile)
character(LEN=*) :: cellfile
character*(128) :: line

open(nfcell,file=cellfile,status='old')
do
    read(nfcell,'(a)',end=999) line
    write(nfout,'(a)') line
enddo
999 continue
close(nfcell)
end subroutine

!----------------------------------------------------------------------------------------
! Save parameters (values hard-coded, not yet in input file) 
!----------------------------------------------------------------------------------------
subroutine save_parameters

write(nfout,*) 'PARAMETER data'
write(nfout,*) '--------------'
write(nfout,*) 'DELTA_X: ',DELTA_X
write(nfout,*) 'BALANCER_INTERVAL: ',BALANCER_INTERVAL
write(nfout,*) 'use_add_count: ',use_add_count
write(nfout,*) 'use_traffic: ',use_traffic
write(nfout,*) 'use_chemotaxis: ',use_chemotaxis
write(nfout,*) 'computed_outflow: ',computed_outflow
write(nfout,*) 'use_cytokines: ',use_cytokines
write(nfout,*) 'use_cognate: ',use_cognate
write(nfout,*) 'random_cognate: ',random_cognate
write(nfout,*) 'use_diffusion: ',use_diffusion
write(nfout,*) 'NGEN_EXIT: ',NGEN_EXIT
write(nfout,*) 'NDIFFSTEPS: ',NDIFFSTEPS
write(nfout,*) 'VEGF_MODEL: ',VEGF_MODEL
write(nfout,*) 'VEGF_alpha: ',VEGF_alpha
write(nfout,*) 'VEGF_beta: ',VEGF_beta
write(nfout,*) 'VEGF_decayrate: ',VEGF_decayrate
write(nfout,*) 'vasc_maxrate: ',vasc_maxrate
write(nfout,*) 'vasc_decayrate: ',vasc_decayrate
write(nfout,*) 'vasc_beta: ',vasc_beta
write(nfout,*) 'vasc_n: ',vasc_n
write(nfout,*) 'fix_avidity: ',fix_avidity
write(nfout,*) 'avidity_nlevels: ',avidity_nlevels
write(nfout,*) 'avidity_logscale: ',avidity_logscale
write(nfout,*) 'avidity_min: ',avidity_min
write(nfout,*) 'avidity_step: ',avidity_step

end subroutine

!-----------------------------------------------------------------------------------------
! When a T cell dies it is removed from the cell list (%ID -> 0)
! and removed from occupancy()%indx()
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine Bcell_death(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2), ctype, stype, region
logical :: cognate

write(logmsg,*) 'Bcell_death: ',kcell
call logger(logmsg)
cognate = (associated(cellist(kcell)%cptr))
if (cognate) then
	call get_region(cellist(kcell)%cptr,region)
endif

cellist(kcell)%ID = 0
totalres%dN_Dead = totalres%dN_Dead + 1
totalres%N_Dead = totalres%N_Dead + 1
ctype = cellist(kcell)%ctype
stype = struct_type(ctype)
if (stype == COG_TYPE_TAG) then
    cognate_list(cellist(kcell)%cptr%cogID) = 0
endif
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,*) 'Bcell_death: ngaps > max_ngaps'
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell

NBcells = NBcells - 1
site = cellist(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(1) = indx(2)
    occupancy(site(1),site(2),site(3))%indx(2) = 0
elseif (indx(2) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(2) = 0
else
    write(logmsg,*) 'ERROR: Bcell_death: cell not at site: ',kcell,site,indx
    call logger(logmsg)
    stop
endif

! Now we need to make a site unavailable
! Remove (make OUTSIDE) a boundary site near to the specified site.
! The criteria for selection of a site to remove are that it is on the blob boundary,
! preferably "sticking out", and that it is vacant.  A site may be made vacant by 
! moving a cell to a nearby available site.
! One way to handle this is to maintain a count of the number of sites to be added(removed).
! At regular intervals the counts can be aggregated and if the total is big enough (+ or -)
! they can all be processed within a gather-scatter cycle.
if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(logmsg,*) 'Bcell_death: No site removal code, use add count'
    call logger(logmsg)
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
! A cognate cell divides.  It has already been determined that there is space for the extra cell.
! One cell retains the cogID, the other gets cogID = 0 for now, to be allocated
! sequentially later, after gather_data() (unless Mnodes = 1, in which case the
! allocation of cogID is done in CreateBcell.)
! The other cell is placed in freeslot of site2 (which may be the same site as kcell).
!-----------------------------------------------------------------------------------------
subroutine cell_division(kcell,site2,freeslot,ok)
integer :: kcell, site2(3), freeslot
logical :: ok
integer :: icnew, ctype, gen, region, site(3), indx(2)
integer :: iseq, tag, kfrom, kto
real :: tnow
type(cog_type), pointer :: p1, p2
!real :: IL_state(CYT_NP)
integer :: kpar = 0
real(DP) :: R
logical :: germinal

write(*,*) 'cell_division: ',kcell
ok = .true.
tnow = istep*DELTA_T
!call show_cognate_cell(kcell)
p1 => cellist(kcell)%cptr
gen = get_generation(p1)
call get_region(p1,region)
!write(*,*) 'cell_division: ',kcell,gen,istep,tnow
if (gen == BC_MAX_GEN) then
    write(logmsg,*) 'cell_division: reached maximum generation: ',kcell
	call logger(logmsg)
    return
endif
if (gen == 1) then
    write(logmsg,*) 'First division at hour: ',tnow/60
    call logger(logmsg)
else
!    write(*,*) 'Time for next division: ',kcell,gen,(tnow-cellist(kcell)%entrytime)/60
endif
ndivided(gen) = ndivided(gen) + 1
tdivided(gen) = tdivided(gen) + (tnow - p1%dividetime)

if (ngaps > 0) then
    icnew = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: cell_division: cell list full: ',nlist
		call logger(logmsg)
		ok = .false.
		return
	endif
    icnew = nlist
endif
cellist(kcell)%lastdir = random_int(1,6,kpar)
gen = gen + 1
call set_generation(p1,gen)
germinal = p1%germinal
if (.not.germinal) then
	R = par_uni(kpar)	
	if (R < GCC_PROB) then
		p1%germinal = .true.
		call set_stage(p1,BCL6_UP)
		p1%stagetime = tnow + T_BCL6_UP
	else
		call set_stage(p1,DIVIDING)
		p1%stagetime = tnow + T_DIVISION
	endif
else
	call set_stage(p1,DIVIDING)
	p1%stagetime = tnow + T_DIVISION
endif
!call set_stage(p1,POST_DIVISION)
!if (TCR_splitting) then
!    p1%stimulation = p1%stimulation/2
!endif
ctype = cellist(kcell)%ctype
p1%dietime = tnow + BClifetime(p1)
p1%dividetime = tnow
p1%stagetime = tnow + dividetime(gen,ctype)
site = cellist(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx

call CreateBcell(icnew,cellist(icnew),site2,ctype,gen,DIVIDING,region,ok)
if (.not.ok) return

p2 => cellist(icnew)%cptr

p2%stimulation = p1%stimulation
p2%status = p1%status
cellist(icnew)%entrytime = tnow
if (.not.germinal) then
	R = par_uni(kpar)	
	if (R < GCC_PROB) then
		p2%germinal = .true.
		call set_stage(p2,BCL6_UP)
		p2%stagetime = tnow + T_BCL6_UP
	else
		call set_stage(p2,DIVIDING)
		p2%stagetime = tnow + T_DIVISION
	endif
else
	call set_stage(p2,DIVIDING)
	p2%stagetime = tnow + T_DIVISION
endif
ndivisions = ndivisions + 1
occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = icnew
if (use_add_count) then
    nadd_sites = nadd_sites + 1
else
    write(logmsg,*) 'cell_division: No site removal code, use add count'
	call logger(logmsg)
	ok = .false.
    return
endif
NBcells = NBcells + 1
end subroutine

!-----------------------------------------------------------------------------------------
! Locate a free slot in a site adjacent to site1: site2 (freeslot)
! Returns freeslot = 0 if there is no free space in an adjacent site.
! The occupancy array occ() can be either occupancy() or big_occupancy(), hence
! the need for xlim (= NXX or NX)
!-----------------------------------------------------------------------------------------
subroutine get_free_slot(occ,xlim,site1,site2,freeslot)
type(occupancy_type) :: occ(:,:,:)
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
	indx2 = occ(site2(1),site2(2),site2(3))%indx
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


!-----------------------------------------------------------------------------------------
! Create a new B cell at site
! We could give a progeny cell (the result of cell division) the same ID as its parent.
!-----------------------------------------------------------------------------------------
subroutine CreateBcell(kcell,cell,site,ctype,gen,stage,region,ok)
type(cell_type) :: cell
integer :: kcell, site(3), ctype, gen, stage, region
logical :: ok
integer :: stype, cogID, i
real :: tnow, param1, param2
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
stype = struct_type(ctype)
!write(*,*) 'CreateBcell: ',kcell,ctype,stype
cell%entrytime = tnow
if (stype == NONCOG_TYPE_TAG) then
    if (dbug) then
        write(*,*) 'CreateBcell: ',kcell,site,associated(cell%cptr)
    endif
    if (associated(cell%cptr)) then
        deallocate(cell%cptr)
    endif
elseif (stype == COG_TYPE_TAG) then
    if (.not.associated(cell%cptr)) then
        allocate(cell%cptr)
    endif
    param1 = log(BC_AVIDITY_MEDIAN)
    param2 = log(BC_AVIDITY_SHAPE)
    cell%cptr%status = 0
    call set_generation(cell%cptr,gen)
    call set_stage_region(cell%cptr,stage,region)
    if (fix_avidity) then
        i = mod(navid,avidity_nlevels)
        navid = navid + 1
        cell%cptr%avidity = avidity_level(i+1)
    else
        cell%cptr%avidity = rv_lognormal(param1,param2,kpar)
    endif
    cell%cptr%stimulation = 0
    cell%cptr%germinal = .false.
!    if (use_cytokines) then
!        call IL2_init_state(cell%cptr%IL_state,cell%cptr%IL_statep)
!    endif
    ! What should the initial CD69 level be?  If 0, can a cognate cell exit immediately?
    ! We would prefer not to impose a time or generation constraint on the exit of
    ! cognate T cells, but otherwise if CD69 is initially 0, a cell will be susceptible
    ! to chemotaxis and exit until it has received enough TCR signal to drive CD69 high.
!    cell%cptr%CD69 = 0
!    cell%cptr%S1P1 = 0
    cell%cptr%dietime = tnow + BClifetime(cell%cptr)
    cell%cptr%dividetime = tnow
    cell%cptr%stagetime = BIG_TIME

! Maintain cognate_list at start or if we are running on a single node
! Otherwise cogID and cognate_list is maintained by make_cognate_list
!    if (istep == 0 .or. Mnodes == 1) then
        lastcogID = lastcogID + 1
        if (lastcogID > MAX_COG) then
            write(logmsg,'(a,i6)') 'Error: CreateBcell: cognate_list dimension exceeded: ',MAX_COG
            call logger(logmsg)
            ok = .false.
            return
        endif
        cogID = lastcogID
        cell%cptr%cogID = cogID
        cognate_list(cogID) = kcell
!    else
!        cell%cptr%cogID = 0
!    endif
else
    write(logmsg,*) 'ERROR: CreateBcell: bad ctype: ',ctype
    call logger(logmsg)
    stop
endif
! Interim measure:
cell%receptor = 1
!if (use_S1P) then
!    cell%receptor(S1P) = 10.0
!endif
!if (use_CCL21) then
!    cell%receptor(CCL21) = 10.0
!endif
!if (use_OXY) then
!    cell%receptor(OXY) = 0.4
!endif
!if (use_CXCL13) then
!    cell%receptor(CXCL13) = 0.4
!endif
lastID = lastID + 1     ! Each node makes its own numbers, with staggered offset
cell%ID = lastID
cell%site = site
cell%ctype = ctype
cell%step = 0
cell%lastdir = random_int(1,6,kpar)
end subroutine

!--------------------------------------------------------------------------------
! Add a cell (kcell) with characteristics (ctype, gen, stage) at site.
!--------------------------------------------------------------------------------
subroutine AddBcell(site,ctype,gen,stage,region,kcell,ok)
integer :: site(3), ctype, gen, stage, region, kcell
logical :: ok
integer :: indx(2)

ok = .true.
if (ngaps > 0) then
    kcell = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: add_Bcell: cell list full: ',nlist
		call logger(logmsg)
		ok = .false.
		return
	endif
    kcell = nlist
endif
if (dbug) then
    write(*,'(a,9i7)') 'AddBcell: ',istep,kcell,site,ctype,gen,stage,region
endif
call CreateBcell(kcell,cellist(kcell),site,ctype,gen,stage,region,ok)
if (.not.ok) return

indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == 0) then
    indx(1) = kcell
elseif (indx(2) == 0) then
    indx(2) = kcell
else
    write(logmsg,*) 'ERROR: add_Bcell: no free slot: ',site,indx
    call logger(logmsg)
    stop
endif
occupancy(site(1),site(2),site(3))%indx = indx
NBcells = NBcells + 1

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine add_site_local(site)
integer :: site(3)

if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(logmsg,*) 'add_site_local: no code to add site, use add count'
    call logger(logmsg)
    stop
endif
end subroutine

!--------------------------------------------------------------------------------
! Add a vacant site at a boundary to account for the T cell added at site
! In the case of an end node (me = 0 or me = Mnodes-1) the added site can be
! at a different x value, but for internal nodes the added site  must go at
! (or near) the same x value.
! For the spherical blob case, try to place the site at a point near the bdry
! on a line drawn through site from the centre.
! NOT USED
!--------------------------------------------------------------------------------
subroutine add_vacant_site(site,kpar)
integer :: site(3),kpar
integer :: k, site0(3),newsite(3)
real(DP) :: R
real :: dxyz(3)
logical :: redo

if (dbug) write(*,'(a,4i6)') 'add_vacant_site: ',site
site0 = site
dxyz = real(site0) - Centre
do k = 2,3
!    call random_number(R)
    R = par_uni(kpar)
    dxyz(k) = dxyz(k) + (R - 0.5)
enddo
call normalize(dxyz)
redo = .false.
k = 0
do
    k = k+1
    newsite = site0 + k*0.5*dxyz
    if (newsite(1) < 1 .or. newsite(1) > NX) then
        site0 = site0 + (k-1)*0.5*dxyz
        dxyz(1) = 0
        call normalize(dxyz)
        redo = .true.
        exit
    endif
    if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
        write(*,*) 'ERROR: add_vacant_site: reached grid limits (a): ',k,site,dxyz
        stop
    endif
    if (occupancy(newsite(1),newsite(2),newsite(3))%indx(1) == OUTSIDE_TAG) then
        exit
    endif
enddo
if (redo) then
    if (dbug) write(*,*) 'redo: ', site0
    k = 0
    do
        k = k+1
        newsite = site0 + k*0.5*dxyz
        if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
            write(*,*) 'ERROR: add_vacant_site: reached grid limits (b): ',k,site,dxyz
            newsite = site0 + (k-1)*0.5*dxyz
            write(*,*) newsite,occupancy(newsite(1),newsite(2),newsite(3))%indx(1)
            stop
        endif
        if (occupancy(newsite(1),newsite(2),newsite(3))%indx(1) == OUTSIDE_TAG) then
            exit
        endif
    enddo
endif
if (dbug) write(*,'(a,4i6)') 'newsite: ',newsite
occupancy(newsite(1),newsite(2),newsite(3))%indx = 0
if (dbug) write(*,'(a,7i6)') 'site, vacant site: ',site,newsite

end subroutine


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine test_exp
integer :: i, N=100
real(DP) :: R, sum
integer :: kpar = 0

sum = 0
write(*,*) 'start'
do i = 1,N
	R = par_rexp(kpar)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
sum = 0
write(*,*) 'start'
do i = 1,N
!	call random_number(R)
	R = par_uni(kpar)
	R = -log(R)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
end subroutine

!-----------------------------------------------------------------------------------------
! Initial cell position data is loaded into occupancy() and cellist().
!-----------------------------------------------------------------------------------------
subroutine placeCells(ok)
logical :: ok
integer :: id, cogid, x, y, z, site(3), ctype
integer :: idc, kdc, k, x2, y2, z2, gen, stage, region
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax, nzlim, nassigned
real(DP) :: R
real :: d, d2, p1, p2, tnow, prox, tmins
integer, allocatable :: permc(:)
logical :: added, done
integer :: kpar = 0

write(logmsg,*) 'placeCells: ',NX,NY,NZ
call logger(logmsg)
ok = .false.

NBcells = 0
IDTEST = 0
tnow = 0
nlist = 0
do x = 1,NX
    do y = 1,NY
	    do z = 1,NZ
            occupancy(x,y,z)%indx = 0
            occupancy(x,y,z)%exitnum = 0
            nullify(occupancy(x,y,z)%bdry)
            site = (/x,y,z/)
			if (.not.insideEllipsoid(site)) then
				occupancy(x,y,z)%indx = OUTSIDE_TAG
			else
			    nlist = nlist+1
			endif
        enddo
    enddo
enddo

nzlim = NZ
allocate(permc(nlist))
do k = 1,nlist
    permc(k) = k
enddo
call permute(permc,nlist,kpar)
istep = 0
ntagged  = 0
id = 0
cogid = 0
done = .false.
do while (.not.done)
do x = 1,NX
    do y = 1,NY
	    do z = 1,nzlim
	        if (occupancy(x,y,z)%indx(1) == 0) then ! vacant site, not OUTSIDE_TAG or DC
                id = id+1
                lastID = id
                site = (/x,y,z/)
                gen = 1
                stage = NAIVE
                region = FOLLICLE
                ctype = 1
                if (evaluate_residence_time) then
                    ctype = 1
                elseif (calibrate_motility) then
                    if (taggable(site)) then
						ctype = TAGGED_CELL
					endif
			    else
                    ctype = select_cell_type(kpar)
                    if (ctype /= NONCOG_TYPE_TAG) then
                        ncogseed = ncogseed + 1
                    endif
                endif
                k = permc(id)
                call CreateBcell(k,cellist(k),site,ctype,gen,stage,region,ok)
                if (.not.ok) return
                occupancy(x,y,z)%indx(1) = k
                NBcells = NBcells + 1
                if (calibrate_motility .and. taggable(site)) then
                    ntagged = ntagged + 1
                    cellist(k)%ctype = TAGGED_CELL
                endif
            endif
	    enddo
	    if (done) exit
    enddo
    if (done) exit
enddo
done = .true.
enddo
deallocate(permc)
call make_cognate_list(ok)
if (.not.ok) return

write(nfout,*) 'nlist,RESIDENCE_TIME: ',nlist,RESIDENCE_TIME
nlist = id	! this is already the case for 3D blob
if (NBcells /= nlist) then
    write(logmsg,*) 'Error: inconsistent B cell and nlist counts: ', NBcells, nlist
    call logger(logmsg)
    stop
endif
Nsites = NBcells 
NBcells0 = NBcells
aRadius = (ELLIPSE_RATIO**2*Nsites*3/(4*PI))**0.33333
bRadius = aRadius/ELLIPSE_RATIO
scale_factor = real(NBC_LN)*NLN_RESPONSE/NBcells0
end subroutine

!-----------------------------------------------------------------------------
! Number of exit portals required for the current cell population, ncells. 
! For SURFACE_PORTALS, roughlyly proportional to surface area, i.e. to
! ncells^(2/3), factor = exit_fraction
! Modified to use an adjustment that was computed with:
!	Tres = 12
!	chemotaxis factor = 1
!	chemotaxis radius = 100 um
! exit_fraction = mN + c
! where
!	N = ncells
!	m = 1.607E-5
!	c = 0.00602
! 
! Better is a quadratic fit (from Excel, steady_chemo_0_1.xls):
! exit_fraction = a.x^2 + b.x + c
! where x = Ncells/1000
!-----------------------------------------------------------------------------
integer function requiredExitPortals(ncells)
integer :: ncells
real :: a, b, c, x, Fe
real, parameter :: pow = 2./3.
real, parameter :: a_nochemo_24h = -4.058E-08
real, parameter :: b_nochemo_24h = 5.590E-05
real, parameter :: c_nochemo_24h = 1.006E-02

! This is to adjust Ne when there is general exit chemotaxis, to generate steady-state
real, parameter :: K_Ke = 1.0		! 0.68 for Ke = 1.0, 1.0 for Ke = 0.0 (Pe = 0.02)

	a = a_nochemo_24h
	b = b_nochemo_24h
	c = c_nochemo_24h
	base_exit_prob = 1.0
if (FIXED_NEXITS) then
	x = NBcells0/1000
	exit_fraction = (a*x**2 + b*x + c)*24/residence_time
	requiredExitPortals = exit_fraction*NBcells0**pow + 0.5
else
	x = ncells/1000
!	Fe = (a*x**2 + b*x + c)
!	Fe = 3.028E-03*x**3.571E-01		! power law fit (7 points)
	Fe = 2.990E-03*x**3.595E-01		! power law fit (8 points, 51k - 1.1m cells)
	exit_fraction = Fe*24.0/residence_time
	requiredExitPortals = K_Ke*exit_fraction*ncells**pow + 0.5 
!	write(logmsg,'(a,f7.4)') 'exit_fraction: ',exit_fraction 
!	call logger(logmsg)
endif
!requiredExitPortals = requiredExitPortals/base_exit_prob
end function

!---------------------------------------------------------------------
! Sites within the SOI of an exit are labelled with %exitnum, unless
! the site is also within the SOI of another exit, in which case the
! site is labelled with the exit index of the nearest exit, or in
! the case of a tie the choice is made randomly.
! When a site is an exit portal location (centre), %exitnum = -(the exit index)
! Note:
! When the blob grows, more exits will be added, and when the blob
! shrinks again exits must be removed.  When an exit is removed from
! the blob we set exitlist(iexit)%ID = 0, and all sites within the SOI
! must have %exitnum either set to 0 or, if within the SOI of another
! exit, set to the ID of the nearest exit.
!---------------------------------------------------------------------
subroutine placeExits
integer :: Nex, iexit, site(3)

Lastexit = 0
Nexits = 0
if (use_traffic) then
	Nex = requiredExitPortals(NBcells0)
    max_exits = 10*Nex
    write(logmsg,*) 'NBcells0, Nex, max_exits: ',NBcells0,Nex,max_exits
    call logger(logmsg)
    allocate(exitlist(max_exits))       ! Set the array size to 10* the initial number of exits
else
    Nex = 0
    write(logmsg,*) 'No exits'
    call logger(logmsg)
    return
endif
do iexit = 1,Nex
	call addExitPortal
enddo
write(logmsg,*) 'Number of exit portals: ',Nex
call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getExitNum(iexit)
integer :: iexit
integer :: i

iexit = 0
! First look for an unused index
do i = 1,Lastexit
	if (exitlist(i)%ID == 0) then
		iexit = i
		exit
	endif
enddo
if (iexit == 0) then
	Lastexit = Lastexit + 1
	iexit = Lastexit
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine addExitPortal
integer :: site(3)
integer :: iexit

call choosePortalSite(site)
call getExitNum(iexit)
Nexits = Nexits + 1
if (Lastexit > max_exits) then
	write(logmsg,*) 'Error: addExitPortal: too many exits: need to increase max_exits: ',max_exits
	call logger(logmsg)
	stop
endif
call placeExitPortal(iexit,site)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getBoundarySite(u,site,ok)
real :: u(3)
integer :: site(3)
logical :: ok
integer :: k, x, y, z, nb
real :: r0, r

ok = .false.		
r0 = aRadius
do k = -10,10
	r = r0 + k*0.5
	x = x0 + u(1)*r
	y = y0 + u(2)*r
	z = z0 + u(3)*r
	nb = blobNeighbours(x,y,z)
	if (nb == 0) then
		exit
	elseif (nb <= 17) then
		ok = .true.		
		exit
	endif
enddo
site = (/x,y,z/)
end subroutine

!-----------------------------------------------------------------------------
! Exit sites (portals) are on a portion of the lower follicle surface, at
! T zone interface.
!-----------------------------------------------------------------------------
subroutine choosePortalSite(site)
integer :: site(3)
integer :: xex, yex, zex
real(DP) :: R
real :: prox, u(3), theta, rx
integer :: kpar = 0
logical :: ok

!call logger('choosePortalSite')
! randomly choose direction in 3D, then locate sites near this vector on the blob boundary
! Need to restrict the sinus interface by removing the follicle interface.
do
	R = par_uni(kpar)
	u(1) = 2*R - 1
	if (u(1) > 0.5) cycle
	rx = sqrt(1-u(1)*u(1))
	R = par_uni(kpar)
	theta = 2*PI*R
	u(2) = rx*cos(theta)
	u(3) = rx*sin(theta)
	call getBoundarySite(u,site,ok)
	if (.not.ok) cycle
	xex = site(1)
	yex = site(2)
	zex = site(3)
	prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
	prox = 0.5*prox
	if (tooNearExit(site,prox)) then	! too near another exit
		cycle
	endif	
	exit
enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine placeExitPortal(iexit,site)
integer :: iexit, site(3)
integer :: kexit, x, y, z, x2, y2, z2, ek(3), vk(3), ns
integer :: xex, yex, zex, xmin, xmax, ymin, ymax, zmin, zmax
real(DP) :: R
real :: d2, d2k
integer :: kpar = 0

if (iexit == Lastexit + 1) then
	Lastexit = iexit
endif
ns = 0
exitlist(iexit)%ID = iexit      ! when an exit site is lost (because the blob retracted) set %ID = 0
exitlist(iexit)%site = site
xex = site(1)
yex = site(2)
zex = site(3)
occupancy(xex,yex,zex)%exitnum = -iexit     ! This site holds an exit 
!write(logmsg,*) 'placeExitPortal: site,exitnum: ',xex,yex,zex,-iexit
!call logger(logmsg)

xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    x2 = (x-xex)*(x-xex)
    do y = ymin,ymax
        y2 = (y-yex)*(y-yex)
        do z = zmin,zmax
	        z2 = (z-zex)*(z-zex)
	        d2 = x2 + y2 + z2
	        if (d2 <= chemo_N*chemo_N) then
		        if (d2 > 0) then   ! don't touch exit site
     ! NOTE!!!  If this site is already marked as within another exit's SOI, we need to
     ! determine which exit is closer, and set %exitnum to this exit's ID
                    kexit = occupancy(x,y,z)%exitnum
                    if (kexit > 0) then
                        ek = exitlist(kexit)%site
                        vk = (/x,y,z/) - ek
                        d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
!                            write(*,*) 'other exit: ',iexit,kexit,d2,d2k
                        if (d2k < d2) then
!	                            write(*,*) 'closer exit: ',iexit,kexit,d2,d2k
                            cycle
                        endif
                        if (d2k == d2) then     ! toss a coin
                            R = par_uni(kpar)
                            if (R < 0.5) cycle
                        endif
                        occupancy(x,y,z)%exitnum = iexit  ! the site is closer to iexit than to kexit
                        ns = ns + 1
                    elseif (kexit == 0) then
                        occupancy(x,y,z)%exitnum = iexit    ! this site is closest to exit iexit
                        ns = ns + 1
                    endif
                endif
	        endif
        enddo
    enddo
enddo
!write(*,*) 'near sites: ',iexit,ns
end subroutine

!---------------------------------------------------------------------
! Remove the last nremex exits in the list.
! Perhaps it would make more sense to remove exits that are close to 
! others, or to adjust the locations of remaining exits.
!---------------------------------------------------------------------
subroutine removeExits1(nremex)
integer :: nremex
integer :: Nex, k, iexit, site(3)

!write(*,*) 'removeExits: ',nremex
Nex = Lastexit
k = 0
do iexit = Nex,1,-1
    if (exitlist(iexit)%ID == 0) cycle
    k = k+1
    site = exitlist(iexit)%site
!    write(logmsg,'(a,4i6)') 'removeExits: removeExitPortal: ',iexit,site
!    call logger(logmsg)
    call removeExitPortal(site)
!    write(logmsg,'(a,4i6)') 'did removeExitPortal'
!    call logger(logmsg)
    if (k == nremex) exit
enddo
!call checkExits("in removeExits")
end subroutine

!---------------------------------------------------------------------
! Remove exits that are close to others.
!---------------------------------------------------------------------
subroutine removeExits(nremex)
integer :: nremex
integer :: Nex, k, iexit, kexit, site1(3), site2(3)
real :: r(3), sum, d2, cmin
real, allocatable :: closeness(:)

!write(*,*) 'removeExits: ',nremex
Nex = Lastexit
allocate(closeness(Nex))
closeness = 0
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	closeness(iexit) = sum
enddo

do k = 1,nremex
	cmin = 1.0e10
	do kexit = 1,Nex
		if (closeness(kexit) == 0) cycle
		if (closeness(kexit) < cmin) then
			cmin = closeness(kexit)
			iexit = kexit
		endif
	enddo
	closeness(iexit) = 0
    site1 = exitlist(iexit)%site
!    write(logmsg,'(a,4i6)') 'removeExits: removeExitPortal: ',iexit,site
!    call logger(logmsg)
    call removeExitPortal(site1)
!    write(logmsg,'(a,4i6)') 'did removeExitPortal'
!    call logger(logmsg)
enddo
!call checkExits("in removeExits")
deallocate(closeness)
end subroutine

!---------------------------------------------------------------------
! Remove the exit portal at site.
!---------------------------------------------------------------------
subroutine removeExitPortal(site)
integer :: site(3)
integer :: iexit,xex,yex,zex,xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,ek(3),vk(3),k,kmin,i
real :: d2k, d2min

xex = site(1)
yex = site(2)
zex = site(3)
iexit = -occupancy(xex,yex,zex)%exitnum
if (iexit <= 0) then
	write(logmsg,*) 'Error: removeExitPortal: no exit at: ',site,iexit
	call logger(logmsg)
	do i = 1,Lastexit
		site = exitlist(i)%site
		write(logmsg,*) 'exit: ',i,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
	enddo
	stop
endif
exitlist(iexit)%ID = 0
occupancy(xex,yex,zex)%exitnum = iexit      ! to ensure that the site is processed in the next section
											! effectively it is treated like a site within the portal's SOI
!write(*,*) 'removeExitPortal: site,exitnum: ',site,iexit

xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    do y = ymin,ymax
        do z = zmin,zmax
            if (occupancy(x,y,z)%exitnum == iexit) then
                occupancy(x,y,z)%exitnum = 0
                ! at this point we need to check to see if (x,y,z) is within the SOI of another exit
                d2min = 9999
                kmin = 0
                do k = 1,Lastexit
                    if (exitlist(k)%ID == 0) cycle
                    ek = exitlist(k)%site
                    vk = (/x,y,z/) - ek
                    d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
		            if (d2k <= chemo_N*chemo_N) then
		                if (d2k < d2min) then
		                    d2min = d2k
		                    kmin = k
		                endif
		            endif
		        enddo
		        if (kmin > 0) then
		            occupancy(x,y,z)%exitnum = kmin
!					write(*,*) '  removeExitPortal: site,exitnum: ',x,y,z,kmin
		        endif
            endif
        enddo
    enddo
enddo
if (iexit == Lastexit) then
	Lastexit = Lastexit - 1
endif
Nexits = Nexits - 1
!write(*,*) 'removeExitPortal: Nexits: ',Nexits
end subroutine

!---------------------------------------------------------------------
! Move the exit portal at site closer to the blob centre.
!---------------------------------------------------------------------
subroutine moveExitPortalInwards(site0)
integer :: site0(3)
integer :: site(3), site1(3), iexit0, iexit1, dx, dy, dz
real :: d0, dc, d2, d2min, r(3)
real, parameter :: dlim = 1.95

iexit0 = -occupancy(site0(1),site0(2),site0(3))%exitnum
! First, remove the exit
call removeExitPortal(site0)
!write(logmsg,*) 'Removed exit within moveExitPortalInwards: ',iexit0,site0
!call logger(logmsg)
!call checkexits("in moveExitPortalInwards")
! Now add the exit back at a closer site.  Choose a new site at
! least a distance of dlim closer to blob centre.  Choose the
! site closest to site0 from among the neighbours within the
! specified distance from the centre.
! NOTE: need to check that there is not a DC or another exit portal nearby.  
! If there is then a completely new portal site will need to be chosen, 
! using choosePortalSite()
d0 = cdistance(site0)
d2min = 1.0e10
site1 = 0
do dx = -3,3
	do dy = -3,3
		do dz = -3,3
			site = site0 + (/dx,dy,dz/)
			dc = cdistance(site)
			if (d0-dc >= dlim) then
				r = site - site0
				d2 = norm2(r)
				if (d2 < d2min) then
					d2min = d2
					site1 = site
				endif
			endif
		enddo
	enddo
enddo
if (site1(1) == 0) then
	write(logmsg,*) 'Error: moveExitPortalInwards: no site to move to: ',iexit0,site0
	call logger(logmsg)
	stop
endif
call getExitNum(iexit1)
!write(logmsg,*) 'moveExitPortalInwards: ',iexit0,' from: ',site0
!call logger(logmsg)
!write(logmsg,*) '                new #: ',iexit1,'   at: ',site1
!call logger(logmsg)
Nexits = Nexits + 1
call placeExitPortal(iexit1,site1)
if (exitlist(iexit1)%site(1) /= site1(1) .or. &
	exitlist(iexit1)%site(2) /= site1(2) .or. &
	exitlist(iexit1)%site(3) /= site1(3)) then
	write(logmsg,*) 'Error: moveExitPortalInwards: bad site: ',iexit1,exitlist(iexit1)%site,site1
	call logger(logmsg)
	stop
endif
!write(logmsg,*) 'now: ',iexit1,' is at: ',exitlist(iexit1)%site
!call logger(logmsg)
!write(logmsg,*) 'and exit #10 is at: ',exitlist(10)%site
!call logger(logmsg)
!write(logmsg,'(a,6i6)') 'moveExitPortalInwards: placeExitPortal: ',istep,iexit,site1,Nexits
!call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
! The exit portal iexit needs to be moved by one site in the direction 
! closest to that given by the unit vector v(:).
!---------------------------------------------------------------------
subroutine moveExitPortal(iexit0,v)
integer :: iexit0
real :: v(3)
integer :: iexit1, site0(3), site1(3), k, kmax
real :: proj, pmax, jump(3)

pmax = 0
do k = 1,27
	if (k == 14) cycle
	jump = jumpvec(:,k)
	proj = dot_product(jump,v)/norm(jump)
	if (proj > pmax) then
		pmax = proj
		kmax = k
	endif
enddo
site0 = exitlist(iexit0)%site
site1 = site0 + jumpvec(:,kmax)
write(logmsg,'(a,i4,a,3i4,a,3i4)') 'moveExitPortal: ',iexit0,' from: ',site0,' to: ',site1
call logger(logmsg)
call removeExitPortal(site0)
call getExitNum(iexit1)
Nexits = Nexits + 1
call placeExitPortal(iexit1,site1)
end subroutine

!---------------------------------------------------------------------
! Display exit portal location info
!---------------------------------------------------------------------
subroutine displayExits
integer :: Nex, k, iexit, kexit, site1(3), site2(3), count
real :: r(3), sum, d2

!write(*,*) 'removeExits: ',nremex
Nex = Lastexit
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	count = neighbourhoodCount(site1)
	write(logmsg,'(5i4,3f10.3)') iexit,site1,count,cdistance(site1),(PI/180)*atan2(site1(2)-Centre(2),site1(3)-Centre(3)),1000*sum
	call logger(logmsg)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! need to ensure that exits do not get too close together
!-----------------------------------------------------------------------------------------
subroutine checkExitSpacing(msg)
character*(*) :: msg
integer :: iexit1, iexit2, site1(3), site2(3), im1=0, im2=0
real :: r(3), d, dmin, f(3), df(3)
logical :: ok = .true.
real :: dlim

dlim = exit_prox*chemo_radius
dmin = 1.0e10
do iexit1 = 1,Lastexit
	if (exitlist(iexit1)%ID == 0) cycle
	site1 = exitlist(iexit1)%site
	f = 0
	do iexit2 = 1,Lastexit
		if (iexit1 == iexit2) cycle
		if (exitlist(iexit2)%ID == 0) cycle
		site2 = exitlist(iexit2)%site
		r = site1 - site2
		d = norm(r)
		if (d < dmin) then
			dmin = d
			im1 = iexit1
			im2 = iexit2
		endif
		if (d < dlim) then
			write(logmsg,'(2i4,f6.2)') iexit1,iexit2,d
			call logger(logmsg)
			df = r/d		! unit vector in direction of force
			df = df/d		! scale by inverse distance
			f = f + df
		endif
	enddo
	if (f(1) /= 0 .or. f(2) /= 0 .or. f(3) /= 0) then
		f = f/norm(f)
		call moveExitPortal(iexit1,f)
	endif
enddo
write(logmsg,*) msg,': Min exit spacing: ',im1,im2,dmin
call logger(logmsg)
end subroutine



!--------------------------------------------------------------------------------------
! Motility behaviour
! There are potentially 4 motility states:
! (1) Naive cell, very short non-cognate DC contacts 3 min (motility level 1)
! (2) Short cognate DC contacts 11 min (motility level 2)
! (3) Long cognate DC contacts > 1 hr clusters (motility level 3)
! (4) Short cognate DC contacts 18 min  swarms (motility level 4)
! It isn't clear what level of motility occurs at each stage. Possibly
! the contact duration is related to the motility, since reduced motility
! leads to a higher probability of rebinding.
! A simple way to vary motility is to keep the persistence parameter rho fixed
! and vary alpha (plongjump), the probability of 2 steps, and possibly beta,
! the probability of moving at all.  In any case the jump parameters for each case
! are stored as pjump(2) and dirprob(0:6).
! In principle: speed & Cm -> rho, beta, delta -> dirprob(), pjump().
!
! Mark Miller says that their motility measurements didn't exclude periods when
! T cells were in contact with DC, therefore we have no info about motility in
! various stages of activation.  For now it is safest to use the same motility
! parameters for all stages.
!--------------------------------------------------------------------------------------

!---------------------------------------------------------------------
! Updates T cell state, and DC %stimulation (for use in modifying
! %density in update_DCstate).
!---------------------------------------------------------------------
subroutine updater(ok)
logical :: ok
integer :: kcell, ctype, stype, region, iseq, tag, kfrom, kto, k, ncog, ntot
integer :: site(3), site2(3), freeslot, indx(2), status, DC(2), idc
!real :: C(N_CYT), mrate(N_CYT)
real :: tnow, dstim, S, cyt_conc, mols_pM, Ctemp, dstimrate, stimrate
logical :: divide_flag, producing, first, dbg, unbound, flag, flag1
!logical, save :: first = .true.
type (cog_type), pointer :: p

!write(*,*) 'updater: ',me
ok = .true.
dbg = .false.
flag = .false.
flag1 = .false.
ntot = 0
ncog = 0
tnow = istep*DELTA_T
! Scaling factor to convert total number of molecules in the region to conc in pM
!mols_pM = L_um3*M_pM/(NBcells*Vc*Navo)

!do kcell = 1,nlist
!    if (cellist(kcell)%ID == 0) cycle
!    ntot = ntot + 1
!    if (dbg) write(*,*) 'kcell: ',kcell
!    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
!    if (stype /= COG_TYPE_TAG) cycle
    ! Only cognate cells considered
do k = 1,lastcogID
    kcell = cognate_list(k)
    ncog = ncog + 1
    p => cellist(kcell)%cptr
    if (.not.associated(p)) then
        write(logmsg,*) 'ERROR: updater: p not associated: ',kcell
	    call logger(logmsg)
        stop
    endif
    call get_region(p,region)
	! Cell death
    if (tnow > p%dietime) then
        call Bcell_death(kcell)
        cycle
    endif

! Stage transition
    call updatestage(kcell, tnow, divide_flag)

! Cell division
    if (divide_flag) then
		indx = occupancy(site(1),site(2),site(3))%indx
		freeslot = 0
		if (indx(1) == 0) then
			site2 = site
			freeslot = 1
		elseif (indx(2) == 0) then
			site2 = site
			freeslot = 2
		else
			call get_free_slot(occupancy,NX,site,site2,freeslot)
		endif
        if (freeslot /= 0) then     ! there is a free slot at site2 (which may be = site)
            call cell_division(kcell,site2,freeslot,ok)
            if (.not.ok) return
        endif
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------------
! Note: USE_STAGETIME(icstage) means use stagetime for the transition from icstage.
! p%stagetime = time that the cell is expected to make transition to the next stage.
!--------------------------------------------------------------------------------------
subroutine updatestage(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, region, gen, ctype, site(3)
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
site = cellist(kcell)%site
p => cellist(kcell)%cptr
!stage = get_stage(p)
call get_stage(p,stage,region)
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)
		if (AntigenEncounter(site)) then
	        call set_stage(p,ANTIGEN_MET)
	        p%stagetime = tnow + T_CCR7_UP		
		endif		
	case (ANTIGEN_MET)    ! possible transition from ANTIGEN_MET to CCR7_UP
        call set_stage(p,CCR7_UP)
		p%stagetime = tnow
	case (CCR7_UP)     ! possible transition from CCR7_UP to TCELL_MET
		if (TCellEncounter(site)) then
	        call set_stage(p,TCELL_MET)
	        p%stagetime = tnow + T_CCR7_DOWN		
		endif
	case (TCELL_MET)
		call set_stage(p,DIVIDING)
        p%stagetime = tnow + max(0.,T_FIRST_DIVISION - T_CCR7_DOWN)
	case (DIVIDING)	
        gen = get_generation(p)
		if (gen == BC_MAX_GEN) then
            call set_stage(p,FINISHED)
			p%stagetime = BIG_TIME
		else
		    divide_flag = .true.
		endif
	case (BCL6_UP)
		call set_stage(p,DIVIDING)
        p%stagetime = tnow + max(0.,T_DIVISION - T_BCL6_UP)		
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------------
! Check for (probabilistic) encounter of a B cell with antigen, near the upper surface.
!--------------------------------------------------------------------------------------
logical function AntigenEncounter(site)
integer :: site(3)
integer :: v(3), x, y, z
real :: r2

AntigenEncounter = .false.
v = site - Centre
x = v(1)
y = v(2)
z = v(3)
if (y < 0) return
r2 = (x/aRadius)**2 + (y**2 + z**2)/bRadius**2
if (r2 < 0.9) then
	AntigenEncounter = .true.
	write(*,*) 'AntigenEncounter: ',site
endif
end function

!--------------------------------------------------------------------------------------
! Check for (probabilistic) encounter of a B cell with a T helper cell, near the lower surface.
!--------------------------------------------------------------------------------------
logical function TCellEncounter(site)
integer :: site(3)
integer :: v(3), x, y, z
real :: r2

TCellEncounter = .false.
v = site - Centre
x = v(1)
y = v(2)
z = v(3)
if (y > 0) return
r2 = (x/aRadius)**2 + (y**2 + z**2)/bRadius**2
if (r2 < 0.9) then
	TCellEncounter = .true.
	write(*,*) 'TCellEncounter: ',site
endif
end function


end module
