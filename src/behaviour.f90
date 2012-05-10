! B cell behaviour
module behaviour

use global
use fields
use FDC
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
real :: T_half_lo, T_half_hi
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

! B cell death rates, from half-lifes
T_half_lo = 10*24*60.   ! 10 days in minutes
T_half_hi = 2*24*60.    ! 1 day in minutes
Kdeath_lo(:) = log(2.)/T_half_lo    ! rate constant
Kdeath_hi(:) = log(2.)/T_half_hi    ! rate constant
end subroutine

!--------------------------------------------------------------------------------------
! If a cognate cell is NAIVE it is long-lived - effectively it does not die.
! When a cognate cell has contact with a cognate DC its lifetime is immediately limited.
! NOTE: This is not used.  For Bcells we are using death rates (see DeathProbability())
!--------------------------------------------------------------------------------------
real function BClifetime(ptr)
type (cog_type), pointer :: ptr
integer :: gen, stage
real :: p1, p2
integer :: kpar = 0

BClifetime = BIG_TIME
return

!call get_stage(ptr,stage,region)
stage = get_stage(ptr)
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
real function DivisionTime()
integer :: kpar = 0
real, parameter :: rndfraction = 0.2
real(DP) :: R

R = par_uni(kpar)
DivisionTime = T_DIVISION*(1 + rndfraction*(2*R-1))
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
integer :: i, ic, shownoncog, ncpu_dummy, iuse(MAX_RECEPTOR), iuse_rate(MAX_CHEMO)
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
read(nfcell,*) iuse(S1PR1)
read(nfcell,*) iuse(S1PR2)
read(nfcell,*) iuse_rate(S1P)
read(nfcell,*) chemo(S1P)%bdry_rate
read(nfcell,*) chemo(S1P)%bdry_conc
read(nfcell,*) chemo(S1P)%diff_coef
read(nfcell,*) chemo(S1P)%halflife
read(nfcell,*) receptor(S1PR1)%strength
read(nfcell,*) receptor(S1PR2)%strength
read(nfcell,*) iuse(CCR7)
read(nfcell,*) iuse_rate(CCL21)
read(nfcell,*) chemo(CCL21)%bdry_rate
read(nfcell,*) chemo(CCL21)%bdry_conc
read(nfcell,*) chemo(CCL21)%diff_coef
read(nfcell,*) chemo(CCL21)%halflife
read(nfcell,*) receptor(CCR7)%strength
read(nfcell,*) iuse(EBI2)
read(nfcell,*) iuse_rate(OXY)
read(nfcell,*) chemo(OXY)%bdry_rate
read(nfcell,*) chemo(OXY)%bdry_conc
read(nfcell,*) chemo(OXY)%diff_coef
read(nfcell,*) chemo(OXY)%halflife
read(nfcell,*) receptor(EBI2)%strength
read(nfcell,*) iuse(CXCR5)
read(nfcell,*) iuse_rate(CXCL13)
read(nfcell,*) chemo(CXCL13)%bdry_rate
read(nfcell,*) chemo(CXCL13)%bdry_conc
read(nfcell,*) chemo(CXCL13)%diff_coef
read(nfcell,*) chemo(CXCL13)%halflife
read(nfcell,*) receptor(CXCR5)%strength
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
chemo(S1P)%name    = 'S1P'
chemo(CCL21)%name  = 'CCL21'
chemo(OXY)%name    = 'OXY'
chemo(CXCL13)%name = 'CXCL13'
receptor(S1PR1)%name = 'S1PR1'
receptor(CCR7)%name = 'CCR7'
receptor(EBI2)%name = 'EBI2'
receptor(CXCR5)%name = 'CXCR5'
receptor(S1PR2)%name = 'S1PR2'
receptor(S1PR1)%chemokine = S1P
receptor(CCR7)%chemokine = CCL21
receptor(EBI2)%chemokine = OXY
receptor(CXCR5)%chemokine = CXCL13
receptor(S1PR2)%chemokine = S1P
chemo(:)%used = .false.
do i = 1,MAX_RECEPTOR
	receptor(i)%used = (iuse(i) == 1)		! interim measure
	if (receptor(i)%used .or. chemo(receptor(i)%chemokine)%used) then
		chemo(receptor(i)%chemokine)%used = .true.
	endif
	receptor(i)%level = receptor_level(i,:)
	if (i == S1PR2) then
		receptor(i)%sign = -1
	else
		receptor(i)%sign = 1
	endif
enddo
do ic = 1,MAX_CHEMO
	chemo(ic)%use_secretion = (iuse_rate(ic) == 1)
	if (chemo(ic)%use_secretion .and. .not.use_ODE_diffusion) then
		write(logmsg,*) 'Error: simulating chemokine ',chemo(ic)%name,' secretion requires use_ODE_diffusion'
		call logger(logmsg)
		ok = .false.
		return
	endif
	chemo(ic)%decay_rate = DecayRate(chemo(ic)%halflife)
enddo

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
write(nfout,*) iuse
call ShowChemokines

chemo_radius = chemo_radius/DELTA_X				! convert from um to lattice grids
chemo_N = max(3,int(chemo_radius + 0.5))	    ! convert from um to lattice grids
chemo_exp = log(1/CHEMO_MIN)/log(chemo_radius)

sigma = log(divide_shape1)
divide_dist1%p1 = log(60*divide_mean1/exp(sigma*sigma/2))
divide_dist1%p2 = sigma
sigma = log(divide_shape2)
divide_dist2%p1 = log(60*divide_mean2/exp(sigma*sigma/2))
divide_dist2%p2 = sigma

!sigma = log(DC_ANTIGEN_SHAPE)
!DC_ANTIGEN_MEDIAN = DC_ANTIGEN_MEAN/exp(sigma*sigma/2)
!sigma = log(DC_LIFETIME_SHAPE)
!DC_LIFETIME_MEDIAN = DC_LIFETIME_MEAN/exp(sigma*sigma/2)


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
write(nfout,*) 'use_cognate: ',use_cognate
write(nfout,*) 'random_cognate: ',random_cognate
write(nfout,*) 'use_ode_diffusion: ',use_ode_diffusion
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

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowChemokines
integer :: i
character*(8) :: usage
type(receptor_type), pointer :: p=>null()
type(chemokine_type), pointer :: q=>null()

write(nfout,'(a)') 'Receptors:'
do i = 1,MAX_RECEPTOR
	p => receptor(i)
	if (p%used) then
		usage = 'USED'
	else
		usage = 'NOT USED'
	endif
	write(nfout,'(a,i2,2x,a,a,a,4x,a)') 'Receptor: ',i,p%name,'  chemokine: ',chemo(p%chemokine)%name,usage
	write(nfout,'(a,f8.2)') 'relative strength: ',p%sign*p%strength
	write(nfout,'(a,4f6.1)') 'stage levels: ',p%level
enddo
write(nfout,'(a)') 'Chemokines:'
do i = 1,MAX_CHEMO
	q => chemo(i)
	if (q%used) then
		usage = 'USED'
	else
		usage = 'NOT USED'
	endif
	write(nfout,'(a,i2,2x,a,4x,a)') 'Chemokine: ',i,q%name,usage
	write(nfout,'(a,f8.1)') 'Boundary conc: ',q%bdry_conc
	write(nfout,'(a,f8.1)') 'Diffusion coeff: ',q%diff_coef
	write(nfout,'(a,f8.1)') 'Halflife: ',q%halflife
	write(nfout,'(a,f8.4)') 'Decay rate: ',q%decay_rate
enddo
	
end subroutine

!-----------------------------------------------------------------------------------------
! When a T cell dies it is removed from the cell list (%ID -> 0)
! and removed from occupancy()%indx()
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine BcellDeath(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2), ctype, stype, region
logical :: cognate

!write(logmsg,*) 'BcellDeath: ',kcell
!call logger(logmsg)
cognate = (associated(cellist(kcell)%cptr))
if (cognate) then
!	call get_region(cellist(kcell)%cptr,region)
	region = get_region(cellist(kcell)%cptr)
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
    write(logmsg,*) 'BcellDeath: ngaps > max_ngaps'
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
    write(logmsg,*) 'ERROR: BcellDeath: cell not at site: ',kcell,site,indx
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
    write(logmsg,*) 'BcellDeath: No site removal code, use add count'
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
integer :: icnew, ctype, gen, region, status, site(3), indx(2)
integer :: iseq, tag, kfrom, kto
real :: tnow, prob_gcc, prob_plasma
type(cog_type), pointer :: p1, p2
!real :: IL_state(CYT_NP)
integer :: kpar = 0
real(DP) :: R

!write(*,*) 'cell_division: ',kcell
ok = .true.
tnow = istep*DELTA_T
!call show_cognate_cell(kcell)
p1 => cellist(kcell)%cptr
gen = get_generation(p1)
!call get_region(p1,region)
region = get_region(p1)
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
status = get_status(p1)
if (status == BCL6_LO) then
	prob_gcc = get_GccProb(gen)
	prob_plasma = get_PlasmaProb(gen) 
	R = par_uni(kpar)	
	if (R > prob_gcc + prob_plasma) then
		call set_stage(p1,DIVIDING)
		p1%stagetime = tnow + DivisionTime()
	elseif (R < prob_gcc) then
		call set_status(p1,BCL6_HI)
		call set_stage(p1,GCC_COMMIT)
		p1%stagetime = tnow
	else
		call set_status(p1,PLASMA)
		call set_stage(p1,PLASMA)
		p1%stagetime = tnow
	endif
else
	call set_stage(p1,DIVIDING)
	p1%stagetime = tnow + DivisionTime()
endif
!call set_stage(p1,POST_DIVISION)
!if (TCR_splitting) then
!    p1%stimulation = p1%stimulation/2
!endif
ctype = cellist(kcell)%ctype
!p1%dietime = tnow + BClifetime(p1)
!p1%dividetime = tnow
!p1%stagetime = tnow + dividetime(gen,ctype)
site = cellist(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx

call CreateBcell(icnew,cellist(icnew),site2,ctype,gen,DIVIDING,region,ok)
if (.not.ok) return

!write(*,*) 'New cell: ',icnew

p2 => cellist(icnew)%cptr

p2%stimulation = p1%stimulation
p2%generation = p1%generation
p2%region = p1%region
cellist(icnew)%entrytime = tnow
if (status == BCL6_LO) then
	prob_gcc = get_GccProb(gen)
	prob_plasma = get_PlasmaProb(gen) 
	R = par_uni(kpar)	
	if (R > prob_gcc + prob_plasma) then
		call set_stage(p2,DIVIDING)
		p2%stagetime = tnow + DivisionTime()
	elseif (R < prob_gcc) then
		call set_status(p2,BCL6_HI)
		call set_stage(p2,GCC_COMMIT)
		p2%stagetime = tnow
	else
		call set_status(p2,PLASMA)
		call set_stage(p2,PLASMA)
		p2%stagetime = tnow
	endif
else
	call set_status(p2,status)
	call set_stage(p2,DIVIDING)
	p2%stagetime = tnow + DivisionTime()
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
!    cell%cptr%status = 0
    call set_generation(cell%cptr,gen)
!    call set_stage_region(cell%cptr,stage,region)
	call set_stage(cell%cptr,stage)
	call set_region(cell%cptr,region)
    if (fix_avidity) then
        i = mod(navid,avidity_nlevels)
        navid = navid + 1
        cell%cptr%avidity = avidity_level(i+1)
    else
        cell%cptr%avidity = rv_lognormal(param1,param2,kpar)
    endif
    cell%cptr%stimulation = 0
    cell%cptr%status = BCL6_LO	! default
!    if (use_cytokines) then
!        call IL2_init_state(cell%cptr%IL_state,cell%cptr%IL_statep)
!    endif
    ! What should the initial CD69 level be?  If 0, can a cognate cell exit immediately?
    ! We would prefer not to impose a time or generation constraint on the exit of
    ! cognate T cells, but otherwise if CD69 is initially 0, a cell will be susceptible
    ! to chemotaxis and exit until it has received enough TCR signal to drive CD69 high.
!    cell%cptr%CD69 = 0
!    cell%cptr%S1P1 = 0
!    cell%cptr%dietime = tnow + BClifetime(cell%cptr)
    cell%cptr%dietime = BIG_TIME
    cell%cptr%dividetime = BIG_TIME
    cell%cptr%stagetime = 0

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
cell%receptor_level = receptor%level(NAIVE_TAG)
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
        write(logmsg,*) 'ERROR: add_vacant_site: reached grid limits (a): ',k,site,dxyz
        call logger(logmsg)
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
            write(logmsg,*) 'ERROR: add_vacant_site: reached grid limits (b): ',k,site,dxyz
            call logger(logmsg)
            newsite = site0 + (k-1)*0.5*dxyz
            write(logmsg,*) newsite,occupancy(newsite(1),newsite(2),newsite(3))%indx(1)
            call logger(logmsg)
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
real function get_GccProb(gen)
integer :: gen

if (gen <= 2) then
	get_GccProb = 0
elseif (gen == 3) then
	get_GccProb = 0.1
else
	get_GccProb = 0.3
endif
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function get_PlasmaProb(gen)
integer :: gen

if (gen <= 3) then
	get_PlasmaProb = 0
else
	get_PlasmaProb = PLASMA_PROB
endif
end function

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
subroutine PlaceCells(ok)
logical :: ok
integer :: id, cogid, x, y, z, site(3), ctype
integer :: idc, kdc, k, x2, y2, z2, gen, stage, region
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax, nzlim, nassigned
real(DP) :: R
real :: d, d2, p1, p2, tnow, prox, tmins
integer, allocatable :: permc(:)
logical :: added, done
integer :: kpar = 0

write(logmsg,*) 'PlaceCells: ',NX,NY,NZ
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
            occupancy(x,y,z)%FDC_nbdry = 0
!            occupancy(x,y,z)%exitnum = 0
            nullify(occupancy(x,y,z)%bdry)
            site = (/x,y,z/)
			if (.not.InsideEllipsoid(site)) then
				occupancy(x,y,z)%indx = OUTSIDE_TAG
			else
			    nlist = nlist+1
			endif
        enddo
    enddo
enddo

call placeFDCs(ok)
if (.not.ok) stop

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
                if (id == nlist) then
					done = .true.
					exit
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
!	prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
	prox = 0.5*prox
	if (tooNearExit(site,prox)) then	! too near another exit
		cycle
	endif	
	exit
enddo

end subroutine



!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine ReceptorLevel(kcell,ifrom,ito,t,t1,t2,level)
integer :: kcell, ifrom, ito
real :: t, t1, t2, level(:)
real :: alpha, rfrom, rto
integer :: ireceptor

if (t < t1) then
	write(logmsg,'(a,i6,3f8.2)') 'Error: ReceptorLevel: t<t1: ',kcell,t,t1,t2
	call logger(logmsg)
	stop
endif
if (t > t2) then
	alpha = 1
else
	alpha = (t-t1)/(t2-t1)
endif
do ireceptor = 1,MAX_CHEMO
	if (receptor(ireceptor)%used) then
		rfrom = receptor(ireceptor)%level(ifrom)
		rto = receptor(ireceptor)%level(ito)
		level(ireceptor) = (1-alpha)*rfrom + alpha*rto
	endif
enddo
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
! For cognate cells, updates T cell state, and DC %stimulation 
! (for use in modifying %density in update_DCstate).
!
! Treatment of cell death
!------------------------
! There are two possibilities:
! (1) Specify a time-to-die for each cell (i.e. on cell division)
! (2) Use a probability of death in every time step (i.e. a rate of death)
! Previously (1) was implemented.  Now implement (2) for activated
! cognate cells, BCL6lo and BCL6hi, with generation-dependent rates.
!---------------------------------------------------------------------
subroutine updater(ok)
logical :: ok
integer :: kcell, stage, iseq, tag, kfrom, kto, k, ncog, ntot
integer :: site(3), site2(3), freeslot, indx(2), status, DC(2), idc
integer :: tmplastcogID
!real :: C(N_CYT), mrate(N_CYT)
real :: tnow, dstim, S, cyt_conc, mols_pM, Ctemp, dstimrate, stimrate
real :: t1, t2, die_prob
logical :: divide_flag, producing, first, dbg, unbound, flag, flag1
!logical, save :: first = .true.
real(DP) :: R
integer :: kpar = 0
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
tmplastcogID = lastcogID
do k = 1,tmplastcogID
    kcell = cognate_list(k)
    if (kcell == 0) then
!        write(*,*) 'updater: kcell = 0: ',k,tmplastcogID
        cycle
    endif
    ncog = ncog + 1
    p => cellist(kcell)%cptr
    if (.not.associated(p)) then
        write(logmsg,*) 'ERROR: updater: p not associated: ',kcell
	    call logger(logmsg)
        stop
    endif

	! Cell death
!    if (tnow > p%dietime) then
!        call BcellDeath(kcell)
!        cycle
!    endif
	die_prob = DeathProbability(p)
	R = par_uni(kpar)
	if (R < die_prob) then
        call BcellDeath(kcell)
        cycle
    endif

	stage = get_stage(p)
	if (stage == ANTIGEN_MET) then
		t2 = p%stagetime
		t1 = t2 - T_CCR7_UP
!		write(*,*) 'updater: tnow,t1,t2: ',kcell,tnow,t1,t2
		call ReceptorLevel(kcell,NAIVE_TAG,ANTIGEN_TAG,tnow,t1,t2,cellist(kcell)%receptor_level)		
	elseif (stage == TCELL_MET) then
		t2 = p%stagetime
		t1 = t2 - T_EBI2_UP
		call ReceptorLevel(kcell,ANTIGEN_TAG,ACTIVATED_TAG,tnow,t1,t2,cellist(kcell)%receptor_level)		
	elseif (stage == GCC_COMMIT) then
		t2 = p%stagetime
		t1 = t2 - T_BCL6_UP
		call ReceptorLevel(kcell,ACTIVATED_TAG,GCC_TAG,tnow,t1,t2,cellist(kcell)%receptor_level)		
	endif
	
! Stage transition
    call updatestage(kcell, tnow, divide_flag)

! Cell division
    if (divide_flag) then
		site = cellist(kcell)%site
		indx = occupancy(site(1),site(2),site(3))%indx
		freeslot = 0
		if (indx(1) == 0) then
			site2 = site
			freeslot = 1
		elseif (indx(2) == 0) then
			site2 = site
			freeslot = 2
		else
			call get_free_slot(NX,site,site2,freeslot)
		endif
        if (freeslot /= 0) then     ! there is a free slot at site2 (which may be = site)
            call cell_division(kcell,site2,freeslot,ok)
            if (.not.ok) return
        endif
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------------
! stage is the current stage of the cell, 
! p%stagetime = time that the cell is expected to make transition to the next stage.
!--------------------------------------------------------------------------------------
subroutine updatestage(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, gen, ctype, site(3)
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
site = cellist(kcell)%site
p => cellist(kcell)%cptr
!call get_stage(p,stage,region)
stage = get_stage(p)
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)
		if (AntigenEncounter(kcell,site)) then
	        call set_stage(p,ANTIGEN_MET)
	        p%stagetime = tnow + T_CCR7_UP	
	        write(logmsg,'(a,i6,2f8.2)') 'updatestage: -> ANTIGEN_MET: ',kcell,tnow,p%stagetime	
	        call logger(logmsg)
		endif		
	case (ANTIGEN_MET)    ! possible transition from ANTIGEN_MET to CCR7_UP
        call set_stage(p,CCR7_UP)
		p%stagetime = tnow
	case (CCR7_UP)     ! possible transition from CCR7_UP to TCELL_MET
		if (TCellEncounter(kcell,site)) then
	        call set_stage(p,TCELL_MET)
	        p%stagetime = tnow + T_EBI2_UP		
		endif
	case (TCELL_MET)
		call set_stage(p,DIVIDING)
        p%stagetime = tnow + max(0.,T_FIRST_DIVISION - T_EBI2_UP)
	case (DIVIDING)		! the next stage for DIVIDING cells is set in cell_division(), call triggered by divide_flag
        gen = get_generation(p)
		if (gen == BC_MAX_GEN) then
            call set_stage(p,FINISHED)
			p%stagetime = BIG_TIME
		else
		    divide_flag = .true.
		endif
	case (GCC_COMMIT)
        call set_stage(p,BCL6_UP)
		p%stagetime = tnow + T_BCL6_UP
	case (BCL6_UP)
		call set_stage(p,DIVIDING)
        p%stagetime = tnow + max(0.,DivisionTime() - T_BCL6_UP)		
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------------
! Check for (probabilistic) encounter of a B cell with antigen, near the upper surface.
!--------------------------------------------------------------------------------------
logical function AntigenEncounter(kcell,site)
integer :: kpar = 0
integer :: kcell,site(3)
integer :: v(3), x, y, z
real :: r2
real :: encounter_prob = 0.1

AntigenEncounter = .false.
v = site - Centre
x = v(1)
y = v(2)
z = v(3)
if (y < 0) return
r2 = (x/aRadius)**2 + (y**2 + z**2)/bRadius**2
if (r2 > 0.9 .and. par_uni(kpar) < encounter_prob) then
	AntigenEncounter = .true.
	write(logmsg,*) 'AntigenEncounter: ',kcell,site
	call logger(logmsg)
endif
end function

!--------------------------------------------------------------------------------------
! Check for (probabilistic) encounter of a B cell with a T helper cell, near the lower surface.
!--------------------------------------------------------------------------------------
logical function TCellEncounter(kcell,site)
integer :: kcell,site(3)
integer :: kpar = 0
integer :: v(3), x, y, z
real :: r2
real :: encounter_prob = 0.1

TCellEncounter = .false.
v = site - Centre
x = v(1)
y = v(2)
z = v(3)
if (y > 0) return
r2 = (x/aRadius)**2 + (y**2 + z**2)/bRadius**2
if (r2 > 0.9 .and. par_uni(kpar) < encounter_prob) then
	TCellEncounter = .true.
	write(logmsg,*) 'TCellEncounter: ',kcell,site
	call logger(logmsg)
endif
end function

!--------------------------------------------------------------------------------------
! The death probability in a time step depends on status (BCL6_LO, BCL6_HI, PLASMA)
! and generation.
!--------------------------------------------------------------------------------------
real function DeathProbability(p)
type(cog_type), pointer :: p
integer :: status, gen

status = get_status(p)
if (status == PLASMA) then
	DeathProbability = 0
elseif (status == BCL6_LO) then
	gen = min(get_generation(p),MMAX_GEN)
	DeathProbability = Kdeath_lo(gen)*DELTA_T
elseif (status == BCL6_HI) then
	gen = min(get_generation(p),MMAX_GEN)
	DeathProbability = Kdeath_hi(gen)*DELTA_T
endif
end function

end module
