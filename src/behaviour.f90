! B cell behaviour
module behaviour

use global
use fields
use FDC
use motility

implicit none

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

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
real :: T_HALF_LO, T_HALF_HI
integer :: i

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
T_HALF_LO = 10*24*60.   ! 10 days in minutes
T_HALF_HI = 2*24*60.    ! 1 day in minutes
Kdeath_lo(:) = log(2.)/T_HALF_LO    ! rate constant
Kdeath_hi(:) = log(2.)/T_HALF_HI    ! rate constant
end subroutine

!--------------------------------------------------------------------------------------
! If a cognate cell is NAIVE it is long-lived - effectively it does not die.
! When a cognate cell has contact with a cognate DC its lifetime is immediately limited.
! NOTE: This is not used.  For Bcells we are using death rates (see DeathProbability())
!--------------------------------------------------------------------------------------
real function BClifetime(p)
type(cog_type), pointer :: p
integer :: gen, stage
real :: p1, p2
integer :: kpar = 0

BClifetime = BIG_TIME
return

stage = get_stage(p)
if (stage == NAIVE) then
    BClifetime = BIG_TIME
    return
endif
gen = get_generation(p)
if (gen < 1 .or. gen > BC_MAX_GEN) then
    write(logmsg,*) 'BClifetime: bad gen: ',gen
    call logger(logmsg)
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
select case (divide_dist(gen)%class)
case (NORMAL_DIST)
	dividetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	dividetime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	dividetime = p1
end select
end function

!--------------------------------------------------------------------------------------
! idist = 1 is the first division - for now just use the median (put p2 = 0)
!       = 2 is later divisions, log-normally distributed
!--------------------------------------------------------------------------------------
real function DivisionTime(idist)
integer :: idist
integer :: kpar = 0
real, parameter :: rndfraction = 0.2
real(DP) :: R

!R = par_uni(kpar)
!DivisionTime = T_DIVISION*(1 + rndfraction*(2*R-1))
if (idist == 1) then
	DivisionTime = rv_lognormal(divide_dist2%p1,0.0,kpar)
else
	DivisionTime = rv_lognormal(divide_dist2%p1,divide_dist2%p2,kpar)
endif
!write(logmsg,*) 'DivisionTime: ',idist,DivisionTime
!call logger(logmsg)
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

!----------------------------------------------------------------------------------------1123
!----------------------------------------------------------------------------------------
subroutine read_Bcell_params(ok)
logical :: ok
real :: sigma, divide_mean1, divide_shape1, divide_mean2, divide_shape2
integer :: i, ic, k, shownoncog, ncpu_dummy, iuse(MAX_RECEPTOR), iuse_rate(MAX_CHEMO)
integer :: usetraffic, usechemo, computedoutflow, itestcase
character(4) :: logstr

ok = .true.
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
do k = 1,5
	read(nfcell,*) receptor_level(S1PR1,k)
enddo
read(nfcell,*) receptor(S1PR1)%saturation_threshold
read(nfcell,*) receptor(S1PR1)%refractory_time
do k = 1,5
	read(nfcell,*) receptor_level(S1PR2,k)
enddo
read(nfcell,*) receptor(S1PR2)%saturation_threshold
read(nfcell,*) receptor(S1PR2)%refractory_time

read(nfcell,*) iuse(CCR7)
read(nfcell,*) iuse_rate(CCL21)
read(nfcell,*) chemo(CCL21)%bdry_rate
read(nfcell,*) chemo(CCL21)%bdry_conc
read(nfcell,*) chemo(CCL21)%diff_coef
read(nfcell,*) chemo(CCL21)%halflife
read(nfcell,*) receptor(CCR7)%strength
do k = 1,5
	read(nfcell,*) receptor_level(CCR7,k)
enddo
read(nfcell,*) receptor(CCR7)%saturation_threshold
read(nfcell,*) receptor(CCR7)%refractory_time

read(nfcell,*) iuse(EBI2)
read(nfcell,*) iuse_rate(OXY)
read(nfcell,*) chemo(OXY)%bdry_rate
read(nfcell,*) chemo(OXY)%bdry_conc
read(nfcell,*) chemo(OXY)%diff_coef
read(nfcell,*) chemo(OXY)%halflife
read(nfcell,*) receptor(EBI2)%strength
do k = 1,5
	read(nfcell,*) receptor_level(EBI2,k)
enddo
read(nfcell,*) receptor(EBI2)%saturation_threshold
read(nfcell,*) receptor(EBI2)%refractory_time

read(nfcell,*) iuse(CXCR5)
read(nfcell,*) iuse_rate(CXCL13)
read(nfcell,*) chemo(CXCL13)%bdry_rate
read(nfcell,*) chemo(CXCL13)%bdry_conc
read(nfcell,*) chemo(CXCL13)%diff_coef
read(nfcell,*) chemo(CXCL13)%halflife
read(nfcell,*) receptor(CXCR5)%strength
do k = 1,5
	read(nfcell,*) receptor_level(CXCR5,k)
enddo
read(nfcell,*) receptor(CXCR5)%saturation_threshold
read(nfcell,*) receptor(CXCR5)%refractory_time

read(nfcell,*) BASE_NFDC			        ! base number of FDCs
read(nfcell,*) BASE_NMRC			        ! base number of MRCs
read(nfcell,*) base_exit_prob               ! base probability of exit at a boundary site
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) shownoncog                   ! display a representative fraction of non-cognate B cells
read(nfcell,*) noncog_display_fraction		! fraction of non-conate B cells to display
read(nfcell,*) fixedfile					! file with "fixed" parameter values
close(nfcell)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
    if (itestcase == 1 .or. itestcase == 2) then
        usetraffic = 0
        BASE_NFDC = 0
        base_exit_prob = 0
        iuse(EBI2) = 0
        iuse(CXCR5) = 0
    endif 
    if (itestcase == 2) then
        shownoncog = 1
    endif
endif

! Chemokines
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
if (BASE_NFDC == 0) then
    use_FDCs = .false.
else
    use_FDCs = .true.
endif
if (BASE_NMRC == 0) then
    use_MRCs = .false.
else
    use_MRCs = .true.
endif
if (.not.use_FDCs .and. .not.use_MRCs) then
    receptor(CXCR5)%used = .false.
    chemo(CXCL13)%used = .false.
endif
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
if (computedoutflow == 1) then
	computed_outflow = .true.
else
	computed_outflow = .false.
endif
if (shownoncog == 1) then
	show_noncognate = .true.
else
	show_noncognate = .false.
endif

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
write(nfout,*) 'computed_outflow: ',computed_outflow
write(nfout,*) 'use_cognate: ',use_cognate
write(nfout,*) 'random_cognate: ',random_cognate
write(nfout,*) 'use_ode_diffusion: ',use_ode_diffusion
write(nfout,*) 'NDIFFSTEPS: ',NDIFFSTEPS
write(nfout,*) 'VEGF_MODEL: ',VEGF_MODEL
write(nfout,*) 'VEGF_alpha: ',VEGF_alpha
write(nfout,*) 'VEGF_beta: ',VEGF_beta
write(nfout,*) 'VEGF_decayrate: ',VEGF_decayrate
write(nfout,*) 'vasc_maxrate: ',vasc_maxrate
write(nfout,*) 'vasc_decayrate: ',vasc_decayrate
write(nfout,*) 'vasc_beta: ',vasc_beta
write(nfout,*) 'vasc_n: ',vasc_n
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
	write(nfout,'(a,i2,f8.2)') 'sign, relative strength: ',p%sign,p%strength
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
! When a T cell dies:
!   %exists -> .false.
!   if non-cognate %ID -> 0 (i.e. it is marked for removal from the cell list)
!   Note that cognate cells are left in the list for tracing lineage/fate.
! In both cases it is removed from occupancy()%indx()
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine BcellDeath(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2)
logical :: cognate

cognate = (associated(cellist(kcell)%cptr))
cellist(kcell)%exists = .false.
totalres%dN_Dead = totalres%dN_Dead + 1
totalres%N_Dead = totalres%N_Dead + 1
if (cognate) then
	call set_stage(cellist(kcell)%cptr,DEAD)
	call set_region(cellist(kcell)%cptr,GONE)
endif
if (use_gaplist .and. .not.cognate) then
	cellist(kcell)%ID = 0
	ngaps = ngaps + 1
	if (ngaps > max_ngaps) then
		write(logmsg,*) 'BcellDeath: ngaps > max_ngaps'
		call logger(logmsg)
		stop
	endif
	gaplist(ngaps) = kcell
endif

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
integer :: kpar = 0
real(DP) :: R

ok = .true.
tnow = istep*DELTA_T
p1 => cellist(kcell)%cptr
gen = get_generation(p1)
region = get_region(p1)
if (gen == BC_MAX_GEN) then
    write(logmsg,*) 'cell_division: reached maximum generation: ',kcell
	call logger(logmsg)
    return
endif
if (gen == 1) then
    write(logmsg,*) 'First division at hour: ',tnow/60
    call logger(logmsg)
else
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
		p1%stagetime = tnow + DivisionTime(2)
	elseif (R < prob_gcc) then
!		call set_status(p1,BCL6_HI)
		call set_stage(p1,GCC_COMMIT)
		p1%stagetime = tnow
	else
		call set_stage(p1,PLASMA)
		call ReceptorLevel(kcell,ACTIVATED_TAG,ANTIGEN_TAG,1.0,0.0,1.0,cellist(kcell)%receptor_level)	! instantaneous change
		p1%stagetime = tnow
	endif
else
	call set_stage(p1,DIVIDING)
	p1%stagetime = tnow + DivisionTime(2)
endif
ctype = cellist(kcell)%ctype
site = cellist(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx

call CreateBcell(icnew,cellist(icnew),site2,ctype,gen,DIVIDING,status,region,ok)
if (.not.ok) return

p2 => cellist(icnew)%cptr

p2%generation = p1%generation
p2%region = p1%region
p2%ID = cellist(kcell)%ID	! progeny cells keep the cog_type ID - this enables lineage to be tracked
cellist(icnew)%receptor_level = cellist(kcell)%receptor_level
cellist(icnew)%entrytime = tnow
if (status == BCL6_LO) then
	prob_gcc = get_GccProb(gen)
	prob_plasma = get_PlasmaProb(gen) 
	R = par_uni(kpar)	
	if (R > prob_gcc + prob_plasma) then
		call set_stage(p2,DIVIDING)
		p2%stagetime = tnow + DivisionTime(2)
	elseif (R < prob_gcc) then
		call set_stage(p2,GCC_COMMIT)
		p2%stagetime = tnow + T_BCL6_UP
	else
		call set_stage(p2,PLASMA)
		call ReceptorLevel(icnew,ACTIVATED_TAG,PLASMA_TAG,1.0,0.0,1.0,cellist(icnew)%receptor_level)		
		p2%stagetime = tnow
	endif
else
	call set_stage(p2,DIVIDING)
	p2%stagetime = tnow + DivisionTime(2)
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
! Currently status just has two values: BCL6_LO, BCL6_HI
!-----------------------------------------------------------------------------------------
subroutine CreateBcell(kcell,cell,site,ctype,gen,stage,status,region,ok)
type(cell_type) :: cell
integer :: kcell, site(3), ctype, gen, stage, status, region
logical :: ok
integer :: stype, cogID, i
real :: tnow, param1, param2
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
stype = struct_type(ctype)
cell%exists = .true.
cell%entrytime = tnow
if (stype == NONCOG_TYPE_TAG) then
    if (associated(cell%cptr)) then
		write(*,*) 'Error: CreateBcell: cptr already associated for non-cognate cell: ',kcell
		stop
    endif
elseif (stype == COG_TYPE_TAG) then
    if (associated(cell%cptr)) then
        deallocate(cell%cptr)
    endif
    allocate(cell%cptr)
    param1 = log(BC_AVIDITY_MEDIAN)
    param2 = log(BC_AVIDITY_SHAPE)
    call set_generation(cell%cptr,gen)
	call set_stage(cell%cptr,stage)
	call set_status(cell%cptr,status)
	call set_region(cell%cptr,region)
    cell%cptr%status = BCL6_LO	! default
    cell%cptr%dietime = BIG_TIME
    cell%cptr%dividetime = BIG_TIME
    cell%cptr%stagetime = 0

! Maintain cognate_list at start or if we are running on a single node
! Otherwise cogID and cognate_list is maintained by make_cognate_list
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
else
    write(logmsg,*) 'ERROR: CreateBcell: bad ctype: ',ctype
    call logger(logmsg)
    stop
endif
cell%receptor_level = receptor%level(NAIVE_TAG)
cell%receptor_saturation_time = 0
lastID = lastID + 1     ! Each node makes its own numbers, with staggered offset
cell%ID = lastID
if (stype == COG_TYPE_TAG) then		! by default the cog_type ID is the same as naive cell ID
	cell%cptr%ID = cell%ID			! this is overridden when a cell is created by cell division
	if (logID == 0) then
		logID = cell%ID
	endif
endif
cell%site = site
cell%ctype = ctype
cell%step = 0
cell%lastdir = random_int(1,6,kpar)
end subroutine

!--------------------------------------------------------------------------------
! Add a cell (kcell) with characteristics (ctype, gen, stage) at site.
!--------------------------------------------------------------------------------
subroutine AddBcell(site,ctype,gen,stage,status,region,kcell,ok)
integer :: site(3), ctype, gen, stage, status, region, kcell
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
call CreateBcell(kcell,cellist(kcell),site,ctype,gen,stage,status,region,ok)
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
integer :: idc, kdc, k, x2, y2, z2, gen, stage, status, region
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
logID = 0
tnow = 0
nlist = 0
do x = 1,NX
    do y = 1,NY
	    do z = 1,NZ
            occupancy(x,y,z)%indx = 0
            occupancy(x,y,z)%FDC_nbdry = 0
            occupancy(x,y,z)%MRC_nbdry = 0
            nullify(occupancy(x,y,z)%bdry)
            site = (/x,y,z/)
			if (.not.InsideEllipsoid(site,Centre,Radius)) then
				occupancy(x,y,z)%indx = OUTSIDE_TAG
			else
			    nlist = nlist+1
			endif
        enddo
    enddo
enddo

call placeMRCs(ok)
call placeFDCs(ok)
if (.not.ok) stop

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
	    do z = 1,NZ
	        if (occupancy(x,y,z)%indx(1) == 0) then ! vacant site, not OUTSIDE_TAG or DC
                id = id+1
                lastID = id
                k = permc(id)
                site = (/x,y,z/)
                gen = 1
                stage = NAIVE
                status = BCL6_LO
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
                call CreateBcell(k,cellist(k),site,ctype,gen,stage,status,region,ok)
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
!call make_cognate_list(ok)		! do not use this now - all cognate cells are kept in the lists (for fate tracking)
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
call SetRadius(Nsites)
scale_factor = real(NBC_LN)*NLN_RESPONSE/NBcells0

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
r0 = Radius%x
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


!--------------------------------------------------------------------------------------
! The receptor levels are ramped over t1 < t < t2 from levels given by:
!   receptor(ireceptor)%level(ifrom)
! to those given by:
!  receptor(ireceptor)%level(ito)
! where ifrom and ito are cognate B cell stages of activation/differentiation:
! NAIVE_TAG
! ANTIGEN_TAG
! ACTIVATED_TAG
! GCC_TAG
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
!write(*,*) 'ReceptorLevel: ',t,t1,t2,alpha
do ireceptor = 1,MAX_RECEPTOR
	if (receptor(ireceptor)%used) then
		rfrom = receptor(ireceptor)%level(ifrom)
		rto = receptor(ireceptor)%level(ito)
		level(ireceptor) = (1-alpha)*rfrom + alpha*rto
!		if (ito == ANTIGEN_TAG) then
!			write(*,*) 'ReceptorLevel: ',ireceptor,rfrom,rto,alpha
!		endif
	endif
enddo
end subroutine

!---------------------------------------------------------------------
! For cognate cells, updates B cell state.
!
! Treatment of cell death
!------------------------
! There are two possibilities:
! (1) Specify a time-to-die for each cell (i.e. on cell division)
! (2) Use a probability of death in every time step (i.e. a rate of death)
! Previously (1) was implemented.  Now implement (2) for activated
! cognate cells, BCL6lo and BCL6hi, with generation-dependent rates.
!---------------------------------------------------------------------
subroutine Updater(ok)
logical :: ok
integer :: kcell, stage, iseq, tag, kfrom, kto, k, ncog, ntot, ndie
integer :: site(3), site2(3), freeslot, indx(2), status, DC(2), idc
integer :: tmplastcogID
real :: tnow, dstim, S, cyt_conc, mols_pM, Ctemp, dstimrate, stimrate
real :: t1, t2, die_prob
logical :: divide_flag, producing, first, dbg, unbound, flag, flag1
real(DP) :: R
integer :: kpar = 0
type(cog_type), pointer :: p

ok = .true.
dbg = .false.
flag = .false.
flag1 = .false.
ntot = 0
ncog = 0
ndie = 0
tnow = istep*DELTA_T
! Only cognate cells considered
tmplastcogID = lastcogID
do k = 1,tmplastcogID
    kcell = cognate_list(k)
    if (kcell == 0) then
		write(*,*) 'This can never happen'
        cycle
    endif
    p => cellist(kcell)%cptr
    if (.not.associated(p)) then
        write(logmsg,*) 'ERROR: updater: p not associated: ',k,kcell
	    call logger(logmsg)
        stop
    endif
    stage = get_stage(p)
    if (stage == LEFT .or. stage == DEAD) cycle
    ncog = ncog + 1

	! Cell death
	if (stage == FINISHED) then
        call BcellDeath(kcell)
        ndie = ndie + 1
        cycle
    endif
	die_prob = DeathProbability(p)
	R = par_uni(kpar)
	if (R < die_prob) then
        call BcellDeath(kcell)
        ndie = ndie + 1
        cycle
    endif

	if (stage == ANTIGEN_MET) then
		t2 = p%stagetime
		t1 = t2 - T_CCR7_UP
!		write(*,*) 'updater: tnow,t1,t2: ',kcell,tnow,t1,t2,NAIVE_TAG,ANTIGEN_TAG
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
if (dbug .and. ndie > 0) write(nfout,*) 'ndie: ',ndie
end subroutine

!--------------------------------------------------------------------------------------
! stage is the current stage of the cell 
! p%nextstagetime = time that the cell is expected to make transition to the next stage.
! When the current time tnow reaches (exceeds) nextstagetime, the case corresponding to
! stage is executed to handle the transition to the next stage.
! The interval between stage transitions can be:
!    Fixed length 
!      ANTIGEN_MET -> CCR7_UP takes T_CCR7_UP
!      TCELL_MET -> EBI2_UP takes T_EBI2_UP
!      GCC_COMMIT -> BCL6_UP takes T_BCL6_UP
!    Variable length determined by stochastic cell motion
!      NAIVE -> ANTIGEN_MET depends on antigen encounter at the upper follicle boundary
!      CCR7_UP -> TCELL_MET depends on helper T cell encounter at the lower follicle boundary
!    Variable length determined by a probability distribution
!      EBI2_UP -> DIVIDING, BCL6_UP -> DIVIDING, DIVIDING -> DIVIDING determined by division 
!      time distribution.
!--------------------------------------------------------------------------------------
subroutine UpdateStage(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, gen, ctype, site(3)
real :: nextstagetime
type(cog_type), pointer :: p

divide_flag = .false.
site = cellist(kcell)%site
p => cellist(kcell)%cptr
stage = get_stage(p)
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
nextstagetime = p%stagetime
if (tnow > nextstagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)
		if (AntigenEncounter(kcell,site)) then
	        call set_stage(p,ANTIGEN_MET)
	        p%stagetime = tnow + T_CCR7_UP		! time of execution of case (ANTIGEN_MET), which sets the stage to CCR7_UP
	        write(logmsg,'(a,i6,2f8.2)') 'updatestage: -> ANTIGEN_MET: ',kcell,tnow,p%stagetime	
	        call logger(logmsg)
		endif		
	case (ANTIGEN_MET)    ! transition from ANTIGEN_MET to CCR7_UP
        call set_stage(p,CCR7_UP)
		p%stagetime = tnow
	case (CCR7_UP)     ! possible transition from CCR7_UP to TCELL_MET
		if (TCellEncounter(kcell,site)) then
	        call set_stage(p,TCELL_MET)
	        p%stagetime = tnow + T_EBI2_UP		! time of execution of case (TCELL_MET), which sets the stage to EBI2_UP
		endif
	case (TCELL_MET)
		call set_stage(p,EBI2_UP)
        p%stagetime = tnow		
	case (EBI2_UP)
		call set_stage(p,DIVIDING)
		p%stagetime = tnow + DivisionTime(1)		! time of execution of case (DIVIDING), which sets divide_flag
	case (DIVIDING)								! the next stage for DIVIDING cells is set in cell_division()
        gen = get_generation(p)
		if (gen == BC_MAX_GEN) then
            call set_stage(p,FINISHED)
			p%stagetime = BIG_TIME
		else
		    divide_flag = .true.				! this triggers a call to cell_division()
		endif
	case (GCC_COMMIT)
        call set_stage(p,BCL6_UP)
		p%stagetime = tnow + T_BCL6_UP
	case (BCL6_UP)
		call set_status(p,BCL6_HI)
		call set_stage(p,DIVIDING)
        p%stagetime = tnow + max(0.,DivisionTime(1) - T_BCL6_UP)		
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
r2 = (x/Radius%x)**2 + (y/Radius%y)**2+ (z/Radius%z)**2
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
r2 = (x/Radius%x)**2 + (y/Radius%y)**2+ (z/Radius%z)**2
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
