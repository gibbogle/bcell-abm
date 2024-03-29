! Chemokine data

module chemokine

use global

implicit none

type receptor_type
	character(8) :: name
	logical :: used
	integer :: chemokine
	integer :: sign
	real(REAL_KIND) :: strength
	real(REAL_KIND) :: level(5)
	real(REAL_KIND) :: saturation_threshold
	real(REAL_KIND) :: refractory_time
end type

type chemokine_type
	character(8) :: name
	logical :: used
	logical :: use_secretion
	real(REAL_KIND) :: bdry_rate
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
	real(REAL_KIND), allocatable :: coef(:,:)
	real(REAL_KIND), allocatable :: conc(:,:,:)
	real(REAL_KIND), allocatable :: grad(:,:,:,:)
end type

type ODEdiff_type
	integer :: ichemo
	integer :: nvars
	integer, allocatable :: ivar(:,:,:)
	integer, allocatable :: varsite(:,:)
	integer, allocatable :: icoef(:,:)
end type

type(chemokine_type), target :: chemo(MAX_CHEMO)
type(receptor_type), target :: receptor(MAX_RECEPTOR)
type(ODEdiff_type) :: ODEdiff

end module