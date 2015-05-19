module probin_energy_module

  use bl_types
  use bl_IO_module
  use bl_constants_module

  implicit none

  integer, parameter  :: MAXSPECIES = 20

  REAL*8 , SAVE :: dia_in(MAXSPECIES)          ! diameters of molecules
  REAL*8 , SAVE :: int_deg_free_in(MAXSPECIES) ! internal degrees of freedom
  REAL*8 , SAVE :: p0_in                       ! initial background pressure, p0
  integer, SAVE :: heating_type                ! form of external heating
  integer, SAVE :: use_fake_diff               ! force diffusion coefficients to constant
  REAL*8 , SAVE :: fake_diff_coeff             ! Adjust transport coefficients artificially
  REAL*8 , SAVE :: fake_soret_factor           ! Adjust transport coefficients artificially
  integer, SAVE :: dpdt_iters                  ! Loops over volume discrepancy correction
  REAL*8 , SAVE :: dpdt_factor                 ! scaling factor for volume discrepancy correction
  integer, SAVE :: deltaT_iters                ! Loops over deltaT iterative solve

  NAMELIST /probin_energy/ dia_in
  NAMELIST /probin_energy/ int_deg_free_in
  NAMELIST /probin_energy/ p0_in
  NAMELIST /probin_energy/ heating_type
  NAMELIST /probin_energy/ use_fake_diff
  NAMELIST /probin_energy/ fake_diff_coeff
  NAMELIST /probin_energy/ fake_soret_factor
  NAMELIST /probin_energy/ dpdt_iters
  NAMELIST /probin_energy/ dpdt_factor
  NAMELIST /probin_energy/ deltaT_iters

contains

  subroutine probin_energy_init()

    integer :: iwrk, nfit, i, ic, ii, j
    integer :: dochem, dostrang
    double precision :: rwrk
    integer ns
    real*8 :: mu,Fij,Fijstar,fact1
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs

    narg = command_argument_count()

    ! default values to be replace from inputs file or command line
    dia_in = 0.d0
    p0_in = 1013250.d0
    heating_type = 0
    int_deg_free_in = 0
    use_fake_diff = 0
    fake_diff_coeff = -1.d0
    fake_soret_factor = 1.d0
    dpdt_iters = 3
    dpdt_factor = 2.d0
    deltaT_iters = 3

    ! read from input file 
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_energy)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--p0_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) p0_in

       case ('--heating_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) heating_type

       case ('--use_fake_diff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_fake_diff

       case ('--fake_diff_coeff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fake_diff_coeff

       case ('--fake_soret_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fake_soret_factor

       case ('--dpdt_iters')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dpdt_iters

       case ('--dpdt_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dpdt_factor

       case ('--deltaT_iters')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) deltaT_iters

       case ('--')
          farg = farg + 1
          exit

       case default

       end select
       farg = farg + 1
    end do

  end subroutine probin_energy_init

end module probin_energy_module
