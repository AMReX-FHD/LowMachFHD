module probin_multispecies_module

  use bl_types
  use bl_space

  implicit none

  integer,save    :: nspecies,max_step
  namelist /probin_multispecies/ nspecies
  namelist /probin_multispecies/ max_step    

contains

  subroutine probin_multispecies_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer    :: narg, farg
    character(len=128) :: fname
    logical :: lexist
    logical :: need_inputs
    integer :: un

    narg = command_argument_count()

    nspecies = 4
    max_step = 10000
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_multispecies)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! stuff that can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--nspecies')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nspecies
      case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step
       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
