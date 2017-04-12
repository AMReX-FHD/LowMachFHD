module restart_module
  
  ! this is just a placeholder - you need to write your own custom restart.f90 and
  ! put the local copy in the example source directory

contains

  subroutine restart()

    use bl_error_module

    call bl_error("restart: write a custom routine and put in local source directory")
    
  end subroutine restart

end module restart_module
