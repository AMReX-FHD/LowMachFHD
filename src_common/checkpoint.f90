module checkpoint_module
  
  ! this is just a placeholder - you need to write your own custom checkpoint.f90 and
  ! put the local copy in the example source directory

  subroutine checkpoint()

    use bl_error_module

    call bl_error("checkpoint: write a custom routine and put in local source directory")
    
  end subroutine checkpoint

end module checkpoint_module
