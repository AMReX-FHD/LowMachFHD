module user_analysis
  ! This module contains user-specific routines for analyzing data
  ! Right not this only analyses fluxes on a planar cut through the domain

  use bl_types
  use bl_error_module
  use parallel
  use probin_common_module
  
  implicit none

  private

  public :: analyze_planar_cut

contains

  subroutine initialize_planar_cut()
    call bl_warn("using default initialize_planar_cut")
    call bl_warn("each application code should have its own copy")
  
  end subroutine

  subroutine analyze_planar_cut(cut)

    real(kind=dp_t), allocatable :: cut(:,:) ! Cut through a plane

    integer :: i
    
    call bl_warn("using default analyze_planar_cut; printing average values to screen")
    call bl_warn("each application code should have its own copy")

    if (parallel_IOProcessor()) then      
      do i=1,size(cut,dim=1)
         write(*,*) cut(i,:)
      end do
    end if

  end subroutine analyze_planar_cut

  subroutine destroy_planar_cut()
  
    call bl_warn("using default destroy_planar_cut")
    call bl_warn("each application code should have its own copy")
  
  end subroutine

end module user_analysis
