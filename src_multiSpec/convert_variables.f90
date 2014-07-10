module convert_variables_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: convert_rho_to_c
  
contains

  subroutine convert_rho_to_c(mla,rho,rho_tot,c,rho_to_c)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(in   ) :: rho_tot(:)
    type(multifab) , intent(inout) ::   c(:)
    logical        , intent(in   ) :: rho_to_c

    ! local
    integer :: n,nlevs,i

    nlevs = mla%nlevel

    if (rho_to_c) then

       ! rho to c - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(c(n),1,rho(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_c(c(n),i,rho_tot(n),1,1,0)
          end do
       end do

    else

       ! c to rho - VALID + GHOST (CAN CHANGE TO DO ONLY GHOST TO SAVE COMPUTATION)
       do n=1,nlevs
          call multifab_copy_c(rho(n),1,c(n),1,nspecies,rho(n)%ng)
          do i=1,nspecies
             call multifab_mult_mult_c(rho(n),i,rho_tot(n),1,1,rho(n)%ng)
          end do
       end do

    end if

  end subroutine convert_rho_to_c

end module convert_variables_module
