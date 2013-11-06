module populate_DbarGama_module

  use probin_multispecies_module
 
  implicit none

  private

  public :: populate_DbarGama
  
contains
  
  subroutine populate_DbarGama(Dbar, Gama, mass)
 
    real(kind=dp_t), intent(inout) :: Dbar(:,:)  ! SM diffusion constants 
    real(kind=dp_t), intent(inout) :: Gama(:,:)  ! non-ideality coefficient 
    real(kind=dp_t), intent(inout) :: mass(:)    ! mass of each species
    integer                        :: n,row,column
 
    ! populate Dbar, Gama & mass; for initial case doesn't change in each cell. 
    n=0; 
    do row=1, nspecies  
       do column=1, row-1
          n=n+1
          Dbar(row, column) = Dbar_in(n)
          Dbar(column, row) = Dbar(row, column) ! symmetric
          Gama(row, column) = 0.d0       
          Gama(column, row) = Gama(row, column) ! symmetric
       enddo
       Dbar(row, row) = 0.d0      ! self-diffusion is zero
       Gama(row, row) = 1.d0      ! set to unit matrix for time being
       mass(row)      = molmass_in(row) ! populate species mass 
    enddo

  end subroutine populate_DbarGama

end module populate_DbarGama_module
