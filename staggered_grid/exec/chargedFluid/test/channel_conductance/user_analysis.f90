module user_analysis
  ! This module contains user-specific routines for analyzing data
  ! Right not this only analyses fluxes on a planar cut through the domain

  use bl_types
  use bl_error_module
  use bl_IO_module
  use parallel
  use probin_common_module
  use probin_charged_module
  
  implicit none

  private

  integer :: flux_unit, velocity_unit

  public :: analyze_planar_cut, initialize_planar_cut, destroy_planar_cut

contains

  subroutine initialize_planar_cut()
    ! One can use this to open files for output, allocate/initialize averages over time, etc.
       
    flux_unit = unit_new()
    open(unit=flux_unit, file = "fluxes.dat", status = 'unknown', action = 'write')
  
  end subroutine

  subroutine analyze_planar_cut(cut,id)
    real(kind=dp_t), intent(in) :: cut(0:,0:,1:) ! Cut through a plane, of dimension (Ny,Nz,nspecies)
    integer, intent(in) :: id ! ID to print, for example, time step

    integer :: i,j,k,step
    real(kind=dp_t) :: species_fluxes(nspecies), mass_flux, current
    real(kind=dp_t) :: velocities(n_cells(2)) ! Assume velocity is a function of y only (slit channel or 2D)
    character(len=32) :: id_string
    
    step=id-1 ! We associate this with the beginning of the step so that the initial condition is written out
    
    if (parallel_IOProcessor() .and. &
       (stats_int > 0) .and. (mod(step,stats_int)==0)) then
    
       do i=1,nspecies ! Compute average species fluxes
          species_fluxes(i)=sum(cut(:,:,i))/size(cut(:,:,i))
       end do   
       mass_flux=sum(species_fluxes)

       if (use_charged_fluid) then
          current = sum(species_fluxes*charge_per_mass(1:nspecies))
          write(flux_unit,'(1000g17.9)') step*dt_saved, current, mass_flux/rho0, species_fluxes
       else
          write(flux_unit,'(1000g17.9)') step*dt_saved, mass_flux/rho0, species_fluxes
       end if

       if(.true.) then
          ! Now compute the velocity as a function of y to check
          velocity_unit = unit_new()
          write(id_string,"(I8.8)") step
          open(unit=velocity_unit, file = "velocities-"//trim(ADJUSTL(id_string))//".dat", status = 'unknown', action = 'write')

          do j=0,size(cut,1)-1
             mass_flux=sum(cut(j,:,:))/size(cut,dim=2)
             write(velocity_unit,'(1000g17.9)') j*dx_saved(2), mass_flux/rho0 
          end do
       end if   
    
    end if

  end subroutine analyze_planar_cut

  subroutine destroy_planar_cut()
    ! One can use this to close files for output, deallocate averages over time, etc.
  
    close(flux_unit)
  
  end subroutine

end module user_analysis
