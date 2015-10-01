module compute_reaction_rates_module

  use bl_types
  use bl_error_module
  use probin_reactdiff_module, only : nspecies, nreactions, stoichiometric_factors, chemical_rates, include_discrete_LMA_correction
  implicit none

  private

  public :: compute_reaction_rates

contains

  subroutine compute_reaction_rates(n_in,reaction_rates,dv)
    
    real(kind=dp_t), intent(in   ) :: n_in(:),dv
    real(kind=dp_t), intent(inout) :: reaction_rates(:)

    integer :: reaction, species
    
    if(include_discrete_LMA_correction) then
       ! Use traditional LMA but correct for the fact that for binary reactions rate ~ N*(N-1) and not N^2, etc.
       do reaction=1, nreactions
          reaction_rates(reaction) = chemical_rates(reaction)
          do species=1, nspecies
             select case(stoichiometric_factors(species,1,reaction))
             case(0)
                ! Species does not participate in reaction
             case(1)
                ! Rate ~ N, where N is number of molecules
                reaction_rates(reaction)=reaction_rates(reaction)*n_in(species)
             case(2)
                ! Rate ~ N*(N-1) where N is number of molecules
                reaction_rates(reaction)=reaction_rates(reaction)*n_in(species)*(n_in(species)-1.0d0/dv) 
             case(3)   
                ! Rate ~ N*(N-1)*(N-2) where N is number of molecules
                reaction_rates(reaction)=reaction_rates(reaction)*n_in(species)*(n_in(species)-1.0d0/dv)*(n_in(species)-2.0d0/dv)
             case default
                ! This is essentially impossible in practice and won't happen
                call bl_error("Stochiometric coefficients larger then 3 not supported")      
             end select
          end do
       end do
    else
       ! Use traditional LMA without accounting for discrete/integer nature of the molecules involved
       do reaction=1, nreactions
          ! This works since raising to the zeroth power gives unity:
          reaction_rates(reaction) = &
               chemical_rates(reaction)*product(n_in(1:nspecies)**stoichiometric_factors(1:nspecies,1,reaction))
       end do
    end if
    
  end subroutine compute_reaction_rates

end module compute_reaction_rates_module
