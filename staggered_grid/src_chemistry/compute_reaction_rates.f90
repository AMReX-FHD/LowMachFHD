module compute_reaction_rates_module

  use bl_types
  use bl_error_module
  use probin_common_module, only: nspecies
  use probin_chemistry_module, only: nreactions, stoichiometric_factors, rate_const, rate_multiplier, &
                                     include_discrete_LMA_correction, exclude_solvent_comput_rates, &
                                     use_mole_frac_LMA

  implicit none

  private

  public :: compute_reaction_rates

contains

  ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
  subroutine compute_reaction_rates(n_in,reaction_rates,dv)

    real(kind=dp_t), intent(in   ) :: n_in(:),dv
    real(kind=dp_t), intent(inout) :: reaction_rates(:)

    integer :: reaction, species

    real(kind=dp_t) :: n_nonneg(nspecies)
    real(kind=dp_t) :: n_sum, mole_frac

    n_nonneg(1:nspecies) = max(0.d0,n_in(1:nspecies))
    n_sum = sum(n_nonneg(1:nspecies))
    if(n_sum < 0.0d0) n_sum=1.0d0/dv 
    
    if (use_mole_frac_LMA .and. include_discrete_LMA_correction) then
      ! Use mole-fraction based LMA (general ideal mixtures) with integer corrections
      
      do reaction=1, nreactions
        reaction_rates(reaction) = rate_multiplier*rate_const(reaction)
        do species=1, nspecies
          ! Donev: Replaced case statement by if here
          ! Donev: Made sure n_sum is never zero for empty cells to avoid division by zero
          
          if(stoichiometric_factors(species,1,reaction)>=1) then
            ! rate ~ N/N_sum
            if(n_nonneg(species)>0.0d0) then ! This species is present in this cell
               reaction_rates(reaction) = reaction_rates(reaction) * n_nonneg(species)/n_sum
            else
               reaction_rates(reaction) = 0.0d0  
            end if            
          end if
          if(stoichiometric_factors(species,1,reaction)>=2) then
            ! rate ~ (N/N_sum)*((N-1)/(N_sum-1))
            ! Donev: Avoid division by zero or negative rates
            if(n_nonneg(species)>1.0d0/dv) then ! There is at least one molecule of this species in this cell
               reaction_rates(reaction) = reaction_rates(reaction) * (n_nonneg(species)-1.0d0/dv)/(n_sum-1.0d0/dv)
            else
               reaction_rates(reaction) = 0.0d0  
            end if  
          end if
          if(stoichiometric_factors(species,1,reaction)>=3) then ! Donev added ternary reactions here
            ! rate ~ (N/N_sum)*((N-1)/(N_sum-1))*((N-2)/(N_sum-2))
            if(n_nonneg(species)>2.0d0/dv) then ! There is at least two molecules of this species in this cell
              reaction_rates(reaction) = reaction_rates(reaction) * (n_nonneg(species)-2.0d0/dv)/(n_sum-2.0d0/dv)
            else
               reaction_rates(reaction) = 0.0d0  
            end if
          end if    
          if(stoichiometric_factors(species,1,reaction)>=4) then
            ! This is essentially impossible in practice and won't happen
            call bl_error("Stochiometric coefficients larger then 3 not supported")
          end if
        end do
      end do
      
    else if ((.not. include_discrete_LMA_correction) .and. (exclude_solvent_comput_rates .eq. 0)) then
       ! Use array syntax to compute the answer faster without all the if statements
    
       ! if mole fraction based LMA is used, convert number densities into mole fractions
       if (use_mole_frac_LMA) then
          n_nonneg(1:nspecies) = n_nonneg(1:nspecies)/n_sum ! Donev: Changed denominator to n_sum here
       end if

       ! Use traditional LMA without accounting for discrete/integer nature of the molecules involved
       ! no species is excluded for calculating rates
       do reaction=1, nreactions
          ! This works since raising to the zeroth power gives unity:
          reaction_rates(reaction) = rate_multiplier*rate_const(reaction)*&
               product(n_nonneg(1:nspecies)**stoichiometric_factors(1:nspecies,1,reaction))
       end do
       
    else ! General case of number-density based LMA is handled by slower code that includes species by species
    
       do reaction=1, nreactions
          !write(*,*) "reaction=", reaction, " rate_const=", rate_const(reaction), &
          !  " stochiometry=", stoichiometric_factors(1:nspecies,1,reaction)
          reaction_rates(reaction) = rate_multiplier*rate_const(reaction)

          do species=1, nspecies
             if (species .eq. exclude_solvent_comput_rates) then
             ! if solvent, skip
             ! for exclude_solvent_comput_rates=0 (default), this never happens
                cycle
             end if

             if (include_discrete_LMA_correction) then
                ! Use traditional LMA but correct for the fact that for binary reactions rate ~ N*(N-1) and not N^2, etc.,
                ! where N is the total number of molecules
                select case(stoichiometric_factors(species,1,reaction))
                case(0)
                   ! Species does not participate in reaction
                case(1)
                   ! Rate ~ N
                   reaction_rates(reaction) = &
                        reaction_rates(reaction)*n_nonneg(species)
                case(2)
                   ! Rate ~ N*(N-1)
                   reaction_rates(reaction) = &
                        reaction_rates(reaction)*n_nonneg(species)*max(n_nonneg(species)-1.0d0/dv,0.d0)
                case(3)   
                   ! Rate ~ N*(N-1)*(N-2)
                   reaction_rates(reaction) = &
                        reaction_rates(reaction)*n_nonneg(species)*max(n_nonneg(species)-1.0d0/dv,0.d0) &
                        *max(n_nonneg(species)-2.0d0/dv,0.d0)
                case default
                   ! This is essentially impossible in practice and won't happen
                   call bl_error("Stochiometric coefficients larger then 3 not supported")      
                end select
             else
                reaction_rates(reaction)= reaction_rates(reaction)*(n_nonneg(species)**stoichiometric_factors(species,1,reaction))
             end if
          end do
          !write(*,*) "reaction=", reaction, " rate=", reaction_rates(reaction)
       end do
       
    end if
    
  end subroutine compute_reaction_rates

end module compute_reaction_rates_module
