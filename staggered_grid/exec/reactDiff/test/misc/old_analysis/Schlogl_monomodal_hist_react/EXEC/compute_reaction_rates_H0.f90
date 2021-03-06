module compute_reaction_rates_module

  use bl_types
  use bl_error_module
  use probin_common_module, only: nspecies
  use probin_chemistry_module, only: nreactions, stoichiometric_factors, rate_const, rate_multiplier, &
                                     include_discrete_LMA_correction, exclude_solvent_comput_rates

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
    real(kind=dp_t) :: tmp

    n_nonneg(1:nspecies) = max(0.d0,n_in(1:nspecies))

    if ((.not. include_discrete_LMA_correction) .and. (exclude_solvent_comput_rates .eq. 0)) then
       ! Use traditional LMA without accounting for discrete/integer nature of the molecules involved
       ! no species is excluded for calculating rates
       do reaction=1, nreactions
          ! This works since raising to the zeroth power gives unity:
          reaction_rates(reaction) = rate_multiplier*rate_const(reaction)*&
               product(n_nonneg(1:nspecies)**stoichiometric_factors(1:nspecies,1,reaction))
       end do
    else
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
                case(0) ! Species does not participate in reaction

                case(1) ! Rate ~ N
                   tmp = n_nonneg(species)*dv-1.d0

                   if (tmp .ge. 1.d0) then
                      reaction_rates(reaction) = reaction_rates(reaction)*n_nonneg(species)
                   else if (tmp .ge. 0.d0) then
                      !tmp = (3.d0-2.d0*tmp)*tmp**2  ! H1
                      !tmp = (10.d0-15.d0*tmp+6.d0*tmp**2)*tmp**3  ! H2
                      reaction_rates(reaction) = reaction_rates(reaction)*n_nonneg(species)*tmp
                   else
                      reaction_rates(reaction) = 0.d0
                   end if

                case(2) ! Rate ~ N*(N-1)
                   tmp = n_nonneg(species)*dv-2.d0

                   if (tmp .ge. 1.d0) then
                      reaction_rates(reaction) = reaction_rates(reaction)*n_nonneg(species)*(n_nonneg(species)-1.0d0/dv)
                   else if (tmp .ge. 0.d0) then
                      !tmp = (3.d0-2.d0*tmp)*tmp**2  ! H1
                      !tmp = (10.d0-15.d0*tmp+6.d0*tmp**2)*tmp**3  ! H2
                      reaction_rates(reaction) = reaction_rates(reaction)*n_nonneg(species)*(n_nonneg(species)-1.0d0/dv)*tmp
                   else
                      reaction_rates(reaction) = 0.d0
                   end if

                case(3) ! Rate ~ N*(N-1)*(N-2)
                   tmp = n_nonneg(species)*dv-3.d0

                   if (tmp .ge. 1.d0) then 
                      reaction_rates(reaction) = &
                           reaction_rates(reaction)*n_nonneg(species)*(n_nonneg(species)-1.0d0/dv) &
                           *(n_nonneg(species)-2.0d0/dv)
                   else if (tmp .ge. 0.d0) then
                      !tmp = (3.d0-2.d0*tmp)*tmp**2  ! H1
                      !tmp = (10.d0-15.d0*tmp+6.d0*tmp**2)*tmp**3  ! H2
                      reaction_rates(reaction) = &
                           reaction_rates(reaction)*n_nonneg(species)*(n_nonneg(species)-1.0d0/dv) &
                           *(n_nonneg(species)-2.0d0/dv)*tmp
                   else
                      reaction_rates(reaction) = 0.d0
                   end if

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
