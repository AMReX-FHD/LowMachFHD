module advance_reaction_diffusion_module

  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use bc_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use chemical_rates_module
  use implicit_diffusion_module
  use probin_common_module, only: variance_coef_mass
  use probin_reactdiff_module, only: nspecies, D_Fick, temporal_integrator, &
       midpoint_stoch_flux_type, nreactions

  implicit none

  private

  public :: advance_reaction_diffusion

contains

  ! Donev: Added this documentation to explain what this does:
  ! this solves dn/dt = div ( D grad (n)) + div (sqrt(2*variance*D*n)*W) + f(n) - g
  !  where f(n) are the chemical production rates (deterministic or stochastic)
  !  and g=ext_src (note minus sign!) is a constant (in time) *deterministic* source term.
  ! To model stochastic particle production (sources) include g in the definition of f instead
  !  or add it as a reaction 0->products

  ! Donev: Here I kept ext_src optional and made sure it works without having to make a temporary multifab
  ! I would like to try to minimize the needless overheads in cases where programming is not made more complicated
  ! In this case it is very easy to do this -- all you do is add if(present(ext_src)) in front of a few lines
  ! Note that the same can be done even in advance_reaction and advance_diffusion but since we always call
  ! those with that term present I kept it there as not optional instead
  subroutine advance_reaction_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src(:)

    ! local
    ! Donev: I deleted here ext_src
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)
    type(multifab) :: rate1(mla%nlevel)
    type(multifab) :: rate2(mla%nlevel) 

    integer :: nlevs,dm,n,i,spec

    type(bl_prof_timer),save :: bpt

    ! Donev: I made this a parameter since it is fixed
    real(kind=dp_t), parameter :: mattingly_lin_comb_coef(1:2) = (/-1.d0, 2.d0/)

    !!!!!!!!
    ! init !
    !!!!!!!!

    ! Donev: I believe that for diffusion there is actually a limit of at least 4 cells, or maybe at least 3
    ! due to issues with ghost cells. Andy should confirm 
    ! If the system is smaller than the minimum for which diffusion works correctly abort here
    ! single cell case? 
    if ((multifab_volume(n_old(1))/nspecies)<=1) then
      ! Donev: There seems to be no point in doing the work, so just skip implementing this
      call bl_error("advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for single cell")
    end if
    ! Donev: This is not technically an error as the code will work, but better tell the user there is a more efficient way:
    if(nreactions<1) then
      call bl_error("advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for diffusion only")
    end if

    call build(bpt,"advance_reaction_diffusion")

    nlevs = mla%nlevel
    dm = mla%dim

    ! build
    do n=1,nlevs
      call multifab_build(diff_fluxdiv(n),mla%la(n),nspecies,0)
      call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0)
      do i=1,dm
        call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
      end do
      call multifab_build(rate1(n),mla%la(n),nspecies,0)

      if (temporal_integrator .eq. -2) then ! explitcit midpoint
        call multifab_build(rate2(n),mla%la(n),nspecies,0)
      end if
    end do

    ! diffusion coefficients (for now just setting each to a different constant)
    ! if one wants a space-dependent D or state-dependent D,
    ! see multispecies code as example
    do n=1,nlevs
      do i=1,dm
        do spec=1,nspecies
          call multifab_setval_c(diff_coef_face(n,i),D_Fick(spec),spec,1,all=.true.)
        end do
      end do
    end do

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! compute stochastic flux divergence
    if (variance_coef_mass .gt. 0.d0) then
      ! fill random flux multifabs with new random numbers
      call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
      ! compute stochastic flux divergence
      call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                the_bc_tower,increment_in=.false.)
    else
      do n=1,nlevs
        call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
      end do
    end if

    !!!!!!!!!!!!!!!!
    ! time advance !
    !!!!!!!!!!!!!!!!

    if (temporal_integrator .eq. -1) then  ! forward Euler

      ! calculate rates
      ! rates could be deterministic or stochastic depending on use_Poisson_rng
      call chemical_rates(mla,n_old,rate1,dx,dt)

      ! Donev: I added some documentation here, please check
      ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^n
      !                   + dt div (sqrt(2 D_k n_k^n dt) Z) ! Gaussian noise
      !                   + 1/dV * P( f(n_k)*dt*dV )        ! Poisson noise
      !                   + dt ext_src
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,rate1(n))
        if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else if (temporal_integrator .eq. -2) then  ! explicit midpoint

      !!!!!!!!!!!!!!!
      ! predictor   !
      !!!!!!!!!!!!!!!

      ! calculate rates from a(n_old)
      call chemical_rates(mla,n_old,rate1,dx,dt/2.d0)

      ! predictor
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt/2.d0,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
        if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt/2.d0,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

      !!!!!!!!!!!!!!!
      ! corrector !
      !!!!!!!!!!!!!!!

      ! Here we do not write this in the form that Mattingly et al do
      !  where we just continue the second half of the time step from where we left
      ! Rather, we compute terms at the midpoint and then add contributions from both halves of the time step to n_old
      ! This works simpler with diffusion but we have to store both rates1 and rates2

      ! compute diffusive flux divergence
      call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

      ! calculate rates from 2*a(n_pred)-a(n_old)
      call chemical_rates(mla,n_old,rate2,dx,dt/2.d0,n_new,mattingly_lin_comb_coef)

      ! compute stochastic flux divergence and add to the ones from the predictor stage
      if (variance_coef_mass .gt. 0.d0) then

        ! first, fill random flux multifabs with new random numbers
        call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

        ! compute n on faces to use in the stochastic flux in the corrector
        ! three possibilities
        select case (midpoint_stoch_flux_type)
        case (1)
          ! use n_old
          call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                    the_bc_tower,increment_in=.true.)
        case (2)
          ! use n_pred 
          call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                    the_bc_tower,increment_in=.true.)
        case (3)
          ! compute n_new=2*n_pred-n_old
          ! here we use n_new as temporary storage since it will be overwritten shortly
          do n=1,nlevs
            call multifab_mult_mult_s_c(n_new(n),1,2.d0,nspecies,n_new(n)%ng)
            call multifab_sub_sub_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
          end do
          ! use n_new=2*n_pred-n_old
          call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                    the_bc_tower,increment_in=.true.)
        case default
          call bl_error("advance_reaction_diffusion: invalid midpoint_stoch_flux_type")
        end select

      ! Donev: I deleted the else clause here since stochfluxdiv is already set to zero
      ! I did this to match what is in advance_diffusion since these two codes are copies of each other
      end if

      ! Donev: I added some documentation here, please check
      ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
      !                   + dt div (sqrt(2 D_k n_k^n dt) Z_1 / sqrt(2) ) ! Gaussian noise
      !                   + dt div (sqrt(2 D_k n_k^? dt) Z_2 / sqrt(2) ) ! Gaussian noise
      !                   + 1/dV * P_1( f(n_k)*dt*dV/2 )                 ! Poisson noise
      !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*dt*dV/2 ) ! Poisson noise
      !                   + dt ext_src
      ! where
      ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
      !       = n_k^pred            (midpoint_stoch_flux_type=2)
      !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
      
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate2(n))
        if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else
      call bl_error("advance_reaction_diffusion: invalid temporal_integrator")
    end if

    !!!!!!!!!!!
    ! destroy !
    !!!!!!!!!!!

    do n=1,nlevs
      call multifab_destroy(diff_fluxdiv(n))
      call multifab_destroy(stoch_fluxdiv(n))
      do i=1,dm
        call multifab_destroy(diff_coef_face(n,i))
      end do
      call multifab_destroy(rate1(n))

      if (temporal_integrator .eq. -2) then  ! explicit midpoint
        call multifab_destroy(rate2(n))
      end if
    end do

    call destroy(bpt)

  end subroutine advance_reaction_diffusion

end module advance_reaction_diffusion_module
