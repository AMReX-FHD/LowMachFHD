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
  use probin_reactdiff_module, only: nspecies, D_Fick, temporal_integrator 

  implicit none

  private

  public :: advance_reaction_diffusion

contains

  subroutine advance_reaction_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src_in(:)

    ! local
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)
    type(multifab) :: ext_src(mla%nlevel)
    type(multifab) :: rate1(mla%nlevel)
    type(multifab) :: rate2(mla%nlevel) 

    integer :: nlevs,dm,n,i,spec

    type(bl_prof_timer),save :: bpt

    ! temporary 
    ! (later, midpoint_stoch_flux_type will be included in probin_reactdiff.)
    ! there are three options for updating the stochastic flux term of the corrector.
    ! midpoint_stoch_flux_type = 1: K_n * W1 + K_n * W2
    ! midpoint_stoch_flux_type = 2: K_n * W1 + K_pred * W2
    ! midpoint_stoch_flux_type = 3: K_n * W1 + K_{2*pred-n} * W2
    integer, parameter :: midpoint_stoch_flux_type = 1

    real(kind=dp_t) :: mattingly_lin_comb_coef(1:2)

    mattingly_lin_comb_coef(1) = -1.d0
    mattingly_lin_comb_coef(2) = 2.d0

    !!!!!!!!
    ! init !
    !!!!!!!!

    ! single cell case? 
    if ((multifab_volume(n_old(1))/nspecies)<=1) then
      call bl_error("advance_reaction_diffusion: single cell case has not been considered yet")
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
      call multifab_build(ext_src(n),mla%la(n),nspecies,0)
      call multifab_build(rate1(n),mla%la(n),nspecies,0)

      if (temporal_integrator .eq. -2) then ! explitcit midpoint
        call multifab_build(rate2(n),mla%la(n),nspecies,0)
      end if
    end do

    ! ext_src
    if (present(ext_src_in)) then
      do n=1,nlevs
        call multifab_copy_c(ext_src(n),1,ext_src_in(n),1,nspecies,0)
      end do
    else
      do n=1,nlevs
        call multifab_setval(ext_src(n),0.d0,all=.true.)
      end do
    end if

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

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,rate1(n))
        call multifab_saxpy_3(n_new(n),dt,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else if (temporal_integrator .eq. -2) then  ! explicit midpoint

      ! calculate rates from a(n_old)
      call chemical_rates(mla,n_old,rate1,dx,dt/2.d0)

      ! predictor
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt/2.d0,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

      !!!!!!!!!!!!!!!
      ! second step !
      !!!!!!!!!!!!!!!

      ! compute diffusive flux divergence
      call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

      ! calculate rates from 2*a(n_pred)-a(n_old)
      call chemical_rates(mla,n_old,rate2,dx,dt/2.d0,n_new,mattingly_lin_comb_coef)

      ! compute stochastic flux divergence
      if (variance_coef_mass .gt. 0.d0) then

        ! first, fill random flux multifabs with new random numbers
        call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

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
          ! compute 2*n_pred-n_old
          do n=1,nlevs
            call multifab_mult_mult_s_c(n_new(n),1,2.d0,nspecies,n_new(n)%ng)
            call multifab_sub_sub_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
          end do
          ! use 2*n_pred-n_old
          call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                    the_bc_tower,increment_in=.true.)
        case default
          call bl_error("advance_reaction_diffusion: invalid midpoint_stoch_flux_type")
        end select

      else

        do n=1,nlevs
          call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
        end do

      end if

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
        call multifab_saxpy_3(n_new(n),dt/2.d0,rate2(n))
        call multifab_saxpy_3(n_new(n),dt,ext_src(n))

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
      call multifab_destroy(ext_src(n))
      call multifab_destroy(rate1(n))

      if (temporal_integrator .eq. -2) then  ! explicit midpoint
        call multifab_destroy(rate2(n))
      end if
    end do

    call destroy(bpt)

  end subroutine advance_reaction_diffusion

end module advance_reaction_diffusion_module
