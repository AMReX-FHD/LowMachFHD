module advance_explicit_midpoint_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use diffusive_n_fluxdiv_module
  use stochastic_n_fluxdiv_module
  use simulate_reaction_module
  use probin_common_module, only: variance_coef_mass
  use probin_reactdiff_module, only: nspecies,D_Fick

  implicit none

  private

  public :: advance_explicit_midpoint

contains

  subroutine advance_explicit_midpoint(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src_in(:)

    ! local
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)
    type(multifab) :: ext_src(mla%nlevel)
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: n_react1(mla%nlevel)
    type(multifab) :: n_react2(mla%nlevel)

    integer :: nlevs,dm,n,i,spec

    type(bl_prof_timer), save :: bpt

    !!! there are three options for stochastic flux updating for the corrector 
    !!! stoch_flux_type = 1: K_n W1 + K_n W2
    !!! stoch_flux_type = 2: K_n W1 + K_pred W2
    !!! stoch_flux_type = 3: K_n W1 + K_{2*pred-n} W2 
    integer, parameter :: stoch_flux_type = 2 

    !!!!!!!!
    ! init !
    !!!!!!!!

    if ((multifab_volume(n_old(1))/nspecies)<=1) then
      call bl_error("one cell case has not been implemented.")
    end if

    call build(bpt,"advance_explicit_midpoint")

    nlevs = mla%nlevel
    dm = mla%dim

    ! build multifabs
    do n=1,nlevs
      do i=1,dm
        call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
      end do
      call multifab_build(ext_src(n),mla%la(n),nspecies,0)
      call multifab_build(diff_fluxdiv(n),mla%la(n),nspecies,0)
      call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0)
      call multifab_build(n_react1(n),mla%la(n),nspecies,0)
      call multifab_build(n_react2(n),mla%la(n),nspecies,0)
    end do

    ! set diffusion coefficients
    ! If one wants a space-dependent D or state=dependent D,
    ! see multispecies code as example
    do n=1,nlevs
      do i=1,dm
        do spec=1,nspecies
          call multifab_setval_c(diff_coef_face(n,i),D_Fick(spec),spec,1,all=.true.)
        end do
      end do
    end do

    ! set ext_src
    if (present(ext_src_in)) then
      do n=1,nlevs
        call multifab_copy_c(ext_src(n),1,ext_src_in(n),1,nspecies,0)
      end do
    else
      do n=1,nlevs
        call multifab_setval(ext_src(n),0.d0,all=.true.)
      end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! advance dt/2 by euler (predictor) !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! set stochastic flux divergence
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

    ! simulate reaction
    call simulate_reaction(mla,n_old,n_react1,dx,dt/2.d0,the_bc_tower)

    ! predictor
    do n=1,nlevs
      call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
      call multifab_saxpy_3(n_new(n),dt/2.d0,diff_fluxdiv(n))
      call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
      call multifab_plus_plus_c(n_new(n),1,n_react1(n),1,nspecies,0)
      call multifab_saxpy_3(n_new(n),dt/2.d0,ext_src(n))

      call multifab_fill_boundary(n_new(n))
      call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                           the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    !!!!!!!!!!!!!
    ! corrector !
    !!!!!!!!!!!!!

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! set stoch_flux here, if 
    ! variance_coef_mass > 0 and stoch_flux_type = 1 or 2,
    ! variance_coef_mass <= 0  

    if (variance_coef_mass .gt. 0.d0) then
        ! fill random flux multifabs with new random numbers
        call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
    else
      do n=1,nlevs
        call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
      end do
    end if

    if (stoch_flux_type .eq. 1 .and. variance_coef_mass .gt. 0.d0) then
      call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                the_bc_tower,increment_in=.true.)
    end if

    if (stoch_flux_type .eq. 2 .and. variance_coef_mass .gt. 0.d0) then
      call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                the_bc_tower,increment_in=.true.)
    end if

    ! compute 2*n_pred-n_old
    do n=1,nlevs
      call multifab_mult_mult_s_c(n_new(n),1,2.d0,nspecies,0)
      call multifab_sub_sub_c(n_new(n),1,n_old(n),1,nspecies,0)
    end do

    ! set stoch_flux here, if 
    ! variance_coef_mass > 0 and stoch_flux_type = 3
    if (stoch_flux_type .eq. 3 .and. variance_coef_mass .gt. 0.d0) then
      call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                the_bc_tower,increment_in=.true.)
    end if

    ! simulate reaction
    call simulate_reaction(mla,n_new,n_react2,dx,dt/2.d0,the_bc_tower)

    ! corrector
    do n=1,nlevs
      call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
      call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
      call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
      call multifab_plus_plus_c(n_new(n),1,n_react1(n),1,nspecies,0)
      call multifab_plus_plus_c(n_new(n),1,n_react2(n),1,nspecies,0)
      call multifab_saxpy_3(n_new(n),dt,ext_src(n))

      call multifab_fill_boundary(n_new(n))
      call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                           the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    !!!!!!!!!!!!!!!!!!!!!
    ! destroy multifabs !
    !!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
      do i=1,dm
        call multifab_destroy(diff_coef_face(n,i))
      end do
      call multifab_destroy(ext_src(n))
      call multifab_destroy(diff_fluxdiv(n))
      call multifab_destroy(stoch_fluxdiv(n))
      call multifab_destroy(n_react1(n))
      call multifab_destroy(n_react2(n))
    end do
   
    call destroy(bpt)

  end subroutine advance_explicit_midpoint

end module advance_explicit_midpoint_module
