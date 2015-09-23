module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use multifab_physbc_module
  use implicit_diffusion_module
  use probin_common_module, only: algorithm_type
  use probin_reactdiff_module, only: nspecies, D_Fick

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)

    integer :: n,nlevs,i,dm,comp

    nlevs = mla%nlevel
    dm = mla%dim
    
    do n=1,nlevs
       call multifab_build(diff_fluxdiv(n) ,mla%la(n),nspecies,0) 
       call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0) 
       do i=1,dm
          call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! fill random flux multifabs with new random numbers
    call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

    ! compute the diffusion coefficients (for now just setting each to a different constant)
    ! If one wants a space-dependent D or state-dependent D see multispecies code as example
    ! We have a routine average_cc_to_face there that is meant to compute face-averaged values
    do n=1,nlevs
       do i=1,dm
          do comp=1,nspecies
             call multifab_setval_c(diff_coef_face(n,i), D_Fick(comp),comp,1,all=.true.)
          end do
       end do
    end do

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! compute stochastic flux divergence
    call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt,the_bc_tower)

    if (algorithm_type .eq. 0) then
       ! explicit predictor-corrector

       ! Euler predictor
       ! n_k^{n+1,*} = n_k^n + dt div (D_k grad n_k)^n
       !                     + dt div (sqrt(2 D_k n_k / dt) Z)^n
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute diffusive flux divergence
       call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

       ! Trapezoidal corrector:
       ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
       !                   + (dt/2) div (D_k grad n_k)^{n+1,*}
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       do n=1,nlevs
          call multifab_plus_plus_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_mult_mult_s_c(n_new(n),1,0.5d0,nspecies,0)
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

    else if (algorithm_type .eq. 1) then
       ! Crank-Nicolson
       ! n_k^{n+1} = n_k^n + (dt/2) div )D_k grad n_k)^n
       !                   + (dt/2) div (D_k grad n_k)^n+1
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       call implicit_diffusion(mla,n_old,n_new,diff_coef_face,diff_fluxdiv,stoch_fluxdiv, &
                               dx,dt,the_bc_tower)

    else if (algorithm_type .eq. 2) then
       ! explicit midpoint scheme

       ! n_k^{n+1/2} = n_k^n + (dt/2) div (D_k grad n_k)^n
       !                     +  dt    div (sqrt(D_k n_k / dt) Z_1)^n
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt/2.d0      ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute diffusive flux divergence at t^{n+1/2}
       call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

       ! fill random flux multifabs with new random numbers
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

       ! compute second-stage stochastic flux divergence
       call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                 the_bc_tower,increment_in=.true.)

       ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
       !                   + dt div (sqrt(D_k n_k / dt) Z_1)^n
       !                   + dt div (sqrt(D_k n_k / dt) Z_2)^n
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt           ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

    else

       call bl_error("invalid algorithm_type")
    end if

    do n=1,nlevs
       call multifab_destroy(diff_fluxdiv(n))
       call multifab_destroy(stoch_fluxdiv(n))
       do i=1,dm
          call multifab_destroy(diff_coef_face(n,i))
       end do
    end do

  end subroutine advance_timestep

end module advance_timestep_module
