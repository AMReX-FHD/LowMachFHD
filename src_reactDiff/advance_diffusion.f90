module advance_diffusion_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use multifab_physbc_module
  use implicit_diffusion_module
  use probin_common_module, only: variance_coef_mass
  use probin_reactdiff_module, only: nspecies, D_Fick, diffusion_type, midpoint_stoch_flux_type

  implicit none

  private

  public :: advance_diffusion

contains

  ! Solves n_t = div ( D grad (n)) + div (sqrt(2*variance*D*n)*W) + g
  ! where g is a constant in time external source
  subroutine advance_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src_in(:) ! The value of g

    ! local
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)
    type(multifab) :: ext_src(mla%nlevel)

    integer :: n,nlevs,i,dm,spec

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim

    ! do not do diffusion if only one cell (well-mixed system)
    if((multifab_volume(n_old(1))/nspecies)<=1) then
       do n=1,nlevs
          ! make sure n_new contains the new state
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
       if (present(ext_src_in)) then
          call multifab_saxpy_3(n_new(n),dt,ext_src_in(n))
       end if
       return
    end if   
    
    call build(bpt,"advance_diffusion")
    
    do n=1,nlevs
       call multifab_build(diff_fluxdiv(n) ,mla%la(n),nspecies,0)
       call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0)
       do i=1,dm
          call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
       end do
       call multifab_build(ext_src(n),mla%la(n),nspecies,0)
    end do

    if (present(ext_src_in)) then
       do n=1,nlevs
          call multifab_copy_c(ext_src(n),1,ext_src_in(n),1,nspecies,0)
       end do
    else
       do n=1,nlevs
          call multifab_setval(ext_src(n),0.d0,all=.true.)
       end do
    end if

    ! fill random flux multifabs with new random numbers
    if (variance_coef_mass .gt. 0.d0) then
       call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
    end if

    ! compute the diffusion coefficients (for now just setting each to a different constant)
    ! If one wants a space-dependent D or state-dependent D see multispecies code as example
    ! We have a routine average_cc_to_face there that is meant to compute face-averaged values
    do n=1,nlevs
       do i=1,dm
          do spec=1,nspecies
             call multifab_setval_c(diff_coef_face(n,i), D_Fick(spec),spec,1,all=.true.)
          end do
       end do
    end do

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! compute stochastic flux divergence
    if (variance_coef_mass .gt. 0.d0) then
       call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                 the_bc_tower,increment_in=.false.)
    else
       do n=1,nlevs
          call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
       end do
    end if

    if (diffusion_type .eq. 0) then
       ! explicit predictor-corrector

       ! Euler predictor
       ! n_k^{n+1,*} = n_k^n + dt div (D_k grad n_k)^n
       !                     + dt div (sqrt(2 D_k n_k / dt) Z)^n
       !                     + dt ext_src
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,ext_src(n))
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
       !                   +  dt    ext_src
       do n=1,nlevs
          call multifab_plus_plus_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,ext_src(n))
          call multifab_mult_mult_s_c(n_new(n),1,0.5d0,nspecies,0)
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

    else if (diffusion_type .eq. 1) then
       ! Crank-Nicolson
       ! n_k^{n+1} = n_k^n + (dt/2) div )D_k grad n_k)^n
       !                   + (dt/2) div (D_k grad n_k)^n+1
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       !                   +  dt    ext_src
       call implicit_diffusion(mla,n_old,n_new,ext_src,diff_coef_face, &
                               diff_fluxdiv,stoch_fluxdiv,dx,dt,the_bc_tower)

    else if (diffusion_type .eq. 2) then
       ! explicit midpoint scheme

       ! n_k^{n+1/2} = n_k^n + (dt/2) div (D_k grad n_k)^n
       !                     + (dt/2) div (sqrt(2 D_k n_k / (dt/2) ) Z_1)^n
       !                     + (dt/2) ext_src
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt/2.d0      ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/2.d0      ,ext_src(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! compute diffusive flux divergence at t^{n+1/2}
       call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

       if (variance_coef_mass .gt. 0.d0) then
          ! fill random flux multifabs with new random numbers
          call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

          ! compute second-stage stochastic flux divergence and
          ! add to first-stage stochastic flux divergence
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
             call bl_error("advance_diffusion: invalid midpoint_stoch_flux_type")
          end select
       end if
       ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
       !                   + dt div (sqrt(2 D_k n_k^n dt) Z_1 / sqrt(2) )
       !                   + dt div (sqrt(2 D_k n_k^? dt) Z_2 / sqrt(2) )
       !                   + dt ext_src
       ! where
       ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
       !       = n_k^pred            (midpoint_stoch_flux_type=2)
       !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt           ,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt           ,ext_src(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

    else

       call bl_error("invalid diffusion_type")
    end if

    do n=1,nlevs
       call multifab_destroy(diff_fluxdiv(n))
       call multifab_destroy(stoch_fluxdiv(n))
       do i=1,dm
          call multifab_destroy(diff_coef_face(n,i))
       end do
       call multifab_destroy(ext_src(n))
    end do

    call destroy(bpt)

  end subroutine advance_diffusion

end module advance_diffusion_module
