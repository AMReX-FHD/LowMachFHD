module advance_timestep_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use bndry_reg_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use ml_solve_module  
  use probin_common_module, only: algorithm_type
  use probin_reactdiff_module, only: nspecies, mg_verbose, cg_verbose, D_Fick

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

    ! for multigrid solver; (alpha - div beta grad) phi = rhs
    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: rhs(mla%nlevel)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)

    ! for diffusion multigrid - not used but needs to be passed in
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    integer :: n,nlevs,i,dm,comp

    nlevs = mla%nlevel
    dm = mla%dim
    
    do n=1,nlevs
       call multifab_build(diff_fluxdiv(n) ,mla%la(n),nspecies,0) 
       call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0) 
       do i=1,dm
          call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
       end do
       ! for multigrid solver; (alpha - div beta grad) phi = rhs
       call multifab_build(alpha(n),mla%la(n),1,0)
       call multifab_build(rhs(n)  ,mla%la(n),1,0)
       call multifab_build(phi(n)  ,mla%la(n),1,1) 
       do i=1,dm
          call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! stores beta*grad phi/dx_fine on coarse-fine interfaces
    ! this gets computed inside of ml_cc_solve
    ! we pass it back out because some algorithms (like projection methods) 
    ! use this information
    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
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
       ! n_k^{n+1,*} = n_k^n + dt(div D_k grad n_k)^n
       !                     + dt(div (sqrt(2 D_k n_k / dt) Z)^n
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_fill_boundary(n_new(n))
       end do

       ! compute diffusive flux diverge
       call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

       ! Trapezoidal corrector:
       ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
       !                   + (dt/2)(div D_k grad n_k)^{n+1,*}
       !                   + (dt)(div (sqrt(2 D_k n_k / dt) Z)^n
       do n=1,nlevs
          call multifab_plus_plus_c(n_new(n),1,n_old(n),1,nspecies,0)
          call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
          call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
          call multifab_mult_mult_s_c(n_new(n),1,0.5d0,nspecies,0)
          call multifab_fill_boundary(n_new(n))
       end do

    else if (algorithm_type .eq. 1) then
       ! Crank-Nicolson

       ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
       !                   + (dt/2)(div D_k grad n_k)^n+1
       !                   +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       ! 
       ! in operator form
       !
       ! (I - (dt/2) div D_k grad)n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
       !                                            +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
       !

         
       ! alpha=1 here for all components since it is simple difusion  
       do n=1,nlevs
          call multifab_setval(alpha(n),1.d0,all=.true.)
       end do   
       do comp=1,nspecies

          ! beta = (dt/2)*D_k
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(beta(n,i),1,diff_coef_face(n,i),comp,1,0)
                call multifab_mult_mult_s(beta(n,i),0.5d0*dt)
             end do
          end do

          ! rhs = n_k^n + (dt/2)(div D_k grad n_k)^n
          !             +  dt    div (sqrt(2 D_k n_k / dt) Z)^n
          do n=1,nlevs
             call multifab_copy_c(rhs(n),1,diff_fluxdiv(n),comp,1,0)
             call multifab_mult_mult_s(rhs(n),0.5d0)
             call multifab_plus_plus_c(rhs(n),1,stoch_fluxdiv(n),comp,1,0)
             call multifab_mult_mult_s(rhs(n),dt)
             call multifab_plus_plus_c(rhs(n),1,n_old(n),comp,1,0)
          end do

          ! initial guess for phi is n_k^n
          do n=1,nlevs
             call multifab_copy_c(phi(n),1,n_old(n),comp,1,1)
          end do

          ! solve the implicit system
          call ml_cc_solve(mla,rhs,phi,fine_flx,alpha,beta,dx, &
               the_bc_tower,scal_bc_comp+comp-1, &
               verbose=mg_verbose, &
               cg_verbose=cg_verbose)

          ! copy solution into n_new
          do n=1,nlevs
             call multifab_copy_c(n_new(n),comp,phi(n),1,1,0)
          end do

       end do

       do n=1,nlevs
          call multifab_fill_boundary(n_new(n))
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
       call multifab_destroy(alpha(n))
       call multifab_destroy(rhs(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(beta(n,i))
       end do
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

  end subroutine advance_timestep

end module advance_timestep_module
