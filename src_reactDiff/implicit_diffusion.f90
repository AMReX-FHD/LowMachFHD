module implicit_diffusion_module

  use ml_layout_module
  use multifab_physbc_module
  use define_bc_module
  use bc_module
  use bndry_reg_module
  use ml_solve_module
  use probin_reactdiff_module, only: nspecies, diffusion_stencil_order, mg_verbose, &
                                     cg_verbose, &
                                     implicit_diffusion_rel_eps, implicit_diffusion_abs_eps

  implicit none

  private

  public :: implicit_diffusion

contains

  ! Note: diff_fluxdiv enters holding (div D_k grad n_k)^n
  ! Note: ext_src may contain stochastic contribution
  !
  ! This solves the following Crank-Nicolson system:
  ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
  !                   + (dt/2)(div D_k grad n_k)^n+1
  !                   +  dt    ext_src
  ! 
  ! using delta formulation, in operator form this becomes:
  !
  ! (I - div_(dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
  !                                       + dt ext_src      
  !
  subroutine implicit_diffusion(mla,n_old,n_new,ext_src,diff_coef_face,diff_fluxdiv, &
                                dx,dt,the_bc_tower)
    
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    type(multifab) , intent(in   ) :: ext_src(:)
    type(multifab) , intent(in   ) :: diff_coef_face(:,:)
    type(multifab) , intent(in   ) :: diff_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! for multigrid solver; (alpha - div beta grad) phi = rhs
    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: rhs(mla%nlevel)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)

    ! for diffusion multigrid - not used but needs to be passed in
    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    integer :: i,dm,n,nlevs,spec

    type(bl_prof_timer),save :: bpt

    call build(bpt,"implicit_diffusion")

    nlevs = mla%nlevel
    dm = mla%dim

    ! for multigrid solver; (alpha - div beta grad) phi = rhs
    do n=1,nlevs
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

    ! copy n_new = n_old
    do n=1,nlevs
       call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
    end do

    ! Crank-Nicolson diffusion solve

    ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
    !                   + (dt/2)(div D_k grad n_k)^n+1
    !                   +  dt    ext_src
    ! 
    ! using delta formulation, in operator form this becomes:
    !
    ! (I - div_(dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
    !                                       + dt ext_src      

    ! alpha=1 here for all components since it is simple difusion  
    do n=1,nlevs
       call multifab_setval(alpha(n),1.d0,all=.true.)
    end do

    do spec=1,nspecies

       ! beta = (dt/2)*D_k
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(beta(n,i),1,diff_coef_face(n,i),spec,1,0)
             call multifab_mult_mult_s(beta(n,i),0.5d0*dt)
          end do
       end do

       ! rhs =   dt div (D_k grad n_k^n)
       !       + dt ext_src    

       do n=1,nlevs
          call multifab_copy_c(rhs(n),1,diff_fluxdiv(n),spec,1,0)
          call multifab_plus_plus_c(rhs(n),1,ext_src(n),spec,1,0)
          call multifab_mult_mult_s(rhs(n),dt)
       end do

       ! initial guess for phi is zero since we are solving delta formulation
       do n=1,nlevs
          call multifab_setval(phi(n),0.d0,all=.true.)
       end do

       ! solve the implicit system
       call ml_cc_solve(mla,rhs,phi,fine_flx,alpha,beta,dx, &
                        the_bc_tower,scal_bc_comp+spec-1, &
                        stencil_order=diffusion_stencil_order, &
                        verbose=mg_verbose, &
                        cg_verbose=cg_verbose, &
                        eps=implicit_diffusion_rel_eps, &
                        abs_eps=implicit_diffusion_abs_eps)

       ! n_new = n_old + delta n
       do n=1,nlevs
          call multifab_plus_plus_c(n_new(n),spec,phi(n),1,1,0)
       end do

    end do

    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    do n=1,nlevs
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

    call destroy(bpt)

  end subroutine implicit_diffusion

end module implicit_diffusion_module
