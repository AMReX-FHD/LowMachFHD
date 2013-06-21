module macproject_module

  use ml_layout_module
  use define_bc_module
  use bl_constants_module
  use bl_error_module
  use bc_module 
  use convert_stag_module
  use mg_module             , only: mg_tower, mg_tower_build, mg_tower_destroy
  use cc_stencil_fill_module, only: stencil_fill_cc_all_mglevels
  use probin_gmres_module   , only: mg_max_vcycles, mg_rel_tol, cg_verbose, &
                                    mg_bottom_solver, mg_max_bottom_nlevels, &
                                    mg_minwidth, mg_nsmooths_bottom, &
                                    mg_nsmooths_down, mg_nsmooths_up, mg_verbose
  use vcycle_counter_module
  use stencil_types_module
  use ml_solve_module       , only : ml_cc_solve
  use bc_module
  use multifab_physbc_module
  use div_and_grad_module

  implicit none

  private

  public :: macproject, subtract_weighted_gradp

contains 

  ! solve L_alpha Phi  = D x_u^* - b_p
  ! does not update any other variables
  subroutine macproject(mla,phi,umac,alphainv_fc,mac_rhs,dx,the_bc_tower,full_solve_in)

    use bndry_reg_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: phi(:)
    type(multifab ), intent(in   ) :: umac(:,:)
    type(multifab ), intent(inout) :: alphainv_fc(:,:)
    type(multifab ), intent(inout) :: mac_rhs(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    logical, intent(in), optional  :: full_solve_in

    ! Local  
    type(multifab)  :: zero_fab(mla%nlevel)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    real(kind=dp_t) :: umac_norm(mla%nlevel)
    real(kind=dp_t) :: rel_solver_eps
    real(kind=dp_t) :: abs_solver_eps
    integer         :: dm,i,n,nlevs
    logical         :: full_solve

    if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
       print*,"Begin call to macproject"
    end if

    full_solve = .false.
    if (present(full_solve_in)) then
       full_solve = full_solve_in
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    do n = 1, nlevs
       call multifab_build(zero_fab(n), mla%la(n),  1, 0)
       call setval(zero_fab(n),ZERO,all=.true.)
       call setval(phi(n)     ,ZERO,all=.true.)
    end do

    ! Compute umac_norm to be used inside the MG solver as part of a stopping criterion
    umac_norm = -1.0_dp_t
    do n = 1,nlevs
       do i = 1,dm
          umac_norm(n) = max(umac_norm(n),norm_inf(umac(n,i)))
       end do
    end do

    ! multiply alphainv_fc by -1 so we solve L_alpha Phi = mac_rhs
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(alphainv_fc(n,i),1,-1.d0,1,alphainv_fc(n,i)%ng)
       end do
    end do

    ! stores (alphainv_fc/dx**2) grad phi at coarse-fine interfaces
    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    rel_solver_eps = mg_rel_tol
    abs_solver_eps = 1.d-16

    if (full_solve) then
       call mac_multigrid(mla,mac_rhs,phi,fine_flx,zero_fab,alphainv_fc,dx,the_bc_tower,pres_bc_comp, &
                          mla%mba%rr,rel_solver_eps,abs_solver_eps, &
                          abort_on_max_iter=.true.)
    else
       call mac_multigrid(mla,mac_rhs,phi,fine_flx,zero_fab,alphainv_fc,dx,the_bc_tower,pres_bc_comp, &
                          mla%mba%rr,rel_solver_eps,abs_solver_eps, &
                          abort_on_max_iter=.false.)
    end if

    ! restore alphainv_fc
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(alphainv_fc(n,i),1,-1.d0,1,alphainv_fc(n,i)%ng)
       end do
    end do

    vcycle_counter = vcycle_counter + mg_max_vcycles

    do n = 1, nlevs
       call multifab_destroy(zero_fab(n))
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
       print*,"End call to macproject"
    end if

  contains

    subroutine mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,pres_bc_comp, &
                             ref_ratio,rel_solver_eps,abs_solver_eps, &
                             abort_on_max_iter)

      type(ml_layout), intent(in   ) :: mla
      type(multifab) , intent(inout) :: rh(:),phi(:)
      type(bndry_reg), intent(inout) :: fine_flx(2:)
      type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
      real(dp_t)     , intent(in   ) :: dx(:,:)
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: pres_bc_comp
      integer        , intent(in   ) :: ref_ratio(:,:)
      real(dp_t)     , intent(in   ) :: rel_solver_eps
      real(dp_t)     , intent(in   ) :: abs_solver_eps
      logical        , intent(in   ) :: abort_on_max_iter

      type(layout  ) :: la
      type(box     ) :: pd

      type(multifab), allocatable :: cell_coeffs(:)
      type(multifab), allocatable :: face_coeffs(:,:)

      type(mg_tower)  :: mgt(mla%nlevel)
      integer         :: dm, ns, nlevs

      ! MG solver defaults
      integer    :: stencil_type, bottom_max_iter, max_nlevel
      integer    :: d, n, nub, gamma, cycle_type, smoother
      integer    :: max_nlevel_in,do_diagnostics
      real(dp_t) :: omega,bottom_solver_eps
      real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
      real(dp_t) :: pxa(mla%dim), pxb(mla%dim)

      type(bl_prof_timer), save :: bpt

      call build(bpt, "mac_multigrid")

      nlevs = mla%nlevel
      dm    = mla%dim

      !! Defaults:

      max_nlevel        = mgt(nlevs)%max_nlevel
      smoother          = mgt(nlevs)%smoother
      nub               = mgt(nlevs)%nub
      gamma             = mgt(nlevs)%gamma
      omega             = mgt(nlevs)%omega
      cycle_type        = mgt(nlevs)%cycle_type
      bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
      bottom_max_iter   = mgt(nlevs)%bottom_max_iter

      ns = 1 + dm*3

      do n = nlevs, 1, -1

         if (n == 1) then
            max_nlevel_in = max_nlevel
         else
            if ( all(ref_ratio(n-1,:) == 2) ) then
               max_nlevel_in = 1
            else if ( all(ref_ratio(n-1,:) == 4) ) then
               max_nlevel_in = 2
            else
               call bl_error("MAC_MULTIGRID: confused about ref_ratio")
            end if
         end if

         pd = layout_get_pd(mla%la(n))

         stencil_type = CC_CROSS_STENCIL

         if (full_solve) then

            call mg_tower_build(mgt(n), mla%la(n), pd, &
                                the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,pres_bc_comp),&
                                stencil_type, &
                                dh = dx(n,:), &
                                ns = ns, &
                                smoother = smoother, &
                                nu1 = mgt(nlevs)%nu1, &
                                nu2 = mgt(nlevs)%nu2, &
                                nuf = mgt(nlevs)%nuf, &
                                nub = nub, &
                                gamma = gamma, &
                                cycle_type = cycle_type, &
                                omega = omega, &
                                bottom_solver = mgt(nlevs)%bottom_solver, &
                                bottom_max_iter = bottom_max_iter, &
                                bottom_solver_eps = bottom_solver_eps, &
                                max_iter = mgt(nlevs)%max_iter, &
                                abort_on_max_iter = abort_on_max_iter, &
                                max_nlevel = max_nlevel_in, &
                                max_bottom_nlevel = mgt(nlevs)%max_bottom_nlevel, &
                                min_width = mgt(nlevs)%min_width, &
                                eps = rel_solver_eps, &
                                abs_eps = abs_solver_eps, &
                                verbose = mg_verbose, &
                                cg_verbose = cg_verbose, &
                                nodal = nodal_flags(rh(nlevs)))

         else

            call mg_tower_build(mgt(n), mla%la(n), pd, &
                                the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,pres_bc_comp),&
                                stencil_type, &
                                dh = dx(n,:), &
                                ns = ns, &
                                smoother = smoother, &
                                nu1 = mg_nsmooths_down, &
                                nu2 = mg_nsmooths_up, &
                                nuf = mg_nsmooths_bottom, &
                                nub = nub, &
                                gamma = gamma, &
                                cycle_type = cycle_type, &
                                omega = omega, &
                                bottom_solver = mg_bottom_solver, &
                                bottom_max_iter = bottom_max_iter, &
                                bottom_solver_eps = bottom_solver_eps, &
                                max_iter = mg_max_vcycles, &
                                abort_on_max_iter = abort_on_max_iter, &
                                max_nlevel = max_nlevel_in, &
                                max_bottom_nlevel = mg_max_bottom_nlevels, &
                                min_width = mg_minwidth, &
                                eps = rel_solver_eps, &
                                abs_eps = abs_solver_eps, &
                                verbose = mg_verbose, &
                                cg_verbose = cg_verbose, &
                                nodal = nodal_flags(rh(nlevs)))

            end if

      end do

      !! Fill coefficient array

      do n = nlevs,1,-1

         allocate(cell_coeffs(mgt(n)%nlevels))
         allocate(face_coeffs(mgt(n)%nlevels,dm))

         la = mla%la(n)

         call multifab_build(cell_coeffs(mgt(n)%nlevels), la, 1, 1)
         call multifab_copy_c(cell_coeffs(mgt(n)%nlevels),1,alpha(n),1, 1,ng=nghost(alpha(n)))

         do d = 1, dm
            call multifab_build_edge(face_coeffs(mgt(n)%nlevels,d),la,1,1,d)
            call multifab_copy_c(face_coeffs(mgt(n)%nlevels,d),1,beta(n,d),1,1,ng=nghost(beta(n,d)))
         end do

         if (n > 1) then
            xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
            xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
         else
            xa = ZERO
            xb = ZERO
         end if

         pxa = ZERO
         pxb = ZERO

         call stencil_fill_cc_all_mglevels(mgt(n), cell_coeffs, face_coeffs, &
              xa, xb, pxa, pxb, 2, &
              the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,pres_bc_comp))

         call destroy(cell_coeffs(mgt(n)%nlevels))
         deallocate(cell_coeffs)

         do d = 1, dm
            call destroy(face_coeffs(mgt(n)%nlevels,d))
         end do
         deallocate(face_coeffs)

      end do

      if (mg_verbose >= 3) then
         do_diagnostics = 1
      else
         do_diagnostics = 0
      end if

      call ml_cc_solve(mla, mgt, rh, phi, fine_flx, ref_ratio,do_diagnostics, &
           rel_solver_eps, abs_solver_eps)

      do n = 1,nlevs
         call multifab_fill_boundary(phi(n))
         call multifab_physbc(phi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:))
      end do

      do n = 1, nlevs
         call mg_tower_destroy(mgt(n))
      end do

    end subroutine mac_multigrid

  end subroutine macproject

  subroutine subtract_weighted_gradp(mla,x_u,alphainv_fc,phi,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: x_u(:,:)
    type(multifab ), intent(in   ) :: alphainv_fc(:,:)
    type(multifab ), intent(in   ) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,dm,n,nlevs

    type(multifab) :: gradp(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(gradp(n,i),mla%la(n),1,1,i)
       end do
    end do

    call compute_grad(mla,phi,gradp,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(gradp(n,i),1,alphainv_fc(n,i),1,1,0)
          call saxpy(x_u(n,i),-1.d0,gradp(n,i))
       end do
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(gradp(n,i))
       end do
    end do

  end subroutine subtract_weighted_gradp
  
end module macproject_module
