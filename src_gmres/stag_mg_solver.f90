module stag_mg_solver_module

  use ml_layout_module
  use stag_applyop_module
  use convert_stag_module
  use define_bc_module
  use bc_module
  use multifab_physbc_stag_module
  use probin_gmres_module , only: stag_mg_omega, stag_mg_max_vcycles, &
                                  stag_mg_nsmooths_bottom, stag_mg_nsmooths_down, &
                                  stag_mg_nsmooths_up, stag_mg_rel_tol, &
                                  stag_mg_smoother, stag_mg_verbosity, &
                                  stag_mg_minwidth, stag_mg_bottom_solver, &
                                  stag_mg_max_bottom_nlevels
  use probin_common_module, only: visc_type, n_cells, max_grid_size
  use vcycle_counter_module

  implicit none

  private

  public :: stag_mg_solver

contains

  ! solve "(alpha*I - L) phi = rhs" using multigrid with Jacobi relaxation
  ! if abs(visc_type) = 1, L = div beta grad
  ! if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
  ! if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
  ! if visc_type > 1 we assume constant coefficients
  ! if visc_type < 1 we assume variable coefficients
  ! beta_cc, and gamma_cc are cell-centered
  ! alpha_fc, phi_fc, and rhs_fc are face-centered
  ! beta_ed is nodal (2d) or edge-centered (3d)
  ! phi_fc must come in initialized to some value, preferably a reasonable guess
  recursive subroutine stag_mg_solver(mla,alpha_fc,beta_cc,beta_ed,gamma_cc,theta, &
                                      phi_fc,rhs_fc,dx,the_bc_tower,do_fancy_bottom_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: alpha_fc(:,:) ! face-centered
    type(multifab) , intent(in   ) ::  beta_cc(:)   ! cell-centered
    type(multifab) , intent(in   ) ::  beta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma_cc(:)   ! cell-centered
    type(multifab) , intent(inout) ::   phi_fc(:,:) ! face-centered
    type(multifab) , intent(in   ) ::   rhs_fc(:,:) ! face-centered
    real(kind=dp_t), intent(in   ) :: theta,dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    logical        , intent(in   ), optional :: do_fancy_bottom_in

    ! local variables
    integer :: dm, i, j, vcycle, m, n, nlevs_mg
    integer :: color,color_start,color_end,nsmooths
    integer :: vcycle_counter_temp

    ! stores initial residual and current residual
    real(kind=dp_t) :: resid0(mla%dim),resid(mla%dim),resid_temp,resid0_l2(mla%dim),resid_l2(mla%dim),temp

    ! hold the problem domain and boxarray at level 1 as well as the current multigrid level
    type(box)      :: pd_base,pd
    type(boxarray) :: ba_base,ba

    ! cell-centered multifabs
    type(multifab), allocatable :: beta_cc_mg(:), gamma_cc_mg(:)

    ! face-centered multifabs
    type(multifab), allocatable :: alpha_fc_mg(:,:), rhs_fc_mg(:,:), phi_fc_mg(:,:)
    type(multifab), allocatable :: Lphi_fc_mg(:,:), resid_fc_mg(:,:)

    ! edge multifab
    type(multifab), allocatable :: beta_ed_mg(:,:)

    ! the layout at each level of multigrid
    type(layout)   , allocatable :: la_mg(:)

    ! grid spacing at each level of multigrid
    real(kind=dp_t), allocatable :: dx_mg(:,:)

    logical :: nodal_temp(mla%dim)

    type(bc_tower) :: the_bc_tower_mg

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! fancy bottom solver stuff
    logical :: do_fancy_bottom
    type(ml_boxarray) :: mba_fancy
    type(ml_layout) :: mla_fancy
    type(multifab) :: alpha_fc_fancy(1,mla%dim)
    type(multifab) :: beta_cc_fancy(1)
    type(multifab), allocatable :: beta_ed_fancy(:,:)
    type(multifab) :: gamma_cc_fancy(1)
    type(multifab) :: phi_fc_fancy(1,mla%dim)
    type(multifab) :: rhs_fc_fancy(1,mla%dim)
    real(kind=dp_t) :: dx_fancy(1,mla%dim)
    type(bc_tower) :: the_bc_tower_fancy
    integer :: lo(mla%dim), hi(mla%dim), bottom_box_size
    type(box) :: bx

    if (present(do_fancy_bottom_in)) then
       do_fancy_bottom = do_fancy_bottom_in
    else
       do_fancy_bottom = (stag_mg_bottom_solver .eq. 4)
    end if

    if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
       print*,""
       print*,"Begin call to stag_mg_solver"
    end if

    ! only works for single level
    if (mla%nlevel .ne. 1) then
       call bl_error("stag_mg_solver only works with mla%nlevel=1")
    end if

    dm = mla%dim
    
    ! get the problem domain from level 1
    pd_base = ml_layout_get_pd(mla,1)

    ! get boxarray from level 1
    ba_base = get_boxarray(mla%la(1))

    ! compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    ! smallest dimension of the smallest grid at the coarsest multigrid level
    call compute_nlevs_mg(nlevs_mg,ba_base)

    ! allocate multifabs used in multigrid coarsening
    allocate( beta_cc_mg(nlevs_mg))    ! cell-centered
    allocate(gamma_cc_mg(nlevs_mg))    ! cell-centered
    allocate(alpha_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(  rhs_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(  phi_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate( Lphi_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(resid_fc_mg(nlevs_mg,dm)) ! face-centered

    if (dm .eq. 2) then
       allocate( beta_ed_mg(nlevs_mg,1))  ! nodal
    else if (dm .eq. 3) then
       allocate( beta_ed_mg(nlevs_mg,3))  ! edge-based
    end if

    allocate(dx_mg(nlevs_mg,dm))
    allocate(la_mg(nlevs_mg))

    call bc_tower_init(the_bc_tower_mg,nlevs_mg,dm,the_bc_tower%domain_bc)

    do n=1,nlevs_mg

       ! compute dx at this level of multigrid
       dx_mg(n,:) = dx(1,:) * 2**(n-1)

       ! create the problem domain for this multigrid level
       pd = coarsen(pd_base,2**(n-1))

       ! create the boxarray for this multigrid level
       call boxarray_build_copy(ba,ba_base)
       call boxarray_coarsen(ba,2**(n-1))

       ! sanity check to make sure level 1 boxarrays match
       if (n .eq. 1 .and. (.not. boxarray_same_q(mla%mba%bas(1),ba) ) ) then
          call print(ba)
          call print(mla%mba%bas(1))
          call bl_error("Finest multigrid level boxarray and coarsest problem boxarrays do not match")
       end if

       ! build the layout, la
       ! force the same processor assignments as mla%la(1).  We can do this since there
       ! are the same number of boxes in the same order in physical space
       if (n .eq. 1) then
          la_mg(1) = mla%la(1)
       else
          call layout_build_ba(la_mg(n),ba,pd,mla%pmask,explicit_mapping=get_proc(mla%la(1)))
       end if

       ! don't need this anymore - free up memory
       call destroy(ba)

       ! build multifabs used in multigrid coarsening
       call multifab_build( beta_cc_mg(n),la_mg(n),1,1)
       call multifab_build(gamma_cc_mg(n),la_mg(n),1,1)

       do i=1,dm
          call multifab_build_edge(alpha_fc_mg(n,i),la_mg(n),1,0,i)
          call multifab_build_edge(  rhs_fc_mg(n,i),la_mg(n),1,0,i)
          call multifab_build_edge(  phi_fc_mg(n,i),la_mg(n),1,1,i)
          call multifab_build_edge( Lphi_fc_mg(n,i),la_mg(n),1,1,i)
          call multifab_build_edge(resid_fc_mg(n,i),la_mg(n),1,0,i)
       end do

       ! build beta_ed_mg
       if (dm .eq. 2) then
          call multifab_build_nodal(beta_ed_mg(n,1),la_mg(n),1,0)
       else
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(beta_ed_mg(n,1),la_mg(n),1,0,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(beta_ed_mg(n,2),la_mg(n),1,0,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(beta_ed_mg(n,3),la_mg(n),1,0,nodal_temp)
       end if

       call bc_tower_level_build(the_bc_tower_mg,n,la_mg(n))

    end do

    ! copy level 1 coefficients into mg array of coefficients
    call multifab_copy_c(beta_cc_mg(1),1,beta_cc(1),1,1,1)
    call multifab_copy_c(beta_ed_mg(1,1),1,beta_ed(1,1),1,1,0)
    if (dm .eq. 3) then
       call multifab_copy_c(beta_ed_mg(1,2),1,beta_ed(1,2),1,1,0)
       call multifab_copy_c(beta_ed_mg(1,3),1,beta_ed(1,3),1,1,0)
    end if
    call multifab_copy_c(gamma_cc_mg(1),1,gamma_cc(1),1,1,1)
    do i=1,dm
       call multifab_copy_c(alpha_fc_mg(1,i),1,alpha_fc(1,i),1,1,0)
       ! multiply alpha_fc_mg by theta
       call multifab_mult_mult_s_c(alpha_fc_mg(1,i),1,theta,1,0)
    end do

    ! coarsen coefficients
    do n=2,nlevs_mg
       call cc_restriction(la_mg(n), beta_cc_mg(n), beta_cc_mg(n-1),the_bc_tower_mg%bc_tower_array(n-1))
       call cc_restriction(la_mg(n),gamma_cc_mg(n),gamma_cc_mg(n-1),the_bc_tower_mg%bc_tower_array(n-1))
       call stag_restriction(la_mg(n),alpha_fc_mg(n-1,:),alpha_fc_mg(n,:),.true.)
       if (dm .eq. 2) then
          call nodal_restriction(la_mg(n),beta_ed_mg(n-1,1),beta_ed_mg(n,1))
       else
          call edge_restriction(la_mg(n),beta_ed_mg(n-1,:),beta_ed_mg(n,:))
       end if
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we wolve the homogeneous problem
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,dm

       ! initialize phi_fc_mg = phi_fc as an initial guess
       call multifab_copy_c(phi_fc_mg(1,i),1,phi_fc(1,i),1,1,0)

       ! set values on physical boundaries
       call multifab_physbc_domainvel(phi_fc_mg(1,i),vel_bc_comp+i-1, &
                                      the_bc_tower_mg%bc_tower_array(1),dx_mg(1,:))
       
       ! fill periodic ghost cells
       call multifab_fill_boundary(phi_fc_mg(1,i))

       ! fill physical ghost cells
       call multifab_physbc_macvel(phi_fc_mg(1,i),vel_bc_comp+i-1, &
                                   the_bc_tower_mg%bc_tower_array(1),dx_mg(1,:))

    end do

    ! set rhs_fc_mg at level 1 by copying in passed-in rhs_fc
    do i=1,dm
       call multifab_copy_c(rhs_fc_mg(1,i),1,rhs_fc(1,i),1,1,0)
    end do

    ! compute norm of initial residual
    ! first compute Lphi
    call stag_applyop_level(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                            Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                            beta_cc_mg(1),beta_ed_mg(1,:),gamma_cc_mg(1),dx_mg(1,:))

    do j=1,dm
       ! compute Lphi - rhs
       call multifab_sub_sub_c(Lphi_fc_mg(1,j),1,rhs_fc_mg(1,j),1,1,0)
       ! compute L0 norm of Lphi - rhs
       resid0(j) = multifab_norm_inf_c(Lphi_fc_mg(1,j),1,1)
       resid0_l2(j) = multifab_norm_l2_c(Lphi_fc_mg(1,j),1,1)
       if (parallel_IOProcessor()  .and. stag_mg_verbosity .ge. 2) then
          print*,"Initial residual",j,resid0(j)
       end if
    end do

    if (all(resid0(1:dm) .eq. 0.d0)) then
       if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
          print*,"Initial residual is zero; exiting staggered multigrid solver"
          return
       end if
    end if

    ! if some (but not all) of the residuals are zero
    ! set the zero residuals to the maximum so the multigrid will begin work
    if (any(resid0(1:dm) .eq. 0.d0)) then
       resid0(1:dm)    = maxval(resid0(1:dm))
       resid0_l2(1:dm) = maxval(resid0_l2(1:dm))
    end if

    if (stag_mg_smoother .eq. 0) then
       color_start = 0
       color_end = 0
    else
       color_start = 1
       color_end = 2*dm
    end if

    do vcycle=1,stag_mg_max_vcycles

       if (parallel_IOProcessor()  .and. stag_mg_verbosity .ge. 2) then
          print*,""
          print*,"Begin V-Cycle",vcycle
       end if

       ! set phi to zero at coarser levels as initial guess for residual equation
       do j=1,dm
          do n=2,nlevs_mg
             call setval(phi_fc_mg(n,j), 0.d0, all=.true.)
          end do
       end do

       ! down the V-cycle
       do n=1,nlevs_mg-1

          ! print out residual
          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

             do j=1,dm
                ! compute Lphi - rhs, and report residual
                call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"before    smooths at level",n,resid_temp
                end if
                
             end do
             
          end if
             
          ! control to do a different number of smooths for the bottom solver
          nsmooths = stag_mg_nsmooths_down

          do m=1,nsmooths

             ! do the smooths
             do color=color_start,color_end

                ! the form of weighted Jacobi we are using is
                ! phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                ! where D is the diagonal matrix containing the diagonal elements of L

                ! compute Lphi
                call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_ed_mg(n,:), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)

                ! update phi = phi + omega*D^{-1}*(rhs-Lphi)
                call stag_mg_update(la_mg(n),phi_fc_mg(n,:),rhs_fc_mg(n,:), &
                                    Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                    beta_cc_mg(n),beta_ed_mg(n,:), &
                                    gamma_cc_mg(n),dx_mg(n,:),color)

                do j=1,dm
                   
                   ! set values on physical boundaries
                   call multifab_physbc_domainvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                                  the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))
                   
                   ! fill periodic ghost cells
                   call multifab_fill_boundary(phi_fc_mg(n,j))

                   ! fill physical ghost cells
                   call multifab_physbc_macvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                               the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

                end do

             end do ! end loop over colors

          end do ! end loop over nsmooths

          !!!!!!!!!!!!!!!!!!!
          ! compute residual

          ! compute Lphi
          call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                  Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                  beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

          do j=1,dm

             ! compute Lphi - rhs, and then multiply by -1
             call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
             call multifab_mult_mult_s_c(Lphi_fc_mg(n,j),1,-1.d0,1,0)
             if (stag_mg_verbosity .ge. 3) then
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"after all smooths at level",n,resid_temp
                end if
             end if

             ! set values on physical boundaries
             call multifab_physbc_domainvel(Lphi_fc_mg(n,j),vel_bc_comp+j-1, &
                                            the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

             ! fill periodic ghost cells
             call multifab_fill_boundary(Lphi_fc_mg(n,j))

             ! fill physical ghost cells
             call multifab_physbc_macvel(Lphi_fc_mg(n,j),vel_bc_comp+j-1, &
                                         the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

          end do

          ! restrict/coarsen residual and put it in rhs_fc
          call stag_restriction(la_mg(n),Lphi_fc_mg(n,:),rhs_fc_mg(n+1,:))

          do j=1,dm
             ! set residual to zero on physical boundaries
             call multifab_physbc_domainvel(rhs_fc_mg(n+1,j),vel_bc_comp+j-1, &
                                            the_bc_tower_mg%bc_tower_array(n+1),dx_mg(n+1,:))
          end do

       end do ! end loop over nlevs_mg

       ! bottom solve
       n = nlevs_mg

       if (do_fancy_bottom .and. all(n_cells(1:dm) / max_grid_size(1:dm) .ge. 2) ) then
          
          !!!!!!!!!!!!!!!!!!!!!!!
          ! fancy bottom solver
          ! puts the entire bottom solve on a single grid and continues to coarsen

          ! tell mba how many levels and dmensionality of problem
          call ml_boxarray_build_n(mba_fancy,1,mla%dim)

          ! set dx equal to dx of the bottom solve
          dx_fancy(1,:) = dx_mg(nlevs_mg,:)
          
          ! create a box containing number of cells of current bottom solve
          lo(1:dm) = 0
          hi(1:dm) = stag_mg_minwidth * (n_cells(1:dm) / max_grid_size(1:dm)) - 1
          bx = make_box(lo,hi)

          ! tell mba about the problem domain
          mba_fancy%pd(1) = bx

          ! initialize the boxarray at level 1 to be one single box
          call boxarray_build_bx(mba_fancy%bas(1),bx)

          ! chop up the box to respect stag_mg_max_bottom_nlevels
          bottom_box_size = maxval(hi(1:dm))+1
          bottom_box_size = min(bottom_box_size,stag_mg_minwidth*2**(stag_mg_max_bottom_nlevels-1))
          call boxarray_maxsize(mba_fancy%bas(1),bottom_box_size)

          ! build the ml_layout, mla_fancy
          call ml_layout_build(mla_fancy,mba_fancy,mla%pmask)

          if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
             print*,'Invoking Fancy Bottom Solver'
             print*,'# of grids of old bottom solve:',nboxes(mla%la(1))
             print*,'# of grids of new bottom solve:',nboxes(mla_fancy%la(1))
          end if

          ! don't need this anymore - free up memory
          call destroy(mba_fancy)

          ! tell the_bc_tower about max_levs, dm, and domain_phys_bc
          call initialize_bc(the_bc_tower_fancy,1,dm,mla_fancy%pmask,2,1)

          ! define level 1 of the_bc_tower
          call bc_tower_level_build(the_bc_tower_fancy,1,mla_fancy%la(1))
          
          call multifab_build(beta_cc_fancy(1),mla_fancy%la(1),1,1)
          call multifab_build(gamma_cc_fancy(1),mla_fancy%la(1),1,1)
          call multifab_copy_c(beta_cc_fancy(1),1,beta_cc_mg(nlevs_mg),1,1,0)
          call multifab_copy_c(gamma_cc_fancy(1),1,gamma_cc_mg(nlevs_mg),1,1,0)
          call multifab_fill_boundary(beta_cc_fancy(1))
          call multifab_fill_boundary(gamma_cc_fancy(1))

          do i=1,dm
             call multifab_build_edge(alpha_fc_fancy(1,i),mla_fancy%la(1),1,0,i)
             call multifab_build_edge(  rhs_fc_fancy(1,i),mla_fancy%la(1),1,0,i)
             call multifab_build_edge(  phi_fc_fancy(1,i),mla_fancy%la(1),1,0,i)

             call multifab_copy_c(alpha_fc_fancy(1,i),1,alpha_fc_mg(nlevs_mg,i),1,1,0)
             call multifab_copy_c(rhs_fc_fancy(1,i),1,rhs_fc_mg(nlevs_mg,i),1,1,0)
             call multifab_copy_c(phi_fc_fancy(1,i),1,phi_fc_mg(nlevs_mg,i),1,1,0)
          end do

          if (dm .eq. 2) then
             allocate( beta_ed_fancy(1,1))  ! nodal
             call multifab_build_nodal(beta_ed_fancy(1,1),mla_fancy%la(1),1,0)
             call multifab_copy_c(beta_ed_fancy(1,1),1,beta_ed_mg(nlevs_mg,1),1,1,0)
          else
             allocate( beta_ed_fancy(1,3))  ! edge-based
             nodal_temp(1) = .true.
             nodal_temp(2) = .true.
             nodal_temp(3) = .false.
             call multifab_build(beta_ed_fancy(1,1),mla_fancy%la(1),1,0,nodal_temp)
             call multifab_copy_c(beta_ed_fancy(1,1),1,beta_ed_mg(nlevs_mg,1),1,1,0)
             nodal_temp(1) = .true.
             nodal_temp(2) = .false.
             nodal_temp(3) = .true.
             call multifab_build(beta_ed_fancy(1,2),mla_fancy%la(1),1,0,nodal_temp)
             call multifab_copy_c(beta_ed_fancy(1,2),1,beta_ed_mg(nlevs_mg,2),1,1,0)
             nodal_temp(1) = .false.
             nodal_temp(2) = .true.
             nodal_temp(3) = .true.
             call multifab_build(beta_ed_fancy(1,3),mla_fancy%la(1),1,0,nodal_temp)
             call multifab_copy_c(beta_ed_fancy(1,3),1,beta_ed_mg(nlevs_mg,3),1,1,0)
          end if

          vcycle_counter_temp = vcycle_counter

          call stag_mg_solver(mla_fancy,alpha_fc_fancy,beta_cc_fancy,beta_ed_fancy, &
                              gamma_cc_fancy,theta,phi_fc_fancy,rhs_fc_fancy,dx_fancy, &
                              the_bc_tower_fancy,.false.)

          vcycle_counter = vcycle_counter_temp

          do i=1,dm
             call multifab_copy_c(phi_fc_mg(nlevs_mg,i),1,phi_fc_fancy(1,i),1,1,0)
          end do

          call multifab_destroy(beta_cc_fancy(1))
          call multifab_destroy(gamma_cc_fancy(1))
          do i=1,dm
             call multifab_destroy(alpha_fc_fancy(1,i))
             call multifab_destroy(rhs_fc_fancy(1,i))
             call multifab_destroy(phi_fc_fancy(1,i))
          end do
          if (dm .eq. 2) then
             call multifab_destroy(beta_ed_fancy(1,1))
          else
             call multifab_destroy(beta_ed_fancy(1,1))
             call multifab_destroy(beta_ed_fancy(1,2))
             call multifab_destroy(beta_ed_fancy(1,3))
          end if

          call destroy(mla_fancy)
          call bc_tower_destroy(the_bc_tower_fancy)

       else

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! just do smooths at the current level as the bottom solve

          ! print out residual
          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

             do j=1,dm
                ! compute Lphi - rhs, and report residual
                call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"before    smooths at level",n,resid_temp
                end if
                
             end do
             
          end if
             
          ! control to do a different number of smooths for the bottom solver
          nsmooths = stag_mg_nsmooths_bottom

          do m=1,nsmooths

             ! do the smooths
             do color=color_start,color_end

                ! the form of weighted Jacobi we are using is
                ! phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                ! where D is the diagonal matrix containing the diagonal elements of L

                ! compute Lphi
                call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_ed_mg(n,:), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)

                ! update phi = phi + omega*D^{-1}*(rhs-Lphi)
                call stag_mg_update(la_mg(n),phi_fc_mg(n,:),rhs_fc_mg(n,:), &
                                    Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                    beta_cc_mg(n),beta_ed_mg(n,:), &
                                    gamma_cc_mg(n),dx_mg(n,:),color)

                do j=1,dm
                   
                   ! set values on physical boundaries
                   call multifab_physbc_domainvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                                  the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))
                   
                   ! fill periodic ghost cells
                   call multifab_fill_boundary(phi_fc_mg(n,j))

                   ! fill physical ghost cells
                   call multifab_physbc_macvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                               the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

                end do

             end do ! end loop over colors

          end do ! end loop over nsmooths

          !!!!!!!!!!!!!!!!!!!
          ! compute residual

          ! compute Lphi
          call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                  Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                  beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

          do j=1,dm

             ! compute Lphi - rhs, and then multiply by -1
             call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
             call multifab_mult_mult_s_c(Lphi_fc_mg(n,j),1,-1.d0,1,0)
             if (stag_mg_verbosity .ge. 3) then
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"after all smooths at level",n,resid_temp
                end if
             end if

             ! set values on physical boundaries
             call multifab_physbc_domainvel(Lphi_fc_mg(n,j),vel_bc_comp+j-1, &
                                            the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

             ! fill periodic ghost cells
             call multifab_fill_boundary(Lphi_fc_mg(n,j))

             ! fill physical ghost cells
             call multifab_physbc_macvel(Lphi_fc_mg(n,j),vel_bc_comp+j-1, &
                                         the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

          end do

       end if

       ! up the V-cycle
       do n=nlevs_mg-1,1,-1

          ! prolongate/interpolate correction to update phi
          call stag_prolongation(la_mg(n),phi_fc_mg(n,:),phi_fc_mg(n+1,:))

          do j=1,dm

             ! set values on physical boundaries
             call multifab_physbc_domainvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                            the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

             ! fill periodic ghost cells
             call multifab_fill_boundary(phi_fc_mg(n,j))

             ! fill physical ghost cells
             call multifab_physbc_macvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                         the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

          end do

          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

             do j=1,dm
                ! compute Lphi - rhs, and report residual
                call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"before    smooths at level",n,resid_temp
                end if
                
             end do
             
          end if

          do m=1,stag_mg_nsmooths_up

             do color=color_start,color_end

                ! the form of weighted Jacobi we are using is
                ! phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                ! where D is the diagonal matrix containing the diagonal elements of L

                ! compute Lphi
                call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_ed_mg(n,:), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)

                ! update phi = phi + omega*D^{-1}*(rhs-Lphi)
                call stag_mg_update(la_mg(n),phi_fc_mg(n,:),rhs_fc_mg(n,:), &
                                    Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                    beta_cc_mg(n),beta_ed_mg(n,:), &
                                    gamma_cc_mg(n),dx_mg(n,:),color)

                do j=1,dm

                   ! set values on physical boundaries
                   call multifab_physbc_domainvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                                  the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

                   ! fill periodic ghost cells
                   call multifab_fill_boundary(phi_fc_mg(n,j))

                   ! fill physical ghost cells
                   call multifab_physbc_macvel(phi_fc_mg(n,j),vel_bc_comp+j-1, &
                                               the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))
                   
                end do

             end do ! end loop over colors

          end do ! end loop over stag_mg_nsmooths_up

          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             call stag_applyop_level(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))

             do j=1,dm
                ! compute Lphi - rhs, and report residual
                call multifab_sub_sub_c(Lphi_fc_mg(n,j),1,rhs_fc_mg(n,j),1,1,0)
                resid_temp = multifab_norm_inf_c(Lphi_fc_mg(n,j),1,1)
                if (parallel_IOProcessor()) then
                   print*,"Residual for comp",j,"after all smooths at level",n,resid_temp
                end if
                
             end do
             
          end if

       end do ! end loop over nlevs_mg

       if (parallel_IOProcessor()  .and. stag_mg_verbosity .ge. 2) then
          print*,"End   V-Cycle",vcycle
       end if

       ! compute norm of residual

       ! compute Lphi
       call stag_applyop_level(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                               Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                               beta_cc_mg(1),beta_ed_mg(1,:),gamma_cc_mg(1),dx_mg(1,:))

       ! compute Lphi - rhs
       do j=1,dm
          call multifab_sub_sub_c(Lphi_fc_mg(1,j),1,rhs_fc_mg(1,j),1,1,0)
       end do

       ! compute L0 norm of Lphi - rhs and determine if problem is solved
       do j=1,dm
          resid(j) = multifab_norm_inf_c(Lphi_fc_mg(1,j),1,1)
          resid_l2(j) = multifab_norm_l2_c(Lphi_fc_mg(1,j),1,1)
          if (parallel_IOProcessor()  .and. stag_mg_verbosity .ge. 2) then
             print*,"Residual    ",j,resid(j)
             print*,"resid/resid0",j,resid(j)/resid0(j)
          end if
       end do
       if (parallel_IOProcessor()  .and. stag_mg_verbosity .ge. 1) then
          write ( *, '(A,I0,100g17.9)' ) 'StagMG: L2 |r|/|r0|: ', vcycle, &
            sqrt(sum(resid_l2(1:dm)**2))/sqrt(sum(resid0_l2(1:dm)**2)), resid_l2(1:dm)/resid0_l2(1:dm)
       end if

       resid(1:dm) = resid(1:dm)/resid0(1:dm)

       if (all(resid(1:dm) .lt. stag_mg_rel_tol)) then
          if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
             print*,"Solved in ", vcycle," staggered V-cycles"
             do j=1,dm
                print*,"resid/resid0",j,resid(j)
             end do
          end if
          exit
       end if

       if (vcycle .eq. stag_mg_max_vcycles) then
          if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
             print*,'Exiting staggered multigrid; maximum number of V-Cycles reached'
             do j=1,dm
                print*,"resid/resid0",j,resid(j)
             end do
          end if
       end if

    end do ! end loop over stag_mg_max_vcycles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Done with multigrid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j=1,dm

       ! copy solution back into phi_fc
       call multifab_copy_c(phi_fc(1,j),1,phi_fc_mg(1,j),1,1,0)

       ! set values on physical boundaries
       call multifab_physbc_domainvel(phi_fc(1,j),vel_bc_comp+j-1, &
                                      the_bc_tower%bc_tower_array(1),dx(1,:))

       ! fill periodic ghost cells
       call multifab_fill_boundary(phi_fc(1,j))

       ! fill physical ghost cells
       call multifab_physbc_macvel(phi_fc(1,j),vel_bc_comp+j-1, &
                                   the_bc_tower%bc_tower_array(1),dx(1,:))

    end do

    vcycle_counter = vcycle_counter + dm*stag_mg_max_vcycles

    ! clean up memory
    call bc_tower_destroy(the_bc_tower_mg)

    do n=1,nlevs_mg
       call multifab_destroy(beta_cc_mg(n))
       call multifab_destroy(gamma_cc_mg(n))
       do i=1,dm
          call multifab_destroy(alpha_fc_mg(n,i))
          call multifab_destroy(rhs_fc_mg(n,i))
          call multifab_destroy(phi_fc_mg(n,i))
          call multifab_destroy(Lphi_fc_mg(n,i))
          call multifab_destroy(resid_fc_mg(n,i))
       end do
       if (dm .eq. 2) then
          call multifab_destroy(beta_ed_mg(n,1))
       else
          call multifab_destroy(beta_ed_mg(n,1))
          call multifab_destroy(beta_ed_mg(n,2))
          call multifab_destroy(beta_ed_mg(n,3))
       end if
       if (n .ne. 1) then
          call destroy(la_mg(n))
       end if
    end do
    deallocate(beta_cc_mg,gamma_cc_mg,alpha_fc_mg,rhs_fc_mg)
    deallocate(phi_fc_mg,Lphi_fc_mg,resid_fc_mg,beta_ed_mg,la_mg,dx_mg)

    if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
       print*,""
       print*,"End call to stag_mg_solver"
    end if

  contains

    ! compute the number of multigrid levels assuming minwidth is the length of the
    ! smallest dimension of the smallest grid at the coarsest multigrid level
    subroutine compute_nlevs_mg(nlevs_mg,ba)
      
      integer       , intent(inout) :: nlevs_mg
      type(boxarray), intent(in   ) :: ba

      ! local
      integer :: i,d,dm,rdir,temp
      integer :: lo(get_dim(ba)),hi(get_dim(ba)),length(get_dim(ba))

      dm = get_dim(ba)

      nlevs_mg = -1

      do i=1,nboxes(ba)

         lo = lwb(get_box(ba,i))
         hi = upb(get_box(ba,i))
         length = hi-lo+1

         do d=1,dm
            temp = length(d)
            rdir = 1
            do while (mod(temp,2) .eq. 0 .and. temp/stag_mg_minwidth .ne. 1)
               temp = temp/2
               rdir = rdir+1
            end do

            if (nlevs_mg .eq. -1) then
               nlevs_mg = rdir
            else          
               nlevs_mg = min(rdir,nlevs_mg)
            end if

         end do

      end do

    end subroutine compute_nlevs_mg

    ! coarsen a cell-centered quantity
    subroutine cc_restriction(la,phi_c,phi_f,the_bc_level)

      use bl_error_module
      use define_bc_module

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(inout) :: phi_c
      type(multifab), intent(in   ) :: phi_f
      type(bc_level), intent(in   ) :: the_bc_level

      ! local
      integer :: i,dm,ng_c
      integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
      integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

      real(kind=dp_t), pointer :: acp(:,:,:,:)
      real(kind=dp_t), pointer :: afp(:,:,:,:)

      dm = get_dim(la)

      ng_c = phi_c%ng

      do i=1,nfabs(phi_c)
         acp => dataptr(phi_c, i)
         afp => dataptr(phi_f, i)
         lo_c = lwb(get_box(phi_c, i))
         hi_c = upb(get_box(phi_c, i))
         lo_f = lwb(get_box(phi_f, i))
         hi_f = upb(get_box(phi_f, i))
         select case (dm)
         case (2)
            call cc_restriction_2d(acp(:,:,1,1),afp(:,:,1,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, &
                                   the_bc_level%adv_bc_level_array(i,:,:,tran_bc_comp))
         case (3)
            call cc_restriction_3d(acp(:,:,:,1),afp(:,:,:,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, &
                                   the_bc_level%adv_bc_level_array(i,:,:,tran_bc_comp))
         end select
      end do

      call multifab_fill_boundary(phi_c)

    end subroutine cc_restriction

    subroutine cc_restriction_2d(phi_c,phi_f,ng_c,lo_c,hi_c,lo_f,hi_f,adv_bc)

      use bc_module

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_c
      real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
      real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_c:,lo_f(2)-ng_c:)
      integer        , intent(in   ) :: adv_bc(:,:)

      ! local
      integer :: i,j

      do j=lo_c(2),hi_c(2)
         do i=lo_c(1),hi_c(1)
            phi_c(i,j) = 0.25d0*(  phi_f(2*i,2*j  ) + phi_f(2*i+1,2*j  ) &
                 + phi_f(2*i,2*j+1) + phi_f(2*i+1,2*j+1) )
         end do
      end do

    end subroutine cc_restriction_2d

    subroutine cc_restriction_3d(phi_c,phi_f,ng_c,lo_c,hi_c,lo_f,hi_f,adv_bc)

      use bc_module

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_c
      real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_c:,lo_f(2)-ng_c:,lo_f(3)-ng_c:)
      integer        , intent(in   ) :: adv_bc(:,:)

      ! local
      integer :: i,j,k

      do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)
            do i=lo_c(1),hi_c(1)
               phi_c(i,j,k) = 0.125d0*(  phi_f(2*i,2*j  ,2*k  ) + phi_f(2*i+1,2*j  ,2*k  ) &
                    + phi_f(2*i,2*j+1,2*k  ) + phi_f(2*i+1,2*j+1,2*k  ) &
                    + phi_f(2*i,2*j  ,2*k+1) + phi_f(2*i+1,2*j  ,2*k+1) &
                    + phi_f(2*i,2*j+1,2*k+1) + phi_f(2*i+1,2*j+1,2*k+1) )
            end do
         end do
      end do

    end subroutine cc_restriction_3d

    ! coarsen a staggered quantity
    subroutine stag_restriction(la,phi_f,phi_c,simple_stencil_in)

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(in   ) :: phi_f(:) ! face-centered
      type(multifab), intent(inout) :: phi_c(:) ! face-centered
      logical, intent(in), optional :: simple_stencil_in

      ! local
      integer :: i,dm,ng_f,ng_c
      integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
      integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

      logical :: simple_stencil

      real(kind=dp_t), pointer :: fpx(:,:,:,:)
      real(kind=dp_t), pointer :: fpy(:,:,:,:)
      real(kind=dp_t), pointer :: fpz(:,:,:,:)
      real(kind=dp_t), pointer :: cpx(:,:,:,:)
      real(kind=dp_t), pointer :: cpy(:,:,:,:)
      real(kind=dp_t), pointer :: cpz(:,:,:,:)

      simple_stencil = .false.
      if (present(simple_stencil_in)) then
         simple_stencil = simple_stencil_in
      end if

      dm = get_dim(la)

      ng_f = phi_f(1)%ng
      ng_c = phi_c(1)%ng

      do i=1,nfabs(phi_c(1))
         fpx => dataptr(phi_f(1), i)
         fpy => dataptr(phi_f(2), i)
         cpx => dataptr(phi_c(1), i)
         cpy => dataptr(phi_c(2), i)
         lo_c = lwb(get_box(phi_c(1), i))
         hi_c = upb(get_box(phi_c(1), i))
         lo_f = lwb(get_box(phi_f(1), i))
         hi_f = upb(get_box(phi_f(1), i))
         select case (dm)
         case (2)
            call stag_restriction_2d(fpx(:,:,1,1),fpy(:,:,1,1),ng_f, &
                                     cpx(:,:,1,1),cpy(:,:,1,1),ng_c, &
                                     lo_c, hi_c, lo_f, hi_f, simple_stencil)
         case (3)
            fpz => dataptr(phi_f(3), i)
            cpz => dataptr(phi_c(3), i)
            call stag_restriction_3d(fpx(:,:,:,1),fpy(:,:,:,1),fpz(:,:,:,1),ng_f, &
                                     cpx(:,:,:,1),cpy(:,:,:,1),cpz(:,:,:,1),ng_c, &
                                     lo_c, hi_c, lo_f, hi_f, simple_stencil)
         end select
      end do

    end subroutine stag_restriction

    subroutine stag_restriction_2d(phix_f,phiy_f,ng_f,phix_c,phiy_c,ng_c, &
                                   lo_c,hi_c,lo_f,hi_f,simple_stencil)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(in   ) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
      real(kind=dp_t), intent(inout) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
      real(kind=dp_t), intent(inout) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
      logical        , intent(in   ) :: simple_stencil

      ! local
      integer :: i,j

      if (simple_stencil) then

         ! 2 point stencils
         do j=lo_c(2),hi_c(2)
            do i=lo_c(1),hi_c(1)+1
               phix_c(i,j) = 0.5d0*(phix_f(2*i,2*j) + phix_f(2*i,2*j+1))
            end do
         end do

         do j=lo_c(2),hi_c(2)+1
            do i=lo_c(1),hi_c(1)
               phiy_c(i,j) = 0.5d0*(phiy_f(2*i,2*j) + phiy_f(2*i+1,2*j))
            end do
         end do

      else

         ! 6 point stencils
         do j=lo_c(2),hi_c(2)
            do i=lo_c(1),hi_c(1)+1
               phix_c(i,j) = 0.25d0*(phix_f(2*i,2*j) + phix_f(2*i,2*j+1)) &
                    + 0.125d0*( phix_f(2*i+1,2*j) + phix_f(2*i+1,2*j+1) &
                    +phix_f(2*i-1,2*j) + phix_f(2*i-1,2*j+1))
            end do
         end do

         do j=lo_c(2),hi_c(2)+1
            do i=lo_c(1),hi_c(1)
               phiy_c(i,j) = 0.25d0*(phiy_f(2*i,2*j) + phiy_f(2*i+1,2*j)) &
                    + 0.125d0*( phiy_f(2*i,2*j+1) + phiy_f(2*i+1,2*j+1) &
                    +phiy_f(2*i,2*j-1) + phiy_f(2*i+1,2*j-1))
            end do
         end do

      end if


    end subroutine stag_restriction_2d

    subroutine stag_restriction_3d(phix_f,phiy_f,phiz_f,ng_f, &
                                   phix_c,phiy_c,phiz_c,ng_c, &
                                   lo_c,hi_c,lo_f,hi_f,simple_stencil)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(in   ) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phiz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(inout) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(inout) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(inout) :: phiz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      logical        , intent(in   ) :: simple_stencil

      ! local
      integer :: i,j,k

      if (simple_stencil) then

         ! 4 point stencils
         do k=lo_c(3),hi_c(3)
            do j=lo_c(2),hi_c(2)
               do i=lo_c(1),hi_c(1)+1
                  phix_c(i,j,k) = 0.25d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                       +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) )
               end do
            end do
         end do

         do k=lo_c(3),hi_c(3)
            do j=lo_c(2),hi_c(2)+1
               do i=lo_c(1),hi_c(1)
                  phiy_c(i,j,k) = 0.25d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                       +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) )
               end do
            end do
         end do

         do k=lo_c(3),hi_c(3)+1
            do j=lo_c(2),hi_c(2)
               do i=lo_c(1),hi_c(1)
                  phiz_c(i,j,k) = 0.25d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                       +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) )
               end do
            end do
         end do

      else

         ! 12 point stencils
         do k=lo_c(3),hi_c(3)
            do j=lo_c(2),hi_c(2)
               do i=lo_c(1),hi_c(1)+1
                  phix_c(i,j,k) = 0.125d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                       +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) ) &
                       + 0.0625* ( phix_f(2*i+1,2*j,2*k  ) + phix_f(2*i+1,2*j+1,2*k  ) &
                       +phix_f(2*i+1,2*j,2*k+1) + phix_f(2*i+1,2*j+1,2*k+1) ) &
                       + 0.0625* ( phix_f(2*i-1,2*j,2*k  ) + phix_f(2*i-1,2*j+1,2*k  ) &
                       +phix_f(2*i-1,2*j,2*k+1) + phix_f(2*i-1,2*j+1,2*k+1) )
               end do
            end do
         end do

         do k=lo_c(3),hi_c(3)
            do j=lo_c(2),hi_c(2)+1
               do i=lo_c(1),hi_c(1)
                  phiy_c(i,j,k) = 0.125d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                       +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) ) &
                       + 0.0625* ( phiy_f(2*i,2*j+1,2*k  ) + phiy_f(2*i+1,2*j+1,2*k  ) &
                       +phiy_f(2*i,2*j+1,2*k+1) + phiy_f(2*i+1,2*j+1,2*k+1) ) &
                       + 0.0625* ( phiy_f(2*i,2*j-1,2*k  ) + phiy_f(2*i+1,2*j-1,2*k  ) &
                       +phiy_f(2*i,2*j-1,2*k+1) + phiy_f(2*i+1,2*j-1,2*k+1) )
               end do
            end do
         end do

         do k=lo_c(3),hi_c(3)+1
            do j=lo_c(2),hi_c(2)
               do i=lo_c(1),hi_c(1)
                  phiz_c(i,j,k) = 0.125d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                       +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) ) &
                       + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k+1) + phiz_f(2*i+1,2*j  ,2*k+1) &
                       +phiz_f(2*i,2*j+1,2*k+1) + phiz_f(2*i+1,2*j+1,2*k+1) ) &
                       + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k-1) + phiz_f(2*i+1,2*j  ,2*k-1) &
                       +phiz_f(2*i,2*j+1,2*k-1) + phiz_f(2*i+1,2*j+1,2*k-1) )
               end do
            end do
         end do

      end if

    end subroutine stag_restriction_3d

    ! coarsen a nodal quantity
    subroutine nodal_restriction(la,phi_f,phi_c)

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(in   ) :: phi_f
      type(multifab), intent(inout) :: phi_c

      ! local
      integer :: i,dm,ng_f,ng_c
      integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
      integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

      real(kind=dp_t), pointer :: fp(:,:,:,:)
      real(kind=dp_t), pointer :: cp(:,:,:,:)

      dm = get_dim(la)

      ng_f = phi_f%ng
      ng_c = phi_c%ng

      if (ng_f .ne. 0 .or. ng_c .ne. 0) then
         call bl_error("nodal_restriction assumes 0 ghost cells")
      end if

      do i=1,nfabs(phi_c)
         fp => dataptr(phi_f, i)
         cp => dataptr(phi_c, i)
         lo_c = lwb(get_box(phi_c, i))
         hi_c = upb(get_box(phi_c, i))
         lo_f = lwb(get_box(phi_f, i))
         hi_f = upb(get_box(phi_f, i))
         select case (dm)
         case (2)
            call nodal_restriction_2d(fp(:,:,1,1),ng_f,cp(:,:,1,1),ng_c, &
                                      lo_c, hi_c, lo_f, hi_f)
         case (3)
            call bl_error("as of 2/25/13, 3D does not require nodal_restriction")
            call nodal_restriction_3d(fp(:,:,:,1),ng_f,cp(:,:,:,1),ng_c, &
                                      lo_c, hi_c, lo_f, hi_f)
         end select
      end do

      call multifab_internal_sync(phi_c)

    end subroutine nodal_restriction

    subroutine nodal_restriction_2d(phi_f,ng_f,phi_c,ng_c,lo_c,hi_c,lo_f,hi_f)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
      real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)

      ! local
      integer :: i,j

      do j=lo_c(2),hi_c(2)+1
         do i=lo_c(1),hi_c(1)+1
            phi_c(i,j) = phi_f(2*i,2*j)
         end do
      end do

    end subroutine nodal_restriction_2d

    subroutine nodal_restriction_3d(phi_f,ng_f,phi_c,ng_c,lo_c,hi_c,lo_f,hi_f)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

      ! local
      integer :: i,j,k

      do k=lo_c(3),hi_c(3)+1
         do j=lo_c(2),hi_c(2)+1
            do i=lo_c(1),hi_c(1)+1
               phi_c(i,j,k) = phi_f(2*i,2*j,2*k)
            end do
         end do
      end do

    end subroutine nodal_restriction_3d

    ! coarsen an edge-based quantity (3d only, no such thing as edges in 2d)
    subroutine edge_restriction(la,phi_f,phi_c)

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(in   ) :: phi_f(:) ! 3 components (xy, xz, yz edges)
      type(multifab), intent(inout) :: phi_c(:) ! 3 components (xy, xz, yz edges)

      ! local
      integer :: i,dm,ng_f,ng_c
      integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
      integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

      real(kind=dp_t), pointer :: fp1(:,:,:,:)
      real(kind=dp_t), pointer :: fp2(:,:,:,:)
      real(kind=dp_t), pointer :: fp3(:,:,:,:)
      real(kind=dp_t), pointer :: cp1(:,:,:,:)
      real(kind=dp_t), pointer :: cp2(:,:,:,:)
      real(kind=dp_t), pointer :: cp3(:,:,:,:)

      dm = get_dim(la)

      ng_f = phi_f(1)%ng
      ng_c = phi_c(1)%ng

      if (ng_f .ne. 0 .or. ng_c .ne. 0) then
         call bl_error("edge_restriction assumes 0 ghost cells")
      end if

      do i=1,nfabs(phi_c(1))
         fp1 => dataptr(phi_f(1), i)
         cp1 => dataptr(phi_c(1), i)
         fp2 => dataptr(phi_f(2), i)
         cp2 => dataptr(phi_c(2), i)
         fp3 => dataptr(phi_f(3), i)
         cp3 => dataptr(phi_c(3), i)
         lo_c = lwb(get_box(phi_c(1), i))
         hi_c = upb(get_box(phi_c(1), i))
         lo_f = lwb(get_box(phi_f(1), i))
         hi_f = upb(get_box(phi_f(1), i))
         select case (dm)
         case (2)
            call bl_error("no such thing as edge_restriction in 2d")
         case (3)
            call edge_restriction_3d(fp1(:,:,:,1),fp2(:,:,:,1),fp3(:,:,:,1),ng_f, &
                                     cp1(:,:,:,1),cp2(:,:,:,1),cp3(:,:,:,1),ng_c, &
                                     lo_c, hi_c, lo_f, hi_f)
         end select
      end do

      do i=1,3
         call multifab_internal_sync(phi_c(i))
      end do

    end subroutine edge_restriction

    subroutine edge_restriction_3d(phixy_f,phixz_f,phiyz_f,ng_f, &
                                   phixy_c,phixz_c,phiyz_c,ng_c,lo_c,hi_c,lo_f,hi_f)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(in   ) :: phixy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phixz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phiyz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(inout) :: phixy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(inout) :: phixz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(inout) :: phiyz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

      ! local
      integer :: i,j,k

      ! xy edges
      do k=lo_c(3),hi_c(3)
         do j=lo_c(2),hi_c(2)+1
            do i=lo_c(1),hi_c(1)+1
               phixy_c(i,j,k) = 0.5d0*(phixy_f(2*i,2*j,2*k)+phixy_f(2*i,2*j,2*k+1))
            end do
         end do
      end do

      ! xz edges
      do k=lo_c(3),hi_c(3)+1
         do j=lo_c(2),hi_c(2)
            do i=lo_c(1),hi_c(1)+1
               phixz_c(i,j,k) =  0.5d0*(phixz_f(2*i,2*j,2*k)+phixz_f(2*i,2*j+1,2*k))
            end do
         end do
      end do

      ! yz edges
      do k=lo_c(3),hi_c(3)+1
         do j=lo_c(2),hi_c(2)+1
            do i=lo_c(1),hi_c(1)
               phiyz_c(i,j,k) =  0.5d0*(phiyz_f(2*i,2*j,2*k)+phiyz_f(2*i+1,2*j,2*k))
            end do
         end do
      end do

    end subroutine edge_restriction_3d

    ! staggered prolongation from coarser grid to finer
    subroutine stag_prolongation(la,phi_f,phi_c)

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(inout) :: phi_f(:) ! face-centered
      type(multifab), intent(in   ) :: phi_c(:) ! face-centered

      ! local
      integer :: i,dm,ng_f,ng_c
      integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
      integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

      real(kind=dp_t), pointer :: fpx(:,:,:,:)
      real(kind=dp_t), pointer :: fpy(:,:,:,:)
      real(kind=dp_t), pointer :: fpz(:,:,:,:)
      real(kind=dp_t), pointer :: cpx(:,:,:,:)
      real(kind=dp_t), pointer :: cpy(:,:,:,:)
      real(kind=dp_t), pointer :: cpz(:,:,:,:)

      dm = get_dim(la)

      ng_f = phi_f(1)%ng
      ng_c = phi_c(1)%ng

      do i=1,nfabs(phi_f(1))
         fpx => dataptr(phi_f(1), i)
         fpy => dataptr(phi_f(2), i)
         cpx => dataptr(phi_c(1), i)
         cpy => dataptr(phi_c(2), i)
         lo_c = lwb(get_box(phi_c(1), i))
         hi_c = upb(get_box(phi_c(1), i))
         lo_f = lwb(get_box(phi_f(1), i))
         hi_f = upb(get_box(phi_f(1), i))
         select case (dm)
         case (2)
            call stag_prolongation_2d(fpx(:,:,1,1),fpy(:,:,1,1),ng_f, &
                                      cpx(:,:,1,1),cpy(:,:,1,1),ng_c, &
                                      lo_c, hi_c, lo_f, hi_f)
         case (3)
            fpz => dataptr(phi_f(3), i)
            cpz => dataptr(phi_c(3), i)
            call stag_prolongation_3d(fpx(:,:,:,1),fpy(:,:,:,1),fpz(:,:,:,1),ng_f, &
                                      cpx(:,:,:,1),cpy(:,:,:,1),cpz(:,:,:,1),ng_c, &
                                      lo_c, hi_c, lo_f, hi_f)
         end select
      end do

    end subroutine stag_prolongation

    subroutine stag_prolongation_2d(phix_f,phiy_f,ng_f,phix_c,phiy_c,ng_c, &
                                    lo_c,hi_c,lo_f,hi_f)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(inout) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
      real(kind=dp_t), intent(inout) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
      real(kind=dp_t), intent(in   ) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)

      ! local
      integer :: i,j,ioff,joff

      do j=lo_f(2),hi_f(2)
         do i=lo_f(1),hi_f(1)+1

            if (mod(j,2) .eq. 0) then
               joff = -1
            else
               joff = 1
            end if

            if (mod(i,2) .eq. 0) then

               ! linear interpolation
               phix_f(i,j) = phix_f(i,j) + 0.75d0*phix_c(i/2,j/2) + 0.25d0*phix_c(i/2,j/2+joff)

            else

               ! bilinear interpolation
               phix_f(i,j) = phix_f(i,j) + 0.375d0*phix_c(i/2  ,j/2) &
                    + 0.125d0*phix_c(i/2  ,j/2+joff) &
                    + 0.375d0*phix_c(i/2+1,j/2) &
                    + 0.125d0*phix_c(i/2+1,j/2+joff)

            end if

         end do
      end do

      do j=lo_f(2),hi_f(2)+1
         do i=lo_f(1),hi_f(1)

            if (mod(i,2) .eq. 0) then
               ioff = -1
            else
               ioff = 1
            end if

            if (mod(j,2) .eq. 0) then

               ! linear interpolation
               phiy_f(i,j) = phiy_f(i,j) + 0.75d0*phiy_c(i/2,j/2) + 0.25d0*phiy_c(i/2+ioff,j/2)

            else

               ! bilinear interpolation
               phiy_f(i,j) = phiy_f(i,j) + 0.375d0*phiy_c(i/2,j/2  ) &
                    + 0.125d0*phiy_c(i/2+ioff,j/2  ) &
                    + 0.375d0*phiy_c(i/2,j/2+1) &
                    + 0.125d0*phiy_c(i/2+ioff,j/2+1)

            end if

         end do
      end do

    end subroutine stag_prolongation_2d

    subroutine stag_prolongation_3d(phix_f,phiy_f,phiz_f,ng_f,phix_c,phiy_c,phiz_c,ng_c, &
                                    lo_c,hi_c,lo_f,hi_f)

      integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
      real(kind=dp_t), intent(inout) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(inout) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(inout) :: phiz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(in   ) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
      real(kind=dp_t), intent(in   ) :: phiz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

      ! local
      integer :: i,j,k,ioff,joff,koff

      real(kind=dp_t), parameter :: nine16 = 9.d0/16.d0
      real(kind=dp_t), parameter :: three16 = 3.d0/16.d0
      real(kind=dp_t), parameter :: one16 = 1.d0/16.d0
      real(kind=dp_t), parameter :: nine32 = 9.d0/32.d0
      real(kind=dp_t), parameter :: three32 = 3.d0/32.d0
      real(kind=dp_t), parameter :: one32 = 1.d0/32.d0

      do k=lo_f(3),hi_f(3)
         do j=lo_f(2),hi_f(2)
            do i=lo_f(1),hi_f(1)+1

               if (mod(j,2) .eq. 0) then
                  joff = -1
               else
                  joff = 1
               end if

               if (mod(k,2) .eq. 0) then
                  koff = -1
               else
                  koff = 1
               end if

               if (mod(i,2) .eq. 0) then
                  ! bilinear in the yz plane
                  phix_f(i,j,k) = phix_f(i,j,k) &
                       + nine16*phix_c(i/2,j/2,k/2) &
                       + three16*phix_c(i/2,j/2+joff,k/2) &
                       + three16*phix_c(i/2,j/2,k/2+koff) &
                       + one16*phix_c(i/2,j/2+joff,k/2+koff)
               else
                  ! bilinear in the yz plane, linear in x
                  phix_f(i,j,k) = phix_f(i,j,k) &
                       + nine32*phix_c(i/2,j/2,k/2) &
                       + three32*phix_c(i/2,j/2+joff,k/2) &
                       + three32*phix_c(i/2,j/2,k/2+koff) &
                       + one32*phix_c(i/2,j/2+joff,k/2+koff)&
                       + nine32*phix_c(i/2+1,j/2,k/2) &
                       + three32*phix_c(i/2+1,j/2+joff,k/2) &
                       + three32*phix_c(i/2+1,j/2,k/2+koff) &
                       + one32*phix_c(i/2+1,j/2+joff,k/2+koff)
               end if

            end do
         end do
      end do

      do k=lo_f(3),hi_f(3)
         do j=lo_f(2),hi_f(2)+1
            do i=lo_f(1),hi_f(1)

               if (mod(i,2) .eq. 0) then
                  ioff = -1
               else
                  ioff = 1
               end if

               if (mod(k,2) .eq. 0) then
                  koff = -1
               else
                  koff = 1
               end if

               if (mod(j,2) .eq. 0) then
                  ! bilinear in the xz plane
                  phiy_f(i,j,k) = phiy_f(i,j,k) &
                       + nine16*phiy_c(i/2,j/2,k/2) &
                       + three16*phiy_c(i/2+ioff,j/2,k/2) &
                       + three16*phiy_c(i/2,j/2,k/2+koff) &
                       + one16*phiy_c(i/2+ioff,j/2,k/2+koff)
               else
                  ! bilinear in the yz plane, linear in y
                  phiy_f(i,j,k) = phiy_f(i,j,k) &
                       + nine32*phiy_c(i/2,j/2,k/2) &
                       + three32*phiy_c(i/2+ioff,j/2,k/2) &
                       + three32*phiy_c(i/2,j/2,k/2+koff) &
                       + one32*phiy_c(i/2+ioff,j/2,k/2+koff)&
                       + nine32*phiy_c(i/2,j/2+1,k/2) &
                       + three32*phiy_c(i/2+ioff,j/2+1,k/2) &
                       + three32*phiy_c(i/2,j/2+1,k/2+koff) &
                       + one32*phiy_c(i/2+ioff,j/2+1,k/2+koff)
               end if

            end do
         end do
      end do

      do k=lo_f(3),hi_f(3)+1
         do j=lo_f(2),hi_f(2)
            do i=lo_f(1),hi_f(1)

               if (mod(i,2) .eq. 0) then
                  ioff = -1
               else
                  ioff = 1
               end if

               if (mod(j,2) .eq. 0) then
                  joff = -1
               else
                  joff = 1
               end if

               if (mod(k,2) .eq. 0) then
                  ! bilinear in the xy plane
                  phiz_f(i,j,k) = phiz_f(i,j,k) &
                       + nine16*phiz_c(i/2,j/2,k/2) &
                       + three16*phiz_c(i/2+ioff,j/2,k/2) &
                       + three16*phiz_c(i/2,j/2+joff,k/2) &
                       + one16*phiz_c(i/2+ioff,j/2+joff,k/2)
               else
                  ! bilinear in the xy plane, linear in z
                  phiz_f(i,j,k) = phiz_f(i,j,k) &
                       + nine32*phiz_c(i/2,j/2,k/2) &
                       + three32*phiz_c(i/2+ioff,j/2,k/2) &
                       + three32*phiz_c(i/2,j/2+joff,k/2) &
                       + one32*phiz_c(i/2+ioff,j/2+joff,k/2)&
                       + nine32*phiz_c(i/2,j/2,k/2+1) &
                       + three32*phiz_c(i/2+ioff,j/2,k/2+1) &
                       + three32*phiz_c(i/2,j/2+joff,k/2+1) &
                       + one32*phiz_c(i/2+ioff,j/2+joff,k/2+1)
               end if

            end do
         end do
      end do

    end subroutine stag_prolongation_3d

    ! finish the Jacobi iteration by multiplying the residual by the inverse
    ! of the diagonal-element-only matrix
    subroutine stag_mg_update(la,phi_fc,rhs_fc,Lphi_fc, &
                              alpha_fc,beta_cc,beta_ed,gamma_cc,dx,color_in)

      type(layout)  , intent(in   ) :: la
      type(multifab), intent(inout) :: phi_fc(:)   ! face-centered
      type(multifab), intent(in   ) :: rhs_fc(:)   ! face-centered
      type(multifab), intent(in   ) :: Lphi_fc(:)  ! face-centered
      type(multifab), intent(in   ) :: alpha_fc(:) ! face-centered
      type(multifab), intent(in   ) :: beta_cc     ! cell-centered
      type(multifab), intent(in   ) :: beta_ed(:)  ! edge-based
      type(multifab), intent(in   ) :: gamma_cc    ! cell-centered
      real(kind=dp_t),intent(in   ) :: dx(:)
      integer        , intent(in   ), optional :: color_in

      ! local
      integer :: i,dm,ng_p,ng_r,ng_l,ng_a,ng_b,ng_n,ng_g,ng_e
      integer :: lo(get_dim(la)), hi(get_dim(la))
      integer :: color

      real(kind=dp_t), pointer :: ppx(:,:,:,:)
      real(kind=dp_t), pointer :: ppy(:,:,:,:)
      real(kind=dp_t), pointer :: ppz(:,:,:,:)
      real(kind=dp_t), pointer :: rpx(:,:,:,:)
      real(kind=dp_t), pointer :: rpy(:,:,:,:)
      real(kind=dp_t), pointer :: rpz(:,:,:,:)

      real(kind=dp_t), pointer :: lpx(:,:,:,:)
      real(kind=dp_t), pointer :: lpy(:,:,:,:)
      real(kind=dp_t), pointer :: lpz(:,:,:,:)
      real(kind=dp_t), pointer :: apx(:,:,:,:)
      real(kind=dp_t), pointer :: apy(:,:,:,:)
      real(kind=dp_t), pointer :: apz(:,:,:,:)
      real(kind=dp_t), pointer ::  bp(:,:,:,:)
      real(kind=dp_t), pointer :: bp1(:,:,:,:)
      real(kind=dp_t), pointer :: bp2(:,:,:,:)
      real(kind=dp_t), pointer :: bp3(:,:,:,:)
      real(kind=dp_t), pointer :: bnp(:,:,:,:)
      real(kind=dp_t), pointer ::  kp(:,:,:,:)

      dm = get_dim(la)

      if (present(color_in)) then
         color = color_in
      else
         color = 0
      end if

      ng_p = phi_fc(1)%ng
      ng_r = rhs_fc(1)%ng
      ng_l = Lphi_fc(1)%ng
      ng_a = alpha_fc(1)%ng
      ng_b = beta_cc%ng
      ng_g = gamma_cc%ng

      do i=1,nfabs(Lphi_fc(1))
         ppx => dataptr(phi_fc(1), i)
         ppy => dataptr(phi_fc(2), i)
         rpx => dataptr(rhs_fc(1), i)
         rpy => dataptr(rhs_fc(2), i)
         lpx => dataptr(Lphi_fc(1), i)
         lpy => dataptr(Lphi_fc(2), i)
         apx => dataptr(alpha_fc(1), i)
         apy => dataptr(alpha_fc(2), i)
         bp  => dataptr(beta_cc, i)
         kp  => dataptr(gamma_cc, i)
         lo = lwb(get_box(Lphi_fc(1), i))
         hi = upb(get_box(Lphi_fc(1), i))
         select case(dm)
         case (2)
            ng_n = beta_ed(1)%ng
            bnp => dataptr(beta_ed(1), i)
            call stag_mg_update_2d(ppx(:,:,1,1),ppy(:,:,1,1),ng_p, &
                                   rpx(:,:,1,1),rpy(:,:,1,1),ng_r, &
                                   lpx(:,:,1,1),lpy(:,:,1,1),ng_l, &
                                   apx(:,:,1,1),apy(:,:,1,1),ng_a, &
                                   bp(:,:,1,1),ng_b,bnp(:,:,1,1),ng_n, &
                                   kp(:,:,1,1), ng_g,lo,hi,dx,color)
         case (3)
            ng_e = beta_ed(1)%ng
            ppz => dataptr(phi_fc(3), i)
            rpz => dataptr(rhs_fc(3), i)
            lpz => dataptr(Lphi_fc(3), i)
            apz => dataptr(alpha_fc(3), i)
            bp1 => dataptr(beta_ed(1), i)
            bp2 => dataptr(beta_ed(2), i)
            bp3 => dataptr(beta_ed(3), i)
            call stag_mg_update_3d(ppx(:,:,:,1),ppy(:,:,:,1),ppz(:,:,:,1),ng_p, &
                                   rpx(:,:,:,1),rpy(:,:,:,1),rpz(:,:,:,1),ng_r, &
                                   lpx(:,:,:,1),lpy(:,:,:,1),lpz(:,:,:,1),ng_l, &
                                   apx(:,:,:,1),apy(:,:,:,1),apz(:,:,:,1),ng_a, &
                                   bp(:,:,:,1),ng_b, &
                                   bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                                   kp(:,:,:,1),ng_g,lo,hi,dx,color)

         end select
      end do

    end subroutine stag_mg_update

    subroutine stag_mg_update_2d(phix,phiy,ng_p,rhsx,rhsy,ng_r, &
                                 Lpx,Lpy,ng_l,alphax,alphay,ng_a,beta,ng_b, &
                                 beta_ed,ng_n,gamma,ng_g,lo,hi,dx,color)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_r,ng_l,ng_a,ng_b,ng_n,ng_g
      real(kind=dp_t), intent(inout) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) ::    rhsx(lo(1)-ng_r:,lo(2)-ng_r:)
      real(kind=dp_t), intent(in   ) ::    rhsy(lo(1)-ng_r:,lo(2)-ng_r:)
      real(kind=dp_t), intent(in   ) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:)
      real(kind=dp_t), intent(in   ) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:)
      real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(in   ) :: beta_ed(lo(1)-ng_n:,lo(2)-ng_n:)
      real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: color

      ! local
      integer :: i,j

      real(kind=dp_t) :: fac, dxsq, fourthirds
      real(kind=dp_t) :: b,c

      ! coloring parameters
      logical :: do_x, do_y
      integer :: offset, ioff

      do_x = .true.
      do_y = .true.
      offset = 1

      if (color .eq. 1 .or. color .eq. 2) then
         do_y = .false.
         offset = 2
      else if (color .eq. 3 .or. color .eq. 4) then
         do_x = .false.
         offset = 2
      end if

      dxsq = dx(1)**2
      fourthirds = 4.d0/3.d0

      if (visc_type .eq. -1) then

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j) + &
                       (beta(i,j)+beta(i-1,j)+beta_ed(i,j)+beta_ed(i,j+1))/dxsq

                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j) + &
                       (beta(i,j)+beta(i,j-1)+beta_ed(i,j)+beta_ed(i+1,j))/dxsq

                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      else if (visc_type .eq. 1) then

         b = beta(lo(1),lo(2))

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j) + 4.d0*b/dxsq
                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j) + 4.d0*b/dxsq
                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      else if (visc_type .eq. -2) then

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j) + &
                       (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_ed(i,j)+beta_ed(i,j+1))/dxsq

                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j) + &
                       (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_ed(i,j)+beta_ed(i+1,j))/dxsq

                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      else if (visc_type .eq. 2) then

         b = beta(lo(1),lo(2))

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j)+6.d0*b/dxsq
                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j)+6.d0*b/dxsq
                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      else if (visc_type .eq. -3) then

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j) + &
                       ( fourthirds*beta(i,j)+gamma(i,j) &
                       +fourthirds*beta(i-1,j)+gamma(i-1,j) &
                       +beta_ed(i,j)+beta_ed(i,j+1))/dxsq

                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j) + &
                       ( fourthirds*beta(i,j)+gamma(i,j) &
                       +fourthirds*beta(i,j-1)+gamma(i,j-1) &
                       +beta_ed(i,j)+beta_ed(i+1,j))/dxsq

                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      else if (visc_type .eq. 3) then

         b = beta(lo(1),lo(2))
         c = gamma(lo(1),lo(2))

         if (do_x) then

            do j=lo(2),hi(2)
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1)+1,offset

                  fac = alphax(i,j) + (14.d0/3.d0*b+2.d0*c)/dxsq
                  phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

               end do
            end do

         end if

         if (do_y) then

            do j=lo(2),hi(2)+1
               ioff = 0
               if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
               do i=lo(1)+ioff,hi(1),offset

                  fac = alphay(i,j) + (14.d0/3.d0*b+2.d0*c)/dxsq
                  phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

               end do
            end do

         end if

      end if

    end subroutine stag_mg_update_2d

    subroutine stag_mg_update_3d(phix,phiy,phiz,ng_p,rhsx,rhsy,rhsz,ng_r, &
                                 Lpx,Lpy,Lpz,ng_l,alphax,alphay,alphaz,ng_a,beta,ng_b, &
                                 beta_xy,beta_xz,beta_yz,ng_e,gamma,ng_g, &
                                 lo,hi,dx,color)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_r,ng_l,ng_a,ng_b,ng_e,ng_g
      real(kind=dp_t), intent(inout) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) ::    phiz(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) ::    rhsx(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      real(kind=dp_t), intent(in   ) ::    rhsy(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      real(kind=dp_t), intent(in   ) ::    rhsz(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      real(kind=dp_t), intent(in   ) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
      real(kind=dp_t), intent(in   ) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
      real(kind=dp_t), intent(in   ) ::     Lpz(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
      real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) ::  alphaz(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(in   ) :: beta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
      real(kind=dp_t), intent(in   ) :: beta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
      real(kind=dp_t), intent(in   ) :: beta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
      real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: color

      ! local
      integer :: i,j,k

      real(kind=dp_t) :: fac, dxsq, fourthirds
      real(kind=dp_t) :: b,c

      ! coloring parameters
      logical :: do_x, do_y, do_z
      integer :: offset, ioff

      do_x = .true.
      do_y = .true.
      do_z = .true.
      offset = 1

      if (color .eq. 1 .or. color .eq. 2) then
         do_y = .false.
         do_z = .false.
         offset = 2
      else if (color .eq. 3 .or. color .eq. 4) then
         do_x = .false.
         do_z = .false.
         offset = 2
      else if (color .eq. 5 .or. color .eq. 6) then
         do_x = .false.
         do_y = .false.
         offset = 2
      end if

      dxsq = dx(1)**2
      fourthirds = 4.d0/3.d0

      if (visc_type .eq. -1) then

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k) + &
                          ( beta(i,j,k)+beta(i-1,j,k) &
                          +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                          +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k) + &
                          ( beta(i,j,k)+beta(i,j-1,k) &
                          +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k) + &
                          ( beta(i,j,k)+beta(i,j,k-1) &
                          +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      else if (visc_type .eq. 1) then

         b = beta(lo(1),lo(2),lo(3))

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k) + 6.d0*b/dxsq
                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k) + 6.d0*b/dxsq
                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k) + 6.d0*b/dxsq
                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      else if (visc_type .eq. -2) then

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k) + &
                          ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                          +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                          +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k) + &
                          ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                          +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k) + &
                          ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                          +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      else if (visc_type .eq. 2) then

         b = beta(lo(1),lo(2),lo(3))

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k) + 8.d0*b/dxsq
                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k) + 8.d0*b/dxsq
                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k) + 8.d0*b/dxsq
                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      else if (visc_type .eq. -3) then

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k) + &
                          ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                          +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                          +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                          +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k) + &
                          ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                          +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                          +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k) + &
                          ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                          +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1) &
                          +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                          +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      else if (visc_type .eq. 3) then

         b = beta(lo(1),lo(2),lo(3))
         c = gamma(lo(1),lo(2),lo(3))

         if (do_x) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1)+1,offset

                     fac = alphax(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                     phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_y) then

            do k=lo(3),hi(3)
               do j=lo(2),hi(2)+1
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphay(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                     phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                  end do
               end do
            end do

         end if

         if (do_z) then

            do k=lo(3),hi(3)+1
               do j=lo(2),hi(2)
                  ioff = 0
                  if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                  do i=lo(1)+ioff,hi(1),offset

                     fac = alphaz(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                     phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                  end do
               end do
            end do

         end if

      end if

    end subroutine stag_mg_update_3d

  end subroutine stag_mg_solver

end module stag_mg_solver_module
