module stag_mg_solver_module

  use ml_layout_module
  use stag_mg_util_module
  use stag_applyop_module
  use convert_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_module

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
  ! alpha_cc, beta_cc, and gamma_cc are cell-centered
  ! phi_fc and rhs_fc are face-centered
  ! phi_fc must come in initialized to some value, preferably a reasonable guess
  subroutine stag_mg_solver(mla,alpha_cc,beta_cc,gamma_cc,theta, &
                            phi_fc,rhs_fc,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: alpha_cc(:) ! cell-centered
    type(multifab) , intent(in   ) :: beta_cc(:) ! cell-centered
    type(multifab) , intent(in   ) :: gamma_cc(:) ! cell-centered
    type(multifab) , intent(inout) :: phi_fc(:,:) ! face-centered
    type(multifab) , intent(in   ) :: rhs_fc(:,:) ! face-centered
    real(kind=dp_t), intent(in   ) :: theta,dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: dm, i, j, vcycle, m, n, nlevs_mg
    integer :: color,color_start,color_end,nsmooths

    ! stores initial residual and current residual
    real(kind=dp_t) :: resid0(mla%dim),resid(mla%dim),resid_temp,resid0_l2(mla%dim),resid_l2(mla%dim),temp

    ! hold the problem domain and boxarray at level 1 as well as the current multigrid level
    type(box)      :: pd_base,pd
    type(boxarray) :: ba_base,ba

    ! cell-centered multifabs
    type(multifab), allocatable :: alpha_cc_mg(:), beta_cc_mg(:), gamma_cc_mg(:)

    ! face-centered multifabs
    type(multifab), allocatable :: alpha_fc_mg(:,:), rhs_fc_mg(:,:), phi_fc_mg(:,:)
    type(multifab), allocatable :: Lphi_fc_mg(:,:), resid_fc_mg(:,:)

    ! nodal multifabs
    type(multifab), allocatable :: beta_nd_mg(:)

    ! edge multifab
    type(multifab), allocatable :: beta_ed_mg(:,:)

    ! the layout at each level of multigrid
    type(layout)   , allocatable :: la_mg(:)

    ! grid spacing at each level of multigrid
    real(kind=dp_t), allocatable :: dx_mg(:,:)

    logical :: nodal_temp(mla%dim)

    type(bc_tower) :: the_bc_tower_mg

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
    ! and stag_mg_maxlevs is the most number of levels allowed
    call compute_nlevs_mg(nlevs_mg,ba_base)

    ! allocate multifabs used in multigrid coarsening
    allocate(alpha_cc_mg(nlevs_mg))    ! cell-centered
    allocate( beta_cc_mg(nlevs_mg))    ! cell-centered
    allocate(gamma_cc_mg(nlevs_mg))    ! cell-centered
    allocate(alpha_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(  rhs_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(  phi_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate( Lphi_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate(resid_fc_mg(nlevs_mg,dm)) ! face-centered
    allocate( beta_nd_mg(nlevs_mg))    ! nodal
    allocate( beta_ed_mg(nlevs_mg,3))  ! edge-based

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
       call multifab_build(alpha_cc_mg(n),la_mg(n),1,1)
       call multifab_build( beta_cc_mg(n),la_mg(n),1,1)
       call multifab_build(gamma_cc_mg(n),la_mg(n),1,1)

       do i=1,dm
          call multifab_build_edge(alpha_fc_mg(n,i),la_mg(n),1,0,i)
          call multifab_build_edge(  rhs_fc_mg(n,i),la_mg(n),1,0,i)
          call multifab_build_edge(  phi_fc_mg(n,i),la_mg(n),1,1,i)
          call multifab_build_edge( Lphi_fc_mg(n,i),la_mg(n),1,1,i)
          call multifab_build_edge(resid_fc_mg(n,i),la_mg(n),1,0,i)
       end do

       if (dm .eq. 2) then
          ! build beta_ng_mg
          call multifab_build_nodal(beta_nd_mg(n),la_mg(n),1,0)
       else
          ! build beta_ed_mg
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

    ! copy level 1 cc coefficients into mg array of coefficients
    call multifab_copy_c(alpha_cc_mg(1),1,alpha_cc(1),1,1,1)
    call multifab_copy_c( beta_cc_mg(1),1,beta_cc(1) ,1,1,1)
    call multifab_copy_c(gamma_cc_mg(1),1,gamma_cc(1),1,1,1)
    ! multiply alpha_cc_mg by theta
    call multifab_mult_mult_s_c(alpha_cc_mg(1),1,theta,1,1)

    ! average to face/edge/nodes at level 1 and then coarsen coefficients

    ! compute alpha at faces at level 1
    call average_cc_to_face(1,alpha_cc_mg,alpha_fc_mg, &
                            1,dm+2,1,the_bc_tower_mg%bc_tower_array)

    if (dm .eq. 2) then
       ! compute beta at nodes at level 1
       call average_cc_to_node(1,beta_cc_mg,beta_nd_mg, &
                               1,dm+2,1,the_bc_tower_mg%bc_tower_array)
    else
       ! compute beta at edges at level 1
       call average_cc_to_edge(1,beta_cc_mg,beta_ed_mg, &
                               1,dm+2,1,the_bc_tower_mg%bc_tower_array)
    end if

    ! coarsen coefficients
    do n=2,nlevs_mg
       call cc_restriction(la_mg(n), beta_cc_mg(n), beta_cc_mg(n-1),the_bc_tower_mg%bc_tower_array(n-1))
       call cc_restriction(la_mg(n),gamma_cc_mg(n),gamma_cc_mg(n-1),the_bc_tower_mg%bc_tower_array(n-1))
       call stag_restriction(la_mg(n),alpha_fc_mg(n-1,:),alpha_fc_mg(n,:),.true.)
       if (dm .eq. 2) then
          call nodal_restriction(la_mg(n),beta_nd_mg(n-1),beta_nd_mg(n))
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
       call multifab_physbc_domainvel(phi_fc_mg(1,i),1,i,1,the_bc_tower_mg%bc_tower_array(1),dx_mg(1,:))
       
       ! fill periodic ghost cells
       call multifab_fill_boundary(phi_fc_mg(1,i))

       ! fill physical ghost cells
       call multifab_physbc_macvel(phi_fc_mg(1,i),1,i,1,the_bc_tower_mg%bc_tower_array(1),dx_mg(1,:))

    end do

    ! set rhs_fc_mg at level 1 by copying in passed-in rhs_fc
    do i=1,dm
       call multifab_copy_c(rhs_fc_mg(1,i),1,rhs_fc(1,i),1,1,0)
    end do

    ! compute norm of initial residual
    ! first compute Lphi
    if (dm .eq. 2) then
       call stag_applyop_2d(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                            Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                            beta_cc_mg(1),beta_nd_mg(1),gamma_cc_mg(1),dx_mg(1,:))
    else
       call stag_applyop_3d(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                            Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                            beta_cc_mg(1),beta_ed_mg(1,:),gamma_cc_mg(1),dx_mg(1,:))
    end if

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
       do n=1,nlevs_mg

          ! print out residual
          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             if (dm .eq. 2) then
                call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_nd_mg(n),gamma_cc_mg(n),dx_mg(n,:))
             else
                call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))
             end if

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
          if (n .eq. nlevs_mg) then
             nsmooths = stag_mg_nsmooths_bottom
          end if

          do m=1,nsmooths

             ! do the smooths
             do color=color_start,color_end

                ! the form of weighted Jacobi we are using is
                ! phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                ! where D is the diagonal matrix containing the diagonal elements of L

                ! compute Lphi
                if (dm .eq. 2) then
                   call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_nd_mg(n), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)
                else
                   call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_ed_mg(n,:), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)
                end if

                ! update phi = phi + omega*D^{-1}*(rhs-Lphi)
                call stag_mg_update(la_mg(n),phi_fc_mg(n,:),rhs_fc_mg(n,:), &
                                    Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                    beta_cc_mg(n),beta_nd_mg(n),beta_ed_mg(n,:), &
                                    gamma_cc_mg(n),dx_mg(n,:),color)

                do j=1,dm
                   
                   ! set values on physical boundaries
                   call multifab_physbc_domainvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))
                   
                   ! fill periodic ghost cells
                   call multifab_fill_boundary(phi_fc_mg(n,j))

                   ! fill physical ghost cells
                   call multifab_physbc_macvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

                end do

             end do ! end loop over colors

          end do ! end loop over nsmooths

          ! compute residual

          ! compute Lphi
          if (dm .eq. 2) then
             call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                  Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                  beta_cc_mg(n),beta_nd_mg(n),gamma_cc_mg(n),dx_mg(n,:))
          else
             call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                  Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                  beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))
          end if

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
             call multifab_physbc_domainvel(Lphi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

             ! fill periodic ghost cells
             call multifab_fill_boundary(Lphi_fc_mg(n,j))

             ! fill physical ghost cells
             call multifab_physbc_macvel(Lphi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

          end do

          if (n .ne. nlevs_mg) then

             ! restrict/coarsen residual and put it in rhs_fc
             call stag_restriction(la_mg(n),Lphi_fc_mg(n,:),rhs_fc_mg(n+1,:))

             do j=1,dm
                ! set residual to zero on physical boundaries
                call multifab_physbc_domainvel(rhs_fc_mg(n+1,j),1,j,1,the_bc_tower_mg%bc_tower_array(n+1),dx_mg(n+1,:))
             end do

          end if

       end do ! end loop over nlevs_mg

       ! up the V-cycle
       do n=nlevs_mg-1,1,-1

          ! prolongate/interpolate correction to update phi
          call stag_prolongation(la_mg(n),phi_fc_mg(n,:),phi_fc_mg(n+1,:))

          do j=1,dm

             ! set values on physical boundaries
             call multifab_physbc_domainvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

             ! fill periodic ghost cells
             call multifab_fill_boundary(phi_fc_mg(n,j))

             ! fill physical ghost cells
             call multifab_physbc_macvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

          end do

          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             if (dm .eq. 2) then
                call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_nd_mg(n),gamma_cc_mg(n),dx_mg(n,:))
             else
                call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))
             end if

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
                if (dm .eq. 2) then
                   call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_nd_mg(n), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)
                else
                   call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                        Lphi_fc_mg(n,:), &
                                        alpha_fc_mg(n,:),beta_cc_mg(n),beta_ed_mg(n,:), &
                                        gamma_cc_mg(n),dx_mg(n,:),color)
                end if

                ! update phi = phi + omega*D^{-1}*(rhs-Lphi)
                call stag_mg_update(la_mg(n),phi_fc_mg(n,:),rhs_fc_mg(n,:), &
                                    Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                    beta_cc_mg(n),beta_nd_mg(n),beta_ed_mg(n,:), &
                                    gamma_cc_mg(n),dx_mg(n,:),color)

                do j=1,dm

                   ! set values on physical boundaries
                   call multifab_physbc_domainvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))

                   ! fill periodic ghost cells
                   call multifab_fill_boundary(phi_fc_mg(n,j))

                   ! fill physical ghost cells
                   call multifab_physbc_macvel(phi_fc_mg(n,j),1,j,1,the_bc_tower_mg%bc_tower_array(n),dx_mg(n,:))
                   
                end do

             end do ! end loop over colors

          end do ! end loop over stag_mg_nsmooths_up

          if (stag_mg_verbosity .ge. 3) then

             ! compute Lphi
             if (dm .eq. 2) then
                call stag_applyop_2d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_nd_mg(n),gamma_cc_mg(n),dx_mg(n,:))
             else
                call stag_applyop_3d(la_mg(n),the_bc_tower_mg%bc_tower_array(n),phi_fc_mg(n,:), &
                                     Lphi_fc_mg(n,:),alpha_fc_mg(n,:), &
                                     beta_cc_mg(n),beta_ed_mg(n,:),gamma_cc_mg(n),dx_mg(n,:))
             end if

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
       if (dm .eq. 2) then
          call stag_applyop_2d(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                               Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                               beta_cc_mg(1),beta_nd_mg(1),gamma_cc_mg(1),dx_mg(1,:))
       else
          call stag_applyop_3d(la_mg(1),the_bc_tower_mg%bc_tower_array(1),phi_fc_mg(1,:), &
                               Lphi_fc_mg(1,:),alpha_fc_mg(1,:), &
                               beta_cc_mg(1),beta_ed_mg(1,:),gamma_cc_mg(1),dx_mg(1,:))
       end if

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
       call multifab_physbc_domainvel(phi_fc(1,j),1,j,1,the_bc_tower%bc_tower_array(1),dx(1,:))

       ! fill periodic ghost cells
       call multifab_fill_boundary(phi_fc(1,j))

       ! fill physical ghost cells
       call multifab_physbc_macvel(phi_fc(1,j),1,j,1,the_bc_tower%bc_tower_array(1),dx(1,:))

    end do

    num_mg_vcycles = num_mg_vcycles + dm*stag_mg_max_vcycles

    ! clean up memory
    call bc_tower_destroy(the_bc_tower_mg)

    do n=1,nlevs_mg
       call multifab_destroy(alpha_cc_mg(n))
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
          call multifab_destroy(beta_nd_mg(n))
       else
          do i=1,3
             call multifab_destroy(beta_ed_mg(n,i))
          end do
       end if
       if (n .ne. 1) then
          call destroy(la_mg(n))
       end if
    end do
    deallocate(alpha_cc_mg,beta_cc_mg,gamma_cc_mg,alpha_fc_mg,rhs_fc_mg)
    deallocate(phi_fc_mg,Lphi_fc_mg,resid_fc_mg,beta_nd_mg,la_mg,dx_mg)

  end subroutine stag_mg_solver

  ! compute the number of multigrid levels assuming minwidth is the length of the
  ! smallest dimension of the smallest grid at the coarsest multigrid level
  subroutine compute_nlevs_mg(nlevs_mg,ba)

    use probin_module, only: stag_mg_minwidth, stag_mg_maxlevs

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

    nlevs_mg = min(nlevs_mg,stag_mg_maxlevs)

  end subroutine compute_nlevs_mg

end module stag_mg_solver_module
