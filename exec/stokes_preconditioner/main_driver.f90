!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! see notes in FluctHydro/varden-staggered/test_precon/doc/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! test_type>0: call gmres
! test_type=0: check discretization accuracy
! test_type<0: call apply_precon only
! abs(test_type) =1: Choose x with prob_sol, set b=Ax
! abs(test_type)!=1: Choose b with prob_sol
subroutine main_driver()

  use ml_layout_module
  use define_bc_module
  use stag_mg_solver_module
  use macproject_module
  use write_plotfile_module
  use convert_stag_module
  use apply_precon_module
  use apply_matrix_module
  use norm_inner_product_module
  use init_module
  use multifab_physbc_stag_module
  use convert_to_homogeneous_module
  use ParallelRNGs
  use bc_module
  use stag_applyop_module
  use gmres_module
  use probin_module       , only: probin_init, test_type, theta_alpha_fac
  use probin_common_module, only: probin_common_init, n_cells, dim_in, fixed_dt, &
                                  max_grid_size, plot_int, seed, &
                                  prob_lo, prob_hi, bc_lo, bc_hi
  use probin_gmres_module , only: probin_gmres_init, scale_factor

  implicit none

  ! will be allocated with dm components
  integer       , allocatable :: lo(:), hi(:)

  ! will be allocated with (nlevs,dm) components
  real(dp_t)    , allocatable :: dx(:,:)

  ! staggered multifabs
  ! will be allocated with (nlevs,dm) components
  type(multifab), allocatable :: umac_exact(:,:) ! temporary storage for exact solution
  type(multifab), allocatable :: umac(:,:)       ! velocity increment
  type(multifab), allocatable :: umac_tmp(:,:)   ! temporary multifab
  type(multifab), allocatable :: rhs_u(:,:)      ! right-hand-side for velocity terms
  type(multifab), allocatable :: rhs_p(:)        ! right-hand-side for pressure term
  type(multifab), allocatable :: grad_pres(:,:)  ! pressure gradient

  ! cell-centered multifabs
  ! will be allocated with nlevs components
  type(multifab), allocatable :: pres_exact(:)    ! temporary storage for exact solution
  type(multifab), allocatable :: pres(:)          ! pressure increment
  type(multifab), allocatable :: pres_tmp(:)      ! temporary multifab
  type(multifab), allocatable :: alpha(:)         ! coefficient for staggered multigrid solver
  type(multifab), allocatable :: alpha_fc(:,:)    ! coefficient for staggered multigrid solver
  type(multifab), allocatable :: beta(:)          ! coefficient for staggered multigrid solver
  type(multifab), allocatable :: beta_ed(:,:)     ! nodal (2d), edge-centered (3d)
  type(multifab), allocatable :: gamma(:)         ! coefficient for staggered multigrid solver

  ! special inhomogeneous boundary condition multifab
  ! vel_bc_n(nlevs,dm) are the normal velocities
  ! in 2D, vel_bc_t(nlevs,2) respresents
  !   1. y-velocity bc on x-faces (nodal)
  !   2. x-velocity bc on y-faces (nodal)
  ! in 3D, vel_bc_t(nlevs,6) represents
  !   1. y-velocity bc on x-faces (nodal in y and x)
  !   2. z-velocity bc on x-faces (nodal in z and x)
  !   3. x-velocity bc on y-faces (nodal in x and y)
  !   4. z-velocity bc on y-faces (nodal in z and y)
  !   5. x-velocity bc on z-faces (nodal in x and z)
  !   6. y-velocity bc on z-faces (nodal in y and z)
  type(multifab), allocatable :: vel_bc_n(:,:)
  type(multifab), allocatable :: vel_bc_t(:,:)

  integer    :: dm,nlevs,n,i
  real(dp_t) :: time,theta_alpha
  
  type(box)         :: bx
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla

  type(bc_tower) :: the_bc_tower

  real(kind=dp_t) :: norm, norm_scale, vol
  real(kind=dp_t), allocatable :: norm_stag(:), norm_u_diff(:)

  real(kind=dp_t) :: mean_val_pres
  real(kind=dp_t), allocatable ::  mean_val_umac(:)

  logical, allocatable :: pmask(:)
  logical :: nodal_temp(3)

  call probin_init()
  call probin_common_init()
  call probin_gmres_init()

  ! Initialize random numbers *after* the global (root) seed has been set:
  call SeedParallelRNG(seed)

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  dm = dim_in

  theta_alpha = theta_alpha_fac/fixed_dt

  time = 0.d0

  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(norm_stag(dm))
  allocate(norm_u_diff(dm))
  allocate(mean_val_umac(dm))

  ! now that we have nlevs and dm, we can allocate these
  allocate(dx(nlevs,dm))
  allocate(umac_exact(nlevs,dm),umac(nlevs,dm),umac_tmp(nlevs,dm),vel_bc_n(nlevs,dm))
  allocate(rhs_u(nlevs,dm),rhs_p(nlevs),grad_pres(nlevs,dm))
  allocate(pres_exact(nlevs),pres(nlevs),pres_tmp(nlevs))
  allocate(alpha(nlevs),beta(nlevs),gamma(nlevs))
  allocate(alpha_fc(nlevs,dm))
  if (dm .eq. 2) then
     allocate(beta_ed(nlevs,1))
     allocate(vel_bc_t(nlevs,2))
  else if (dm .eq. 3) then
     allocate(beta_ed(nlevs,3))
     allocate(vel_bc_t(nlevs,6))
  end if

  ! tell mba how many levels and dmensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dm)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  enddo

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  dx(1,1:dm) = (prob_hi(1)-prob_lo(1)) / n_cells(1:dm)
  select case (dm) 
    case(2)
      if (dx(1,1) .ne. dx(1,2)) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy')
      end if    
    case(3)
      if ((dx(1,1) .ne. dx(1,2)) .or. (dx(1,1) .ne. dx(1,3))) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy=dz')
      end if    
    case default
      call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end select
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! create a box from (0,0) to (n_cells-1,n_cells-1)
  lo(1:dm) = 0
  hi(1:dm) = n_cells(1:dm)-1
  bx = make_box(lo,hi)

  ! tell mba about the problem domain at each level
  mba%pd(1) = bx
  do n=2,nlevs
     mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
  enddo

  ! initialize the boxarray at level 1 to be one single box
  call boxarray_build_bx(mba%bas(1),bx)

  ! overwrite the boxarray at level 1 to respect max_grid_size
  call boxarray_maxsize(mba%bas(1),max_grid_size)

  ! now build the boxarray at other levels
  if (nlevs .ge. 2) then
     call bl_error("Need to build boxarray for n>1")
  end if

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  ! build the ml_layout, mla
  call ml_layout_build(mla,mba,pmask)

  deallocate(pmask)

  ! don't need this anymore - free up memory
  call destroy(mba)

  ! tell the_bc_tower about max_levs, dm, and domain_phys_bc
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask,1,1)
  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  do n=1,nlevs
     do i=1,dm
        call multifab_build_edge(umac_exact(n,i) ,mla%la(n),1,1,i)
        call multifab_build_edge(umac(n,i)       ,mla%la(n),1,1,i)
        call multifab_build_edge(umac_tmp(n,i)   ,mla%la(n),1,1,i)
        call multifab_build_edge(rhs_u(n,i)      ,mla%la(n),1,0,i)
        call multifab_build_edge(grad_pres(n,i)  ,mla%la(n),1,0,i)
     end do
     call multifab_build(rhs_p(n)         ,mla%la(n),1,0)
     call multifab_build(pres_exact(n)    ,mla%la(n),1,1)
     call multifab_build(pres(n)          ,mla%la(n),1,1)
     call multifab_build(pres_tmp(n)      ,mla%la(n),1,1)
     call multifab_build(alpha(n)         ,mla%la(n),1,1)
     call multifab_build(beta(n)          ,mla%la(n),1,1)
     call multifab_build(gamma(n)         ,mla%la(n),1,1)

     do i=1,dm
        call multifab_build_edge(alpha_fc(n,i),mla%la(n),1,0,i)
     end do

     if (dm .eq. 2) then
        call multifab_build_nodal(beta_ed(n,1),mla%la(n),1,0)
     else
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(beta_ed(n,1),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(beta_ed(n,2),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(beta_ed(n,3),mla%la(n),1,0,nodal_temp)
     end if

     ! boundary conditions
     do i=1,dm
        call multifab_build_edge(vel_bc_n(n,i),mla%la(n),1,0,i)
     end do
     if (dm .eq. 2) then
        ! y-velocity bc on x-faces (nodal)
        call multifab_build_nodal(vel_bc_t(n,1),mla%la(n),1,0)
        ! x-velocity bc on y-faces (nodal)
        call multifab_build_nodal(vel_bc_t(n,2),mla%la(n),1,0)
     else
        ! y-velocity bc on x-faces (nodal in y and x)
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
        ! z-velocity bc on x-faces (nodal in z and x)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
        ! x-velocity bc on y-faces (nodal in x and y)
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
        ! z-velocity bc on y-faces (nodal in z and y)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
        ! x-velocity bc on z-faces (nodal in x and z)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
        ! y-velocity bc on z-faces (nodal in y and z)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
     end if

  end do

  ! set inhomogeneous bc condition
  call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,beta_ed,dx,0.d0,the_bc_tower%bc_tower_array)

  ! provide an initial value (or guess) for umac and pres 
  call init_solution(mla,umac_exact,pres_exact,dx,time,the_bc_tower%bc_tower_array, &
                     vel_bc_n,vel_bc_t)

  ! for testing purposes, ensure mean zero for pressure (should not be required any more)
  if (.false.) then        
     do n=1,nlevs
         mean_val_pres = multifab_sum_c(pres_exact(n),1,1,all=.false.) / &
                         multifab_volume(pres_exact(n))
         call multifab_sub_sub_s(pres_exact(n),mean_val_pres,pres_exact(n)%ng) ! Ensure mean zero            
     end do   
  end if
  
  ! initialize alpha, beta, and gamma
  call init_mat(mla,alpha,beta,gamma,dx,time,the_bc_tower%bc_tower_array)

  ! compute alpha_fc and beta_ed
  call average_cc_to_face(nlevs,alpha,alpha_fc,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  if (dm .eq. 2) then
     call average_cc_to_node(nlevs,beta,beta_ed(:,1),1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  else if (dm .eq. 3) then
     call average_cc_to_edge(nlevs,beta,beta_ed,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  end if

  if (abs(test_type) .eq. 1) then
     ! initialize rhs_u and rhs_p by explicitly computing rhs = A x

     ! try scaling the pressure to make the solution well-scaled
     call multifab_mult_mult_s(pres_exact(1),1.0d0/scale_factor,pres_exact(1)%ng)

     call apply_matrix(mla,rhs_u,rhs_p,umac_exact,pres_exact,alpha_fc,beta,beta_ed, &
                       gamma,theta_alpha,dx,the_bc_tower,vel_bc_n,vel_bc_t)
  else
     ! initialize rhs_u and rhs_p with a subroutine
     call init_rhs(mla,rhs_u,rhs_p,dx,time,the_bc_tower%bc_tower_array)
  end if 

  if (plot_int .gt. 0) then
     ! write a plotfile of umac_exact and pres_exact
    call write_plotfile(mla,umac_exact,pres_exact,alpha,beta,gamma,dx,time,0)
  end if

TestType: if (test_type==0) then ! Test the order of accuracy of the stencils
    ! Here we keep the inhomogeneous form of the BCs to test them
    
    ! calculate A*x and save it in umac_tmp
    call apply_matrix(mla,umac_tmp,pres_tmp,umac_exact,pres_exact,alpha_fc, &
                      beta,beta_ed,gamma,theta_alpha,dx,the_bc_tower,vel_bc_n,vel_bc_t)
    
    ! calculate f - (A*u)
    do n=1,nlevs
      do i=1,dm
        call multifab_sub_sub_c(rhs_u(n,i),1,umac_tmp(n,i),1,1,0)
      end do
      call multifab_sub_sub_c(rhs_p(n),1,pres_tmp(n),1,1,0)
    end do 

    if (plot_int .gt. 0) then
      ! write a plotfile of rhs_u 
      call write_plotfile(mla,rhs_u,rhs_p,alpha,beta,gamma,dx,time,100)
    end if

    ! calculate the norms 
    norm = multifab_norm_inf_c(rhs_u(1,1),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"The discretization error: L0 U   =",norm
    norm = multifab_norm_inf_c(rhs_u(1,2),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"The discretization error: L0 V   =",norm
    if (dm .eq. 3) then 
      norm = multifab_norm_inf_c(rhs_u(1,3),1,1,all=.false.)
      if (parallel_IOProcessor()) print*,"The discretization error: L0 W   =",norm
    end if 

    norm = multifab_norm_inf_c(rhs_p(1),1,1,all=.false.)
    if (parallel_IOProcessor()) print*,"The discretization error: L0 P   =",norm

    ! Dirichlet BC case, we measure the global error
    if ((bc_lo(1) .eq. 17) .or. (bc_lo(2) .eq. 17) .or. (bc_hi(1) .eq. 17) .or. (bc_hi(2) .eq. 17))  then 
      ! rhs_u is the local truncation error, we put it as rhs and solve for the global error
      do n=1,nlevs
        call setval(rhs_p(n),0.d0,all=.true.)
      end do 
      call gmres(mla,the_bc_tower,dx,rhs_u,rhs_p,umac,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha)
      ! calculate the norms of the global error
      norm = multifab_norm_inf_c(umac(1,1),1,1,all=.false.)
      if (parallel_IOProcessor()) print*,"The global error: L0 U   =",norm
      norm = multifab_norm_inf_c(umac(1,2),1,1,all=.false.)
      if (parallel_IOProcessor()) print*,"The global error: L0 V   =",norm
      if (dm .eq. 3) then 
        norm = multifab_norm_inf_c(umac(1,3),1,1,all=.false.)
        if (parallel_IOProcessor()) print*,"The global error: L0 W   =",norm
      end if 
    end if 

else TestType ! Actually try to solve the linear system by gmres or pure multigrid

  ! for inhomogeneous boundary conditions, convert problem to homogeneous by
  ! subtracting from the RHS the result of the operator applied to a solution
  ! vector with zeros everywhere in the problem domain, and ghost cells filled to
  ! respect the boundary conditions
  call convert_to_homogeneous(mla,rhs_u,rhs_p,alpha_fc,beta,beta_ed, &
                              gamma,theta_alpha,dx,the_bc_tower,vel_bc_n,vel_bc_t)

  ! compute the average value of umac_exact and pres_exact
  call sum_umac_press(mla,pres_exact,umac_exact,mean_val_pres,mean_val_umac) 
  ! And now set the initial guess to the average value (cannot hurt when there is no null space)
  do n=1,nlevs
     do i=1,dm
        call setval(umac(n,i),mean_val_umac(i),all=.true.)
     end do
     call setval(pres(n),mean_val_pres,all=.true.)
  end do

  call sum_umac_press(mla,pres,umac,mean_val_pres,mean_val_umac) ! Check the initial sum
  if (parallel_IOProcessor()) then 
     write(*,*) "CONSERVE: Initial: <u>=", mean_val_umac, " <p>=", mean_val_pres
  end if
   
  if (plot_int .gt. 0) then
     ! write a plotfile of umac and pres
     call write_plotfile(mla,umac,pres,alpha,beta,gamma,dx,time,1)
  end if
      
  if (test_type>0) then
  
     ! use gmres
     call gmres(mla,the_bc_tower,dx,rhs_u,rhs_p,umac,pres,alpha_fc,beta,beta_ed,gamma,theta_alpha)
  
  else
  
     ! use a single application of the preconditioner
     ! the viscous solver uses 0 as an initial guess so we need to solve in residual
     ! correction form.

     ! Since we want to solve A*x=b with an initial guess x0, we instead solve
     ! A*dx=(b-A*x0) and set x=x0+dx.

     call apply_matrix(mla,umac_tmp,pres_tmp,umac,pres,alpha_fc, &
                       beta,beta_ed,gamma,theta_alpha,dx,the_bc_tower)

     do n=1,nlevs
        do i=1,dm
           call multifab_sub_sub_c(rhs_u(n,i),1,umac_tmp(n,i),1,1,0)
        end do
           call multifab_sub_sub_c(rhs_p(n),1,pres_tmp(n),1,1,0)
     end do

     call apply_precon(mla,rhs_u,rhs_p,umac_tmp,pres_tmp,alpha_fc,beta,beta_ed,gamma, &
                       theta_alpha,dx,the_bc_tower)

     do n=1,nlevs
        do i=1,dm
           call multifab_plus_plus_c(umac(n,i),1,umac_tmp(n,i),1,1,1)
        end do
        call multifab_plus_plus_c(pres(n),1,pres_tmp(n),1,1,1)
     end do

  end if

  call sum_umac_press(mla,pres,umac,mean_val_pres,mean_val_umac) ! Check the final sum
  if (parallel_IOProcessor()) then 
     write(*,*) "CONSERVE: Final: <u>=", mean_val_umac, " <p>=", mean_val_pres
  end if
  
  ! Now we need to convert back to non-homogeneous form by correcting velocities at the boundaries
  ! There is no need to fill ghost cells here
  do n=1,nlevs
     do i=1,dm
        call multifab_physbc_domainvel(umac(n,i),i,the_bc_tower%bc_tower_array(n), &
                                       dx(n,:),vel_bc_n(n,:))
     end do
  end do

  if (plot_int .gt. 0) then
     ! write a plotfile of umac and pres
     call write_plotfile(mla,umac,pres,alpha,beta,gamma,dx,time,2)
  end if

  ! Due to scaling, this may be very different from 1
  norm_scale = multifab_norm_inf_c(pres_exact(1),1,1,all=.false.)
  
  ! for checking the error, umac_exact-umac
  do n=1,nlevs  
     do i=1,dm 
        call multifab_sub_sub_c(umac_exact(n,i),1,umac(n,i),1,1,0)
     end do
     call multifab_sub_sub_c(pres_exact(n),1,pres(n),1,1,0)
  end do

  ! check that the error is almost zero
  if (plot_int .gt. 0) then
     ! write a plotfile of error
     call write_plotfile(mla,umac_exact,pres_exact,alpha,beta,gamma,dx,time,3)
  end if

  vol = multifab_volume(umac_exact(1,1))

  if (parallel_IOProcessor()) print*,""

  ! compute L0 norms for staggered data
  norm = multifab_norm_inf_c(umac_exact(1,1),1,1,all=.false.)
  if (parallel_IOProcessor()) print*,"L0 U   =",norm
  norm = multifab_norm_inf_c(umac_exact(1,2),1,1,all=.false.)
  if (parallel_IOProcessor()) print*,"L0 V   =",norm
  if (dm .eq. 3) then
     norm = multifab_norm_inf_c(umac_exact(1,3),1,1,all=.false.)
     if (parallel_IOProcessor()) print*,"L0 W   =",norm
  end if
  norm = multifab_norm_inf_c(pres_exact(1),1,1,all=.false.)
  if (parallel_IOProcessor()) print*,"L0 P   =",norm
  if (parallel_IOProcessor()) print*,"L0 P normalized   =", norm/norm_scale
  if (parallel_IOProcessor()) print*,""

  ! compute L1 norms for staggered data
  call stag_l1_norm(mla,umac_exact,norm_stag)
  if (parallel_IOProcessor()) print*,"L1 U   =",norm_stag(1) / dble(vol)
  if (parallel_IOProcessor()) print*,"L1 V   =",norm_stag(2) / dble(vol)
  if (dm .eq. 3) then
     if (parallel_IOProcessor()) print*,"L1 W   =",norm_stag(3) / dble(vol)
  end if
  norm = multifab_norm_l1_c(pres_exact(1),1,1,all=.false.)
  if (parallel_IOProcessor()) print*,"L1 P   =",norm / dble(vol)
  if (parallel_IOProcessor()) print*,""

  ! compute L2 norms for staggered data
  call stag_inner_prod(mla,umac_exact,1,umac_exact,1,norm_stag)
  if (parallel_IOProcessor()) print*,"L2 U   =",sqrt(norm_stag(1)/dble(vol))
  if (parallel_IOProcessor()) print*,"L2 V   =",sqrt(norm_stag(2)/dble(vol))
  if (dm .eq. 3) then
     if (parallel_IOProcessor()) print*,"L2 W   =",sqrt(norm_stag(3)/dble(vol))
  end if
  call cc_inner_prod(mla,pres_exact,1,pres_exact,1,norm)
  if (parallel_IOProcessor()) print*,"L2 P   =",sqrt(norm/dble(vol))
  if (parallel_IOProcessor()) print*,""

end if TestType

  ! free memory
  do n=1,nlevs
     do i=1,dm
        call multifab_destroy(umac_exact(n,i))
        call multifab_destroy(umac(n,i))
        call multifab_destroy(umac_tmp(n,i))
        call multifab_destroy(rhs_u(n,i))
        call multifab_destroy(grad_pres(n,i))
        call multifab_destroy(alpha_fc(n,i))
        call multifab_destroy(vel_bc_n(n,i))
     end do
     call multifab_destroy(rhs_p(n))
     call multifab_destroy(pres_exact(n))
     call multifab_destroy(pres(n))
     call multifab_destroy(pres_tmp(n))
     call multifab_destroy(alpha(n))
     call multifab_destroy(beta(n))
     call multifab_destroy(gamma(n))
     if (dm .eq. 2) then
        call multifab_destroy(beta_ed(n,1))
        call multifab_destroy(vel_bc_t(n,1))
        call multifab_destroy(vel_bc_t(n,2))
     else if (dm .eq. 3) then
        call multifab_destroy(beta_ed(n,1))
        call multifab_destroy(beta_ed(n,2))
        call multifab_destroy(beta_ed(n,3))
        call multifab_destroy(vel_bc_t(n,1))
        call multifab_destroy(vel_bc_t(n,2))
        call multifab_destroy(vel_bc_t(n,3))
        call multifab_destroy(vel_bc_t(n,4))
        call multifab_destroy(vel_bc_t(n,5))
        call multifab_destroy(vel_bc_t(n,6))
     end if
  end do
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)
  deallocate(lo,hi,norm_stag,dx)
  deallocate(umac_exact,umac,umac_tmp,rhs_u,grad_pres)
  deallocate(pres_exact,pres,pres_tmp,alpha,beta,gamma)
  deallocate(mean_val_umac)
  deallocate(alpha_fc,beta_ed)

end subroutine main_driver
