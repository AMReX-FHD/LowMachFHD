subroutine main_driver()

  use boxlib
  use bl_IO_module
  use ml_layout_module
  use init_lowmach_module
  use init_temp_module
  use compute_mixture_properties_module
  use initial_projection_module
  use write_plotfileLM_module
  use advance_timestep_overdamped_module
  use advance_timestep_inertial_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use analysis_module
  use analyze_spectra_module
  use div_and_grad_module
  use eos_check_module
  use estdt_module
  use stag_mg_layout_module
  use macproject_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use fill_umac_ghost_cells_module
  use fill_rho_ghost_cells_module
  use ParallelRNGs 
  use mass_flux_utilities_module
  use compute_HSE_pres_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_m_to_umac_module
  use sum_momenta_module
  use checkpoint_module
  use project_onto_eos_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, hydro_grid_int, &
                                  max_grid_size, n_steps_save_stats, n_steps_skip, &
                                  plot_int, chk_int, seed, stats_int, bc_lo, bc_hi, restart, &
                                  probin_common_init, print_int, project_eos_int, &
                                  fixed_dt, max_step, cfl
  use probin_multispecies_module, only: nspecies, Dbar, &
                                        start_time, &
                                        probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init

  use fabio_module

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2,Dbar_max,dt_diffusive
  integer                      :: n,nlevs,i,dm,istep,ng_s,init_step,n_Dbar
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: n_old(:)
  type(multifab), allocatable  :: n_new(:)
  type(multifab), allocatable  :: diff_fluxdiv(:)
  type(multifab), allocatable  :: stoch_fluxdiv(:)

  ! For HydroGrid
  integer :: narg, farg, un, n_rngs
  character(len=128) :: fname
  logical :: lexist
  logical :: nodal_temp(3)
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  
  ! Initialize random numbers *after* the global (root) seed has been set:
  call SeedParallelRNG(seed)

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(n_old(nlevs),n_new(nlevs))
  allocate(diff_fluxdiv(nlevs),stoch_fluxdiv(nlevs))

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  allocate(dx(nlevs,dm))
  dx(1,1:dm) = (prob_hi(1:dm)-prob_lo(1:dm)) / n_cells(1:dm)
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

  ! use refined dx for next level
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  ng_s = 2

  if (restart .ge. 0) then

     call bl_error("restart not supported yet")

  else

     init_step = 1
     time = start_time
     
     ! tell mba how many levels and dimensionality of problem
     call ml_boxarray_build_n(mba,nlevs,dm)

     ! tell mba about the ref_ratio between levels
     ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
     ! we use refinement ratio of 2 in every direction between all levels
     do n=2,nlevs
        mba%rr(n-1,:) = 2
     end do

     ! create a box from (0,0) to (n_cells-1,n_cells-1)
     lo(1:dm) = 0
     hi(1:dm) = n_cells(1:dm)-1
     bx = make_box(lo,hi)

     ! tell mba about the problem domain at every level
     mba%pd(1) = bx
     do n=2,nlevs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     end do

     ! initialize the boxarray at level 1 to be one single box
     call boxarray_build_bx(mba%bas(1),bx)

     ! overwrite the boxarray at level 1 to respect max_grid_size
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! now build the boxarray at other levels
     if (nlevs .ge. 2) then
        call bl_error("Need to build boxarray for n>1")
     end if

     ! build the ml_layout, mla
     call ml_layout_build(mla,mba,pmask)

     ! don't need this anymore - free up memory
     call destroy(mba)

     do n=1,nlevs
        call multifab_build(n_new(n)        ,mla%la(n),nspecies,ng_s)
        call multifab_build(n_old(n)        ,mla%la(n),nspecies,ng_s) 
        call multifab_build(diff_fluxdiv(n) ,mla%la(n),nspecies,0) 
        call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0) 
     end do

  end if

  deallocate(pmask)

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  ! dm+2 = scal_bc_comp = rhotot
  ! scal_bc_comp+1 = c_i
  ! scal_bc_comp+nspecies+1 = molfrac or massfrac (dimensionless fractions)
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = tran_bc_comp = diffusion coefficients (eta,kappa,chi)
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies+2, &
                     num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! these quantities are populated here and defined in bc.f90
  c_bc_comp   = scal_bc_comp + 1
  mol_frac_bc_comp   = scal_bc_comp + nspecies + 1
  temp_bc_comp       = scal_bc_comp + 2*nspecies + 1

  ! allocate and build multifabs that will contain random numbers
  n_rngs = 1
  call init_mass_stochastic(mla,n_rngs)

  !=====================================================================
  ! Initialize values
  !=====================================================================

  if (restart .lt. 0) then

     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call bl_error("Need to define what we mean by CFL time step")
     end if
     
  end if

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     runtime1 = parallel_wtime()

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! advance the solution by dt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! fill random flux multifabs with new random numbers
     call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

     ! compute diffusive flux divergence

     ! compute stochastic flux divergence

     ! n_k^{n+1} = n_k^n + (dt/2)(div D_k grad n_k)^n
     !                   + (dt/2)(div D_k grad n_k)^n+1
     !                   +  dt    div (sqrt(2 D_k n_k) Z)

     ! form RHS for implicit system

     ! solve the implicit system


      time = time + dt

      if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
      end if

      runtime2 = parallel_wtime() - runtime1
      call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
      if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,'Time to advance timestep: ',runtime1,' seconds'
      end if

      ! set old state to new state
      do n=1,nlevs
         call multifab_copy_c(n_old(n),1,n_new(n),1,nspecies,n_old(n)%ng)
      end do

  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  call destroy_mass_stochastic(mla)

  do n=1,nlevs
  end do

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
