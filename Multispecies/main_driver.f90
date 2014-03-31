subroutine main_driver()

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use init_module
  use write_plotfile_module
  use advance_module
  use define_bc_module
  use bc_module
  use analysis_module
  use ParallelRNGs 
  use convert_variables_module
  use probin_common_module
  use probin_multispecies_module
 
  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time
  integer                      :: n,nlevs,i,j,dm,istep,step_count
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: rho(:)
  type(multifab), allocatable  :: rho_exact(:)
  real(kind=dp_t),allocatable  :: molmass(:) 
  real(kind=dp_t),allocatable  :: covW(:,:) 
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 

  ! for time being, we fix nlevs to be 1. for adaptive simulations where the grids 
  ! change, cells at finer resolution don't necessarily exist depending on your 
  ! tagging criteria, so max_levs isn't necessary equal to nlevs
  nlevs = 1
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(dx(nlevs,dm))
  allocate(rho(nlevs))
  allocate(rho_exact(nlevs))
  allocate(molmass(nspecies))
  allocate(covW(nspecies,nspecies))

  !==============================================================
  ! Setup parallelization: Create boxes and layouts for multifabs
  !==============================================================
  
  ! tell mba how many levels and dimensionality of problem
  call ml_boxarray_build_n(mba,nlevs,dm)

  ! tell mba about the ref_ratio between levels
  ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
  ! we use refinement ratio of 2 in every direction between all levels
  do n=2,nlevs
     mba%rr(n-1,:) = 2
  end do

  ! set grid spacing at each level; presently grid spacing is same in all direction
  dx(1,1:dm) = (prob_hi(1)-prob_lo(1)) / n_cells(1:dm)
  
  ! check whether dimensionality & grid spacing passed correct 
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

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  ! dm+2 = scal_bc_comp = rho_tot
  ! scal_bc_comp+1 = rho_i
  ! scal_bc_comp+nspecies+1 = mol_frac
  ! scal_bc_comp+2*nspecies+1 = tran_bc_comp = diff_coef
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies,num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! these quantities are populated here and defined in probin_multispecies 
  rho_part_bc_comp   = scal_bc_comp + 1
  mol_frac_bc_comp   = scal_bc_comp + nspecies + 1
  diff_coeff_bc_comp = tran_bc_comp

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  ! build multifab with nspecies component and one ghost cell
  do n=1,nlevs
     call multifab_build(rho(n),      mla%la(n),nspecies,1)
     call multifab_build(rho_exact(n),mla%la(n),nspecies,1)
  end do

  ! Initialize random numbers *after* the global (root) seed has been set:
  if(use_stoch) call SeedParallelRNG(seed)

  !=====================================================================
  ! Read molar mass from input file (constant throughout space and time)
  !=====================================================================
  molmass(1:nspecies) = molmass_in(1:nspecies)  

  ! initialize the time 
  time = start_time    

  ! initialize rho
  call init_rho(rho,dx,prob_lo,prob_hi,time,the_bc_tower%bc_tower_array)

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  istep = 0

  ! write initial plotfile
  if (plot_int .gt. 0) then
     call write_plotfile(mla,"plt_rho",rho,istep,dx,time,prob_lo,prob_hi)
  end if

  ! print out the total masses
  call sum_mass(rho, step=0)
 
  ! choice of time step with a diffusive CFL of 0.1; CFL=minimum[dx^2/(2*chi)]; 
  ! chi is the largest eigenvalue of diffusion matrix to be input for n-species
  dt = cfl1*dx(1,1)**2/chi
  if (parallel_IOProcessor()) then
     write(*,*) "Using time step dt=", dt
  end if

  ! set the time counter for the time-average
  step_count = 0
  covW       = 0 

  do istep=1,max_step

     if (parallel_IOProcessor()) then
        !print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     ! advance the solution by dt
     call advance(mla,rho,molmass,dx,dt,time,prob_lo,prob_hi,the_bc_tower%bc_tower_array)

     ! print out the total mass to check conservation
     !call sum_mass(rho, istep)

     ! compute error norms
     if (print_error_norms) then
        call print_errors(rho,rho_exact,dx,prob_lo,prob_hi,time,the_bc_tower%bc_tower_array)
     end if

     ! compute coavariance and variances 
     if(max_step .gt. 1) then
        call compute_cov(mla,rho,covW)    
        !call meanvar_W(mla,rho,covW)   
        step_count = step_count + 1 
     end if 

     ! write plotfile at specific intervals
     if ((plot_int.gt.0 .and. mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) then
        
        call write_plotfile(mla,"plt_rho",      rho,istep,dx,time,prob_lo,prob_hi)
        call write_plotfile(mla,"plt_exa",rho_exact,istep,dx,time,prob_lo,prob_hi)

        ! difference between rho and rho_exact
        do n=1,nlevs
           call saxpy(rho_exact(n),-1.0d0,rho(n))
        end do
        
        ! check error with visit
        call write_plotfile(mla,"plt_err",rho_exact,istep,dx,time,prob_lo,prob_hi)
 
     end if
     
     ! increment simulation time
     time = time + dt
        
  end do

  ! print out the total mass to check conservation
  call sum_mass(rho, istep)
 
  ! print out the standard deviation
  if (parallel_IOProcessor()) then
     do i=1,nspecies
        do j=1,nspecies
           print*, ' cov of Wij',i,j,covW(i,j)/step_count
        end do
     end do
           
     print*, ' analytic std of W =', (molmass(1)*rho_in(1,2) + molmass(2)*rho_in(1,1))*&
                                     rho_in(1,1)*rho_in(1,2)/(product(dx(1,1:dm))*thickness*&
                                     (rho_in(1,1)+rho_in(1,2))**4)
  end if

 
  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  deallocate(molmass)
  deallocate(covW)
  do n=1,nlevs
     call multifab_destroy(rho(n))
     call multifab_destroy(rho_exact(n))
  end do
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
