subroutine main_driver()

   use boxlib
   use bl_IO_module
   use bl_space
   use ml_layout_module
   use define_bc_module
   use bc_module
   use init_n_module
   use write_plotfile_module
   use stochastic_n_fluxdiv_module
   use advance_timestep_module
   use analyze_spectra_module
   use compute_n_steady_module
   use restart_module
   use checkpoint_module
   use ParallelRNGs 
   use bl_rng_module
   use multifab_physbc_module
   use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, max_grid_size, &
                                   plot_int, chk_int, print_int, seed, bc_lo, bc_hi, &
                                   restart, use_bl_rng, nspecies, &
                                   probin_common_init, fixed_dt, max_step, n_steps_skip, &
                                   hydro_grid_int, stats_int, n_steps_save_stats, &
                                   cfl, initial_variance_mass
   use probin_reactdiff_module, only: probin_reactdiff_init, D_Fick, &
                                      inhomogeneous_bc_fix, temporal_integrator, &
                                      n_steps_write_avg, &
                                      model_file_init, model_file, integer_populations, &
                                      volume_factor
   use probin_chemistry_module, only: probin_chemistry_init, nreactions, include_discrete_LMA_correction, &
                                      exclude_solvent_comput_rates, use_mole_frac_LMA

   use fabio_module

   implicit none

   ! quantities will be allocated with dm components
   integer, allocatable :: lo(:), hi(:)

   ! nlevs is always 1 and so it is redundant but to be consistent with BoxLib we need 
   ! to carry it around.
   ! quantities will be allocated with (nlevs,dm) components
   real(kind=dp_t), allocatable :: dx(:,:)
   real(kind=dp_t)              :: dt,time,runtime1,runtime2,cellvolume
   integer                      :: n,nlevs,i,j,k,dm,istep,ng_s,init_step,spec
   type(box)                    :: bx
   type(ml_boxarray)            :: mba
   type(ml_layout)              :: mla
   type(bc_tower)               :: the_bc_tower
   logical, allocatable         :: pmask(:)

   type(multifab), allocatable :: n_old(:)
   type(multifab), allocatable :: n_new(:)
   type(multifab), allocatable :: n_steady(:)

   integer :: n_rngs

   ! For HydroGrid
   integer :: narg, farg, un
   character(len=128) :: fname
   logical :: lexist

   ! to test "conservation"
   real(kind=dp_t), allocatable :: n_sum(:)

   real(kind=dp_t), allocatable :: input_array(:,:,:)

   !==============================================================
   ! Initialization
   !==============================================================

   call probin_common_init()
   call probin_reactdiff_init() 
   call probin_chemistry_init()

   if (use_bl_rng) then
      ! Build the random number engine and give initial distributions for the
      ! F_BaseLib/bl_random RNG module
      call rng_init()
   else
      ! Initialize random numbers *after* the global (root) seed has been set:
      ! This is for the RNG module that sits in Hydrogrid
      call SeedParallelRNG(seed)
   end if

   if (nreactions > 0 .and. use_mole_frac_LMA) then
      if (include_discrete_LMA_correction) then
         call bl_error('Error: currently use_mole_frac_LMA can be used only with include_discrete_LMA_correction=F')
      end if
      if (exclude_solvent_comput_rates .ne. 0) then
         call bl_error('Error: currently use_mole_frac_LMA can be used only with exclude_solvent_comput_rates=0')
      end if
   end if

   ! in this example we fix nlevs to be 1
   ! for adaptive simulations where the grids change, cells at finer
   ! resolution don't necessarily exist depending on your tagging criteria,
   ! so max_levs isn't necessary equal to nlevs
   nlevs = 1

   ! dimensionality is set in inputs file
   dm = dim_in

   ! now that we have dm, we can allocate these
   allocate(lo(dm),hi(dm))
   allocate(n_old(nlevs),n_new(nlevs),n_steady(nlevs))

   allocate(n_sum(nspecies))

   ! build pmask
   allocate(pmask(dm))
   pmask = .false.
   do i=1,dm
      if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
         pmask(i) = .true.
      end if
   end do

   ng_s = 1 

   if (restart .ge. 0) then

      init_step = restart + 1

     ! build the ml_layout
     ! read in time and dt from checkpoint
     ! build and fill n_old
     call initialize_from_restart(mla,time,dt,n_old,pmask,ng_s)

   else

      init_step = 1
      time = 0.d0

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
         call multifab_build(n_old(n),mla%la(n),nspecies,ng_s) 
      end do

   end if

   do n=1,nlevs
      call multifab_build(n_new(n),mla%la(n),nspecies,ng_s)
      if (inhomogeneous_bc_fix .and. temporal_integrator .ge. 0) then
         call multifab_build(n_steady(n),mla%la(n),nspecies,0)
      end if
   end do

   deallocate(pmask)

   ! allocate and build multifabs that will contain random numbers
   ! in this case, we only use one rng because we are using an Ito interpretation
   ! so that we can compute the term involving W_1 first, and then compute the one involving W_2 only
   n_rngs = 1
   call init_mass_stochastic(mla,n_rngs)

  ! set grid spacing at each level
  allocate(dx(nlevs,MAX_SPACEDIM))
  dx(1,1:MAX_SPACEDIM) = (prob_hi(1:MAX_SPACEDIM)-prob_lo(1:MAX_SPACEDIM)) &
       / n_cells(1:MAX_SPACEDIM)
  ! check that the grid spacing is the same in each direction for the first dm dimensions
  select case (dm) 
  case(2)
     if (dx(1,1) .ne. dx(1,2)) then
        call bl_error('ERROR: main_driver.f90, in 2D we only support dx=dy')
     end if
  case(3)
     if ((dx(1,1) .ne. dx(1,2)) .or. (dx(1,1) .ne. dx(1,3))) then
        call bl_error('ERROR: main_driver.f90, in 3D we only support dx=dy=dz')
     end if
  case default
     call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end select

  ! use refined dx for next level
  ! assume refinement ratio is the same in each direction
  ! we do this because dx is allocated over MAX_SPACEDIM dimensions,
  ! whereas mba is only allocated over dm dimensions
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,1)
  end do

   !=======================================================
   ! Setup boundary condition bc_tower
   !=======================================================

   ! bc_tower structure in memory
   ! 1:dm = velocity
   ! dm+1 = pressure
   ! dm+2 = scal_bc_comp (for n_i)
   ! scal_bc_comp+nspecies = tran_bc_comp = diffusion coefficients
   call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                      num_scal_bc_in=nspecies, &
                      num_tran_bc_in=1)

   do n=1,nlevs
      ! define level n of the_bc_tower
      call bc_tower_level_build(the_bc_tower,n,mla%la(n))
   end do

   !=====================================================================
   ! Initialize values
   !=====================================================================

   if (restart .ge. 0) then

      ! fill ghost cells
      do n=1,nlevs
         call multifab_fill_boundary(n_old(n))
         call multifab_physbc(n_old(n),1,scal_bc_comp,nspecies, &
                              the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

   else

      if (model_file_init == 0) then
         ! initialize with a subroutine

         call init_n(mla,n_old,dx,the_bc_tower)

      else
         ! initialize from model files
         
         if(dm==2) then
            if (model_file_init>0) then ! Fortran (column-major) order
               allocate(input_array(n_cells(1),n_cells(2),1))
            else ! C (row-major) order   
               allocate(input_array(n_cells(2),n_cells(1),1))
            end if 
         else
            if (model_file_init>0) then ! Fortran (column-major) order
               allocate(input_array(n_cells(1),n_cells(2),n_cells(3)))
            else ! C (row-major) order   
               allocate(input_array(n_cells(3),n_cells(2),n_cells(1)))
            end if   
         end if            

         do n=1,nspecies

            if (parallel_IOProcessor()) then

               ! read in model file for species n into IOProc
               print*,'reading in model_file: ',model_file(n)
               open(unit=100, file=model_file(n), status='old', action='read')
               read(100,*) input_array
               close(unit=100)

            end if

            ! broadcast input_array to all processors
            ! there is no higher-rank routine at present so just do this for now
            do k=1,n_cells(3)
            do j=1,n_cells(2)
               call parallel_bcast_dv(input_array(:,j,k))
            end do
            end do

            ! copy data from input_array into multifab
            call init_n_model(mla,n_old,dx,the_bc_tower,input_array,n)

         end do

         deallocate(input_array)

      end if
      
      if (abs(initial_variance_mass) .gt. 0.d0) then
         ! For integer populations fluctuations have already been added at initialization via the Poisson distribution
         if(.not.integer_populations) call add_init_n_fluctuations(mla,n_old,dx,the_bc_tower)
      end if

   end if

   if (fixed_dt .gt. 0.d0) then
      dt = fixed_dt
      if (parallel_IOProcessor() ) then
         print*,''
         write(*,*) "Specified time step gives diff CFLs=", real(dt*D_Fick(1:nspecies)/dx(1,1)**2)
      end if
   else
      dt = cfl * dx(1,1)**2 / (maxval(D_Fick(1:nspecies)))
   end if

   ! compute n_steady for inhomogeneous_bc_fix
   if (inhomogeneous_bc_fix .and. temporal_integrator .ge. 0) then
      call compute_n_steady(mla,n_steady,dx,dt,the_bc_tower)
   end if

   if (temporal_integrator .lt. 0) then  ! unsplit schemes
      ! Donev: The code will work for a single cell also but may not be the most efficient, so issue warning:
      if ((multifab_volume(n_old(1))/nspecies)<=1) then
        if (parallel_IOProcessor() ) write(0,*), &
           "WARNING in advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for single cell"
      end if

      if(nreactions<1) then
        if (parallel_IOProcessor() ) write(0,*), &
           "WARNING in advance_reaction_diffusion: use splitting based schemes (temporal_integrator>=0) for diffusion only"
      end if
   end if

   !=====================================================================
   ! Initialize HydroGrid for analysis
   !=====================================================================
   if((hydro_grid_int>0) .or. (stats_int>0)) then
      narg = command_argument_count()
      farg = 1
      if (narg >= 1) then
         call get_command_argument(farg, value = fname)
         inquire(file = fname, exist = lexist )
         if ( lexist ) then
            un = unit_new()
            open(unit=un, file = fname, status = 'old', action = 'read')

            ! We will also pass temperature here but no additional scalars
            call initialize_hydro_grid(mla,n_old,dt,dx,namelist_file=un, &
                                       nspecies_in=nspecies, &
                                       nscal_in=0, &
                                       exclude_last_species_in=.false., &
                                       analyze_velocity=.false., &
                                       analyze_density=.true., &
                                       analyze_temperature=.false., &
                                       structFactMultiplier = volume_factor)

            close(unit=un)
         end if
      end if
   end if

   !=====================================================================
   ! Hydrogrid analysis and output for initial data
   !=====================================================================

   if (restart .lt. 0) then
      istep=0
      
      if (plot_int .gt. 0) then
         ! write initial plotfile
         if (parallel_IOProcessor()) then
            write(*,*), 'writing initial plotfile 0'
         end if
         call write_plotfile(mla,n_old,dx,time,istep)
      end if

      if (chk_int .ge. 0) then
         ! write initial checkpoint
         if (parallel_IOProcessor()) then
            write(*,*), 'writing initial checkpoint 0'
         end if
         call checkpoint_write(mla,n_old,time,dt,istep)
      end if

      if ( stats_int .gt. 0) then
         ! write initial vertical and horizontal averages (hstat and vstat files)   
         call print_stats(mla,dx,istep,time,rho=n_old)
      end if
   
      ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
      if (n_steps_skip .eq. 0) then

         ! Add the initial snapshot to the average in HydroGrid
         if (hydro_grid_int > 0) then
            call analyze_hydro_grid(mla,dt,dx,istep,rho=n_old)
         end if

         if ((hydro_grid_int > 0) .and. (n_steps_save_stats > 0)) then
            call save_hydro_grid(id=0, step=0)
         end if

      end if
   
   else ! Here we always write the starting point to check restarts are working
   
      istep = restart
      call write_plotfile(mla,n_old,dx,time,restart)   
      
   end if   

   !=======================================================
   ! Begin time stepping loop
   !=======================================================

   if (parallel_IOProcessor()) then
      print*,"BEGIN time loop istep =",istep,"dt =",dt,"time =",time
   end if
   runtime1 = parallel_wtime()

   do istep=init_step,max_step      

      ! Since we are still at the beginning of the time step here we use istep-1
      if ( (print_int>0) .and. (mod(istep-1,print_int)==0) ) then

         if (parallel_IOProcessor() ) then
            print*,''
            write(*,*) "At istep =",istep-1,"dt =",dt,"time =",time
         end if

          do spec=1,nspecies
             n_sum(spec) = multifab_sum_c(n_old(1),spec,1)
          end do
          if (parallel_IOProcessor() ) then
             !print*,time,' n_sum=',n_sum(:)
             write(*,*) istep-1, time, ' n_avg=', n_sum(:)/(multifab_volume(n_old(1))/nspecies)
          end if
          
          runtime2 = parallel_wtime()-runtime1
          call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
          if (parallel_IOProcessor() ) then
             print*,'Time to advance per timestep: ', runtime1/print_int,' seconds'
             print*,''
          end if
          runtime1=parallel_wtime()
                
       end if
       
       ! Donev: Temporary logging to analyze dynamics of mean
       if ((abs(n_steps_write_avg) > 0) .and. (mod(istep-1,abs(n_steps_write_avg))==0)) then
          do spec=1,nspecies
             n_sum(spec) = multifab_sum_c(n_old(1),spec,1)
          end do
          cellvolume=product(dx(1,1:MAX_SPACEDIM))*volume_factor ! Total system volume
          if (parallel_IOProcessor()) then
          !!if (n_steps_write_avg==2) then ! Write *only* the first species number density
          !!   ! Useful to save space in files where there is only one independent reaction
          !!   write(9,*) real(time), real(n_sum(1)/(multifab_volume(n_old(1))/nspecies))
          !!else if (n_steps_write_avg>0) then ! Write number densities
            if (n_steps_write_avg>0) then ! Write number densities
               write(9,*) real(time), real(n_sum(:)/(multifab_volume(n_old(1))/nspecies)) 
            else ! Here we write total number of molecules in the system instead of number densities             
               ! A. Donev: This is to match the particle code, write the same file it does
               write(21,*) real(time), real(sum(n_sum(:))*cellvolume), real(n_sum(:)*cellvolume)
            end if
          end if
       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! advance the solution by dt
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call advance_timestep(mla,n_old,n_new,n_steady,dx,dt,the_bc_tower)
       time = time + dt

       ! write a plotfile
       if (plot_int .gt. 0 .and. mod(istep,plot_int) .eq. 0) then
          call write_plotfile(mla,n_new,dx,time,istep)
       end if

       ! write a checkpoint
       if (chk_int .gt. 0 .and. mod(istep,chk_int) .eq. 0) then
          call checkpoint_write(mla,n_new,time,dt,istep)
       end if

       ! print out projection (average) and variance
       if ( (stats_int > 0) .and. &
            (mod(istep,stats_int) .eq. 0) ) then
          ! Compute vertical and horizontal averages (hstat and vstat files)   
          call print_stats(mla,dx,istep,time,rho=n_new)
       end if

       if (istep .ge. n_steps_skip) then

          ! Add this snapshot to the average in HydroGrid
          if ( (hydro_grid_int > 0) .and. &
               ( mod(istep,hydro_grid_int) .eq. 0 ) ) then
             call analyze_hydro_grid(mla,dt,dx,istep,rho=n_new)
          end if

          if ( (hydro_grid_int > 0) .and. &
               (n_steps_save_stats > 0) .and. &
               ( mod(istep,n_steps_save_stats) .eq. 0 ) ) then
             call save_hydro_grid(id=istep/n_steps_save_stats,step=istep)
          end if
       
       end if

       ! set old state to new state
       do n=1,nlevs
          call multifab_copy_c(n_old(n),1,n_new(n),1,nspecies,n_old(n)%ng)
       end do

   end do

   if (parallel_IOProcessor()) then
      print*,"END time loop istep =",istep,"dt =",dt,"time =",time
   end if

   !=======================================================
   ! Destroy multifabs and layouts
   !=======================================================

   if((hydro_grid_int>0) .or. (stats_int>0)) then
      call finalize_hydro_grid()
   end if

   call destroy_mass_stochastic(mla)

   do n=1,nlevs
      call multifab_destroy(n_new(n))
      call multifab_destroy(n_old(n))
      if (inhomogeneous_bc_fix .and. temporal_integrator .ge. 0) then
         call multifab_destroy(n_steady(n))
      end if
   end do

   if (use_bl_rng) then
      call rng_destroy()
   end if

   call destroy(mla)
   call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
