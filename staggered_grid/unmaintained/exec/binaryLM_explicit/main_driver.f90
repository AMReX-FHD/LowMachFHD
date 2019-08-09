subroutine main_driver()

  use advance_module
  use compute_dt_module
  use initialize_module
  use probin_module

  ! A. Donev:
  !------------------
  use ParallelRNGs
  
  implicit none


  integer    :: dm,n,istep
  real(dp_t) :: time,dt,dtold,dt_lev
  integer    :: init_step

  real(dp_t)  , pointer     :: dx(:,:)
  type(ml_layout)           :: mla

  type(multifab), pointer     ::     mold(:,:) ! face-based
  type(multifab), pointer     ::     sold(:)   ! cell-centered

  type(multifab), allocatable ::     mnew(:,:) ! face-based
  type(multifab), allocatable ::     snew(:)   ! cell-centered

  type(multifab), allocatable ::     umac(:,:) ! face-based
  type(multifab), allocatable :: rho_face(:,:) ! face-based

  type(bc_tower) ::  the_bc_tower

  integer :: namelist_file ! A. Donev: Where to read inputs from

  logical :: exit_timeloop_after_analysis

  exit_timeloop_after_analysis = .false.

  dt = 1.d99
  dtold = 1.d99

  call probin_init(namelist_file)

  ! Initialize random numbers *after* the global (root) seed has been set:
  if(use_stochastic_forcing) call SeedParallelRNG(seed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the grids and the data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (restart >= 0) then

     call initialize_from_restart(mla,restart,time,dt,dx,sold,mold,the_bc_tower)
     init_step = restart+1

  else if (fixed_grids /= '') then

     call initialize_with_fixed_grids(mla,dx,mold,sold,the_bc_tower)
     init_step = 1
     time = 0.d0

  else  ! Adaptive gridding

     call initialize_with_adaptive_grids(mla,dx,mold,sold,the_bc_tower)
     init_step = 1
     time = 0.d0

  end if

  nlevs = mla%nlevel
  dm = mla%dim

  ! now that visc_coef is set we can set visc_type
  ! this logic can/should probably be moved elsewhere...
  if (visc_coef .ge. 0.d0) then
     visc_type = 2
  else
     visc_type = -2
  end if

  ! try to set time step to something reasonable
  if (fixed_dt .gt. 0.d0) then

     dt = fixed_dt
     if (parallel_IOProcessor() .and. verbose .ge. 1) then
        print*,'Setting dt to fixed_dt =',fixed_dt
     end if

  else

     allocate(umac(nlevs,dm))
     allocate(rho_face(nlevs,dm))
     do n=1,nlevs
        do i=1,dm
           call multifab_build_edge(umac(n,i),mla%la(n),1,ng_mom,i)
           call multifab_build_edge(rho_face(n,i),mla%la(n),1,1,i)
        end do
     end do

     ! compute rho on faces
     call average_cc_to_face(mla%nlevel,sold,rho_face,1,dm+1,1,the_bc_tower%bc_tower_array)

     call convert_m_to_umac(mla,rho_face,mold,umac,.true.)

     do n = 1,nlevs
        call compute_dt(n,umac(n,:),dx(n,:),dtold,dt_lev)
        dt = min(dt,dt_lev)
     end do

     do n=1,nlevs
        do i=1,dm
           call multifab_destroy(umac(n,i))
           call multifab_destroy(rho_face(n,i))
        end do
     end do
     deallocate(umac,rho_face)

     if (dt .eq. 1.d99) then
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
           print*,'Initial Velocity too small to set a time step'
           print*,'Using dt = dx'
        end if
        dt = dx(nlevs,1)
     end if

     if (init_shrink .gt. 0.d0) then
        dt = dt * init_shrink
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
           print*,'Multiplying dt by init_shrink; dt=',dt
        end if
     end if


  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate new-time state and temp variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(mnew(nlevs,dm),snew(nlevs))

  ! build mnew and snew
  do n=1,nlevs
     do i=1,dm
        call multifab_build_edge(mnew(n,i), mla%la(n), 1, ng_mom, i)
     end do
     call multifab_build(snew(n), mla%la(n), nscal, ng_scal)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Complete the initialization of the analysis routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(abs(hydro_grid_int)>0 .or. stats_int>0) then
     call initialize_hydro_grid(mla,sold,mold,dt,dx,namelist_file)
  end if
  close(namelist_file)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the real integration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0)) then

     TimeLoop: do istep = init_step, max_step+n_steps_skip+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Advance the time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (stop_time >= 0.d0 .and. time >= stop_time) then
           exit_timeloop_after_analysis = .true.
        end if

        if (istep .eq. max_step+n_steps_skip+1) then
           exit_timeloop_after_analysis = .true.
        end if

        call advance_timestep(mla,sold,mold,snew,mnew,the_bc_tower,time, &
             dt,dx,istep,exit_timeloop_after_analysis)

        if (exit_timeloop_after_analysis) then
           exit TimeLoop
        end if

        time = time + dt

        ! copy new solution into old solution - valid region only
        do n = 1,nlevs
           do i=1,dm
              call multifab_copy_c(mold(n,i),1,mnew(n,i),1,1,0)
           end do
           call multifab_copy_c(sold(n),1,snew(n),1,nscal,0)
        end do

        if ( (print_int > 0) .and. (mod(istep,print_int) .eq. 0) ) then
           if ( parallel_IOProcessor() ) then
              write(*, "(' END OF STEP = ',i0,1x,' TIME = ',G17.9,1x,'DT = ',G17.9)") &
                   istep,time,dt
              write(*,*) "-------------------------------------"     
           end if

           call print_and_reset_fab_byte_spread()
        end if

     end do TimeLoop ! istep loop

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalization/cleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(abs(hydro_grid_int)>0 .or. stats_int>0) then
     call finalize_hydro_grid()
  end if

  do n = 1,nlevs
     do i=1,dm
        call multifab_destroy(mold(n,i))
        call multifab_destroy(mnew(n,i))
     end do
     call multifab_destroy(sold(n))
     call multifab_destroy(snew(n))
  end do

  deallocate(mold,sold)
  deallocate(mnew,snew)

  deallocate(dx)

  call bc_tower_destroy(the_bc_tower)

  call destroy(mla)

  call probin_close()

  if ( verbose > 0 ) then
     if ( parallel_IOProcessor() ) then
        print *, 'MEMORY STATS AT END OF RUN '
        print*, ' '
     end if
     call print(multifab_mem_stats(),    "    multifab")
     call print(fab_mem_stats(),         "         fab")
     call print(boxarray_mem_stats(),    "    boxarray")
     call print(layout_mem_stats(),      "      layout")
     call print(boxassoc_mem_stats(),    "    boxassoc")
     call print(fgassoc_mem_stats(),     "     fgassoc")
     call print(syncassoc_mem_stats(),   "   syncassoc")
     call print(copyassoc_mem_stats(),   "   copyassoc")
     call print(fluxassoc_mem_stats(),   "   fluxassoc")
     if ( parallel_IOProcessor() ) print*, ''
  end if

end subroutine main_driver
