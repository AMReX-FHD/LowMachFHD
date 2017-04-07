module advance_module 

  use probin_module
  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use macproject_module
  use convert_module
  use mk_advective_fluxdiv_module
  use mk_diffusive_fluxdiv_module
  use mk_external_force_module
  use mk_stochastic_fluxdiv_module
  use analysis_module
  use analyze_spectra_module
  use exact_solutions_module
  use project_onto_eos_module

  integer, save :: last_plt_written = -1
  integer, save :: last_chk_written = -1

contains

subroutine advance_timestep(mla,sold,mold,snew,mnew,the_bc_tower, &
                            time,dt,dx,istep,exit_timeloop_after_analysis)

   implicit none

   type(ml_layout), intent(in   ) :: mla
   type(multifab) , intent(inout) ::    sold(:)
   type(multifab) , intent(inout) ::    mold(:,:)
   type(multifab) , intent(inout) ::    snew(:)
   type(multifab) , intent(inout) ::    mnew(:,:)
   type(bc_tower) , intent(in   ) :: the_bc_tower
   real(dp_t)     , intent(in   ) :: time
   real(dp_t)     , intent(inout) :: dt
   real(dp_t)     , intent(in   ) :: dx(:,:)
   integer        , intent(in   ) :: istep
   logical        , intent(in   ) :: exit_timeloop_after_analysis

   ! rho and rho1 averaged to faces
   type(multifab) :: s_face(mla%nlevel,mla%dim)

   ! face-based (MAC) velocities
   type(multifab) :: umac(mla%nlevel,mla%dim)

   ! cell-centered primitive variables
   type(multifab) :: prim(mla%nlevel)

   ! right-hand-side for projections
   type(multifab) :: divu_rhs(mla%nlevel)

   ! stochastic terms
   type(multifab) :: stoch_m_force(mla%nlevel,mla%dim)
   type(multifab) :: stoch_s_force(mla%nlevel)

   ! updates to m and s
   type(multifab) :: m_update(mla%nlevel,mla%dim)
   type(multifab) :: s_update(mla%nlevel)

   ! diffusion coefficients for velocity and concentration
   type(multifab) ::       eta(mla%nlevel)
   type(multifab) :: eta_nodal(mla%nlevel)         ! averaged to nodes (2D only)
   type(multifab) ::  eta_edge(mla%nlevel,3)       ! averaged to edges (3D only; components are xy/xz/yz edges)
   type(multifab) ::       chi(mla%nlevel)
   type(multifab) ::  chi_face(mla%nlevel,mla%dim) ! averaged to faces
   type(multifab) ::     kappa(mla%nlevel)

   logical :: nodal_temp(mla%dim)

   integer :: i,n,dm,nlevs,n_rngs

   real(kind=dp_t) :: dtold,dt_lev,stage_time
   real(kind=dp_t) :: S_fac,advance_by,stochastic_w1,stochastic_w2

   real(kind=dp_t), allocatable :: weights(:) ! Weights for RNGs

   logical :: skipping_macproj_on_restart, initialize_fluctuations

   stage_time = time

   stochastic_w1 = 0.d0
   stochastic_w2 = 0.d0

   S_fac = 1.d0/rhobar(1) - 1.d0/rhobar(2)

   dm    = mla%dim
   nlevs = mla%nlevel

   do n=1,nlevs

      ! staggered fields
      do i=1,dm
        call multifab_build_edge(         umac(n,i), mla%la(n), 1    , ng_mom, i)
        call multifab_build_edge(     m_update(n,i), mla%la(n), 1    , 0,      i)
        call multifab_build_edge(stoch_m_force(n,i), mla%la(n), 1    , 0,      i)
        call multifab_build_edge(       s_face(n,i), mla%la(n), nscal, 1,      i)
        call multifab_build_edge(     chi_face(n,i), mla%la(n), 1    , 0,      i)

        ! Make sure these have valid values:
        call setval(chi_face(n,i), abs(diff_coef), all=.true.)

      end do

      ! cell-centered fields
      call multifab_build(         prim(n), mla%la(n), nscal, ng_scal)
      call multifab_build(     s_update(n), mla%la(n), nscal, 0)
      call multifab_build(stoch_s_force(n), mla%la(n), nscal, 0)
      call multifab_build(     divu_rhs(n), mla%la(n),     1, 0)
      call multifab_build(          eta(n), mla%la(n),     1, 1)
      call multifab_build(          chi(n), mla%la(n),     1, 1)
      call multifab_build(        kappa(n), mla%la(n),     1, 1)

      ! Make sure these have valid values:
      call setval(eta(n)  , abs(visc_coef), all=.true.)
      call setval(chi(n)  , abs(diff_coef), all=.true.)
      call setval(kappa(n), abs(bulk_visc), all=.true.)

      ! These are multiplied by zero in the Euler step
      ! but they should be initialized to avoid signaling NaNs
      do i=1,dm
         call setval(mnew(n,i), 0.d0, all=.true.)
      end do   
      call setval(snew(n), 0.d0, all=.true.)

   end do

   ! nodal (in 2D) and edge-based (in 3D) eta
    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_build_nodal(eta_nodal(n),mla%la(n),1,0)
          call setval(eta_nodal(n),abs(visc_coef),all=.true.)
       end do
    else
       do n=1,nlevs
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(eta_edge(n,1),mla%la(n),1,0,nodal_temp)
          call setval(eta_edge(n,1),abs(visc_coef),all=.true.)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(eta_edge(n,2),mla%la(n),1,0,nodal_temp)
          call setval(eta_edge(n,2),abs(visc_coef),all=.true.)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(eta_edge(n,3),mla%la(n),1,0,nodal_temp)
          call setval(eta_edge(n,3),abs(visc_coef),all=.true.)
       end do
    end if

   ! if stoch_forcing_stop is positive, stop the stochastic forcing after that step
   if (use_stochastic_forcing .and. stoch_forcing_stop .ge. 0 .and. &
        istep .gt. stoch_forcing_stop) then
      use_stochastic_forcing = .false.
   end if

   if (use_stochastic_forcing) then

      if (temporal_scheme <= 0) then
        ! We always use the same random increments (this gives order 1 in the weak sense)
        n_rngs=0
      else if (temporal_scheme == 1) then ! Trapezoidal
        ! We need a single set of random numbers per step
        n_rngs=1
      else ! Midpoint or RK3
        n_rngs=2  
      end if
      allocate(weights(n_rngs))

      ! compute and store RANDOM NUMBERS for stochastic forcing terms
      call create_random_increments(mla,n_rngs)

   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN EULER STEP (PREDICTOR)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The first step may have different lengths:        
   if (abs(temporal_scheme) .eq. 2) then
      advance_by=0.5d0 ! Only move to the middle of the time-step
   else
      advance_by=1.0d0 ! Take a full time step
   end if

   if (use_stochastic_forcing) then
      if (temporal_scheme == 4) then ! RK3 with new experimental weights
         stochastic_w1 = 1.d0
         stochastic_w2 = (2*sqrt(2.0d0)+sqrt(3.0d0))/5 ! Either minus or plus sign can be chosen, plus is better
      else if (temporal_scheme == 3) then ! RK3 scheme
         stochastic_w1 = 1.d0
         stochastic_w2 = -sqrt(3.d0)
      else if (temporal_scheme > 0) then
         ! The stochastic flux is weighted not by dt, but by sqrt(dt):
         stochastic_w1 = sqrt(1.0d0/advance_by)
         stochastic_w2 = 0.d0
      else   
         ! Pretend that the stochastic flux is constant and fixed during the time step:
         stochastic_w1 = advance_by
         stochastic_w2 = 0.d0
      end if  
   end if

   ! don't call macproject upon restart - it was done before writing the checkpoint
   skipping_macproj_on_restart = .false.
   if ((restart .ge. 0) .and. (istep .eq. restart+1)) then
      skipping_macproj_on_restart = .true.
   end if
   ! Should we complete the initialization by adding equilibrium fluctuations?
   initialize_fluctuations = .false.
   if ((abs(initial_variance)>0) .and. (restart <= 0) .and. (istep == 1)) then
      ! Note that one can add fluctuations even if use_stochastic_forcing=F to get random initial conditions
      initialize_fluctuations = .true.
   end if  

   ! Note that this both finishes the previous time step (in sold/mold)
   ! and computes fluxes at sold/mold:
   call compute_update(sold,mold,.true.)

   ! now call mac projection in every subsequent call to compute_update
   skipping_macproj_on_restart = .false.
   initialize_fluctuations = .false.

   call analyze_snapshot() ! Examine sold, mold, umac and prim

   if (exit_timeloop_after_analysis) then       
      call delete_temp_multifabs()
      return
   end if

   ! advance using forward Euler advection, diffusion, and forcing
   ! new = old + advance_by*dt*update
   call weighted_advance(1.0d0, 0.0d0, advance_by)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN CORRECTOR
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!    

   if (abs(temporal_scheme) .eq. 1) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! TRAPEZOIDAL CORRECTOR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      stage_time = time + dt

      if (use_stochastic_forcing) then
         stochastic_w1 = 1.d0
         stochastic_w2 = 0.d0
      end if

      call compute_update(snew,mnew,.false.)

      ! advance using explicit trapezoidal advection, diffusion, and forcing
      ! new = (old + new + dt*update)/2
      call weighted_advance(0.5d0, 0.5d0, 0.5d0)

   else if (abs(temporal_scheme) .eq. 2) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MIDPOINT CORRECTOR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      stage_time = time + 0.5d0*dt

      if (use_stochastic_forcing) then
         stochastic_w1 = sqrt(0.5d0)
         stochastic_w2 = sqrt(0.5d0)
      end if

      call compute_update(snew,mnew,.false.)

      ! advance using explicit midpoint advection, diffusion, and forcing
      ! new = old + dt*update
      call weighted_advance(1.0d0, 0.0d0, 1.0d0)

   else if (abs(temporal_scheme) >= 3) then ! Three-stage TVD RK3 scheme

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! RK3 STAGE 2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      stage_time = time + 1.d0/3.d0*dt

      if (use_stochastic_forcing) then
         if(abs(temporal_scheme)==4) then
            stochastic_w1 = 1.d0
            stochastic_w2 = (-4*sqrt(2.0d0)+3*sqrt(3.0d0))/5 ! Plus sign for sqrt(3) term is better
         else   
            stochastic_w1 = 1.d0
            stochastic_w2 = sqrt(3.d0)
         end if
      end if

      call compute_update(snew,mnew,.false.)

      ! advance to midpoint (n+1/2) using explicit 3/4 + 1/4 weighting:
      ! new = 3/4*old + 1/4*(new + dt*update)
      call weighted_advance(0.75d0, 0.25d0, 0.25d0)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! RK3 STAGE 3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      stage_time = time + 2.d0/3.d0*dt

      if (use_stochastic_forcing) then
         if(abs(temporal_scheme)==4) then
            stochastic_w1 = 1.d0
            stochastic_w2 = (sqrt(2.0d0)-2*sqrt(3.0d0))/10 ! Minus sign for sqrt(3.0) term is better
         else   
            stochastic_w1 = 1.d0
            stochastic_w2 = 0.d0
         end if
      end if

      call compute_update(snew,mnew,.false.)

      ! advance to end point (n+1) using explicit 1/3 + 2/3 weighting:
      ! new = 1/3*old + 2/3*(new + dt*update)
      call weighted_advance(1.0d0/3.d0, 2.0d0/3.d0, 2.0d0/3.d0)

   end if

   if (enforce_eos) then
      call project_onto_eos(mla,snew)
   end if

   call delete_temp_multifabs()
    
contains

   subroutine delete_temp_multifabs()

      if (use_stochastic_forcing) then ! Get rid of the random number multifabs
         call destroy_random_increments(mla,n_rngs)
      end if

      do n = 1,nlevs
         do i = 1,dm
            call multifab_destroy(umac(n,i))
            call multifab_destroy(m_update(n,i))
            call multifab_destroy(stoch_m_force(n,i))
            call multifab_destroy(s_face(n,i))
            call multifab_destroy(chi_face(n,i))
         end do
         call multifab_destroy(prim(n))
         call multifab_destroy(s_update(n))
         call multifab_destroy(stoch_s_force(n)) 
         call multifab_destroy(divu_rhs(n))
         call multifab_destroy(eta(n))
         call multifab_destroy(chi(n))
         call multifab_destroy(kappa(n))
      end do


    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_destroy(eta_nodal(n))
       end do
    else
       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(eta_edge(n,i))
          end do
       end do
    end if

   end subroutine

   ! Computes primitive variables, projects velocities, and computes fluxes
   ! at reference (current) state s_in and m_in
   subroutine compute_update(s_in,m_in,is_first_stage)

      type(multifab) , intent(inout) :: s_in(:)
      type(multifab) , intent(inout) :: m_in(:,:)
      logical        , intent(in)    :: is_first_stage
      
      real(dp_t) :: c_min, c_max

      ! set s_update and m_update to zero
      call reset_m_s_updates()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Scalars: Non-advective fluxes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! convert s_in to prim in valid region
      ! ----------------------
      call convert_cons_to_prim(mla,s_in,prim,.true.)

      if(initialize_fluctuations) then      
         if ( parallel_IOProcessor() ) then
            write(*,*) "Initializing concentration fluctuations at step ", istep, &
                 " with variance ", initial_variance*variance_coeff*conc_scal
         end if
         ! The densities ought to follow the EOS after adding fluctuations!
         call add_concentration_fluctuations(mla, dx, initial_variance*variance_coeff*conc_scal, &
                                             prim, s_in) ! In valid region only
      end if

      ! fill periodic/interior and physical ghost cells for prim
      do n=1,nlevs
         call multifab_fill_boundary(prim(n))
         call multifab_physbc(prim(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n))
      end do
      
      if(print_cminmax) then
         do n=1,nlevs
            c_min=min_val(prim(n), 2, 1, all=.true.)
            c_max=max_val(prim(n), 2, 1, all=.true.)
            if ( parallel_IOProcessor() .and. (c_min<-0.25_dp_t) .or. (c_max>1.25_dp_t) ) then
               write(*,*) istep," MIN/MAX CONC", c_min, c_max
            end if
         end do            
      end if

      ! convert prim to s_in in valid and ghost region
      ! now both prim and s_in are valid and consistent everywhere
      call convert_cons_to_prim(mla,s_in,prim,.false.)

      ! compute rho and rho1 on faces and put it into s_face
      call average_cc_to_face(mla%nlevel,s_in,s_face,1,dm+1,2,the_bc_tower%bc_tower_array)

      ! ----------------------
      if(initialize_fluctuations) then      
         if ( parallel_IOProcessor() ) then
            write(*,*) "Initializing momentum fluctuations at step ", istep, " with variance ", initial_variance*variance_coeff
         end if
         ! This adds fluctuations everywhere but ghost cells and boundaries will be fixed shortly
         call add_momentum_fluctuations(mla, dx, initial_variance*variance_coeff, s_in, s_face, m_in, umac)
      end if

      ! ----------------------
      ! convert m_in to umac in valid region
      call convert_m_to_umac(mla,s_face,m_in,umac,.true.)

      do n=1,nlevs
         do i=1,dm
            ! fill periodic and interior ghost cells
            call multifab_fill_boundary(umac(n,i))
            ! set normal velocity on physical domain boundaries to zero
            call multifab_physbc_domainvel(umac(n,i),1,i,1,the_bc_tower%bc_tower_array(n))
            ! set the remaining physical domain boundary ghost cells
            call multifab_physbc_macvel(umac(n,i),1,i,1, &
                                        the_bc_tower%bc_tower_array(n),dx(n,:))
         end do
      end do
      ! ----------------------            

      if (diff_coef < 0.d0) then
         ! compute spatially varying chi
         call compute_chi(mla,chi,prim,dx)
         ! average chi to faces
         call average_cc_to_face(mla%nlevel,chi,chi_face,1,dm+1,1,the_bc_tower%bc_tower_array)
      end if

      ! compute spatially varying eta
      if (visc_coef < 0.d0) then
         call compute_eta(mla,eta,prim,dx)
         if (dm .eq. 2) then
            call average_cc_to_node(mla%nlevel,eta,eta_nodal,1,dm+1,1,the_bc_tower%bc_tower_array)
         else
            call average_cc_to_edge(mla%nlevel,eta,eta_edge,1,dm+1,1,the_bc_tower%bc_tower_array)
         end if
      end if

      ! compute spatially varying kappa
      if (bulk_visc < 0.d0) then
         call compute_kappa(mla,kappa,prim,dx)
      end if

      ! compute stochastic forcing for s/m and add to s_update/m_update
      ! this could be a function of s and/or prim, but not m
      if (use_stochastic_forcing) then
         if (is_first_stage .or. temporal_scheme .gt. 0) then
            call add_stochastic_fluxes()
         else  
            ! Reuse the previously-generated stochastic increment:
            call reuse_stochastic_fluxes()
         end if
      end if

      ! compute del dot (rhoD grad c) and add it to s_update
      if (diff_coef .ne. 0.d0) then
         call mk_diffusive_rhoc_fluxdiv(mla,s_update,prim,s_face,chi_face,umac,dx, &
                                        the_bc_tower%bc_tower_array)
      end if

      ! compute external forcing for s and add to s_update
      ! This can depend on the scalars but not on velocity!
      call mk_external_s_force(mla,s_update,s_in,dx,stage_time)

      ! compute rhs for projection, S = S_fac * s_update(rho*c)
      do n=1,nlevs
         call multifab_copy_c(divu_rhs(n), 1, s_update(n), 2, 1, 0)
         call multifab_mult_mult_s_c(divu_rhs(n), 1, S_fac, 1, 0)    
      end do

      ! project velocity subject to div(u)=S
      ! Note that the projection depends on density so we pass s_in here
      ! Also note that the projection does not use or touch ghost cells for umac
      ! Boundary conditions are handled entirely via pressure in the_bc_tower
      if (.not. skipping_macproj_on_restart) then
         call macproject(mla,umac,s_in,dx,the_bc_tower,divu_rhs)
      end if

      ! fill periodic/interior and physical ghost cells for umac
      do n=1,nlevs
         do i=1,dm
            call multifab_fill_boundary(umac(n,i))
            call multifab_physbc_macvel(umac(n,i),1,i,1, &
                                        the_bc_tower%bc_tower_array(n),dx(n,:))
         end do
      end do

      ! convert umac to m_in in valid and ghost region
      ! now both umac and m_in are valid and consistent everywhere
      call convert_m_to_umac(mla,s_face,m_in,umac,.false.)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Velocities: Non-advective fluxes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Note that we do NOT compute fluxes for the normal velocities at physical boundaries
      ! since these are not real variables but rather fixed by the BCs

      ! compute del dot (eta*grad(u) + eta*grad(u)^T) and add it to m_update
      if (visc_coef .ne. 0.d0) then
         call mk_diffusive_m_fluxdiv(mla,m_update,umac,eta,eta_nodal,eta_edge, &
                                     kappa,dx,the_bc_tower%bc_tower_array)
      endif

      ! compute external forcing for m and add to m_update
      ! This can depend on the scalars but for now not on velocity
      call mk_external_m_force(mla,m_update,s_in,dx,stage_time)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Advective fluxes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! compute advection terms for s/m and add to s_update/m_update
      if(.not. fake_projection) then
         ! add scalar advection to s_update
         call mk_advective_s_fluxdiv(mla,umac,s_face,s_update,dx)

         ! add momentum advection to m_update
         call mk_advective_m_fluxdiv(mla,umac,m_in,m_update,dx)
      end if

      do n=1,nlevs
         do i=1,dm
            ! set m_update on physical domain boundaries to zero
            call multifab_physbc_domainvel(m_update(n,i),1,i,1, &
                                           the_bc_tower%bc_tower_array(n))
         end do
      end do
      
   end subroutine compute_update

   subroutine weighted_advance(w1,w2,w3)

      ! snew = w1*sold + w2*snew + dt*w3*s_update
      ! mnew = w1*mold + w2*mnew + dt*w3*m_update

      real(kind=dp_t), intent(in   ) :: w1,w2,w3

      integer :: i,n

      do n=1,nlevs
         call multifab_mult_mult_s(snew(n),w2)
         call saxpy(snew(n),w1,sold(n))
         call saxpy(snew(n),dt*w3,s_update(n))
         do i=1,dm
            call multifab_mult_mult_s(mnew(n,i),w2)
            call saxpy(mnew(n,i),w1,mold(n,i))
            call saxpy(mnew(n,i),dt*w3,m_update(n,i))
         end do
      end do

   end subroutine weighted_advance

   subroutine reset_m_s_updates()

      do n=1,nlevs
         call setval(s_update(n), 0.d0, all=.true.)
         do i=1,dm
            call setval(m_update(n,i), 0.d0, all=.true.)
         end do
      end do    

   end subroutine

   ! The stochastic fluxes are computed based on the present state
   ! But the random numbers might be stored and reused from previous stages
   subroutine add_stochastic_fluxes()

      ! mk_stochastic_fluxdiv *increments* stoch_m_force and stoch_s_force, so we initialize:
      do n=1,nlevs
         call setval(stoch_s_force(n), 0.d0, all=.true.)
         do i=1,dm
            call setval(stoch_m_force(n,i), 0.d0, all=.true.)
         end do
      end do    

      if(n_rngs>=1) weights(1)=stochastic_w1
      if(n_rngs>=2) weights(2)=stochastic_w2
      call mk_stochastic_fluxdiv(mla,the_bc_tower%bc_tower_array,stoch_m_force,stoch_s_force, &
                                 s_face,eta,eta_nodal,eta_edge,chi_face,umac,dx,dt,weights)

      if(n_rngs==0) then ! There is no weighting in mk_stochastic_fluxdiv, so we do it here
         do n=1,nlevs
            call saxpy(s_update(n), stochastic_w1, stoch_s_force(n))
            do i=1,dm
               call saxpy(m_update(n,i), stochastic_w1, stoch_m_force(n,i))
            end do
         end do
      else ! The weighting was already done in mk_stochastic_fluxdiv, so just add this to the sum
         do n=1,nlevs
            call multifab_plus_plus(s_update(n), stoch_s_force(n))
            do i=1,dm
               call multifab_plus_plus(m_update(n,i), stoch_m_force(n,i))
            end do
         end do
      end if   

   end subroutine

   ! This does not regenerate random fluxes but simply reuses the existing ones:
   subroutine reuse_stochastic_fluxes()

      do n=1,nlevs
         call multifab_plus_plus(s_update(n), stoch_s_force(n))
         do i=1,dm
            call multifab_plus_plus(m_update(n,i), stoch_m_force(n,i))
         end do
      end do

   end subroutine

   ! if this is the first Euler stage, record a snapshot:
   ! -compute dt
   ! -write plot/checkfile
   ! -print error norms
   ! -print eos drift
   ! -call analysis routines
   subroutine analyze_snapshot()

      use compute_dt_module
      use write_plotfile_module
      use checkpoint_module

      ! compute dt
      if (fixed_dt .gt. 0.d0) then

         dt = fixed_dt
         if (parallel_IOProcessor() .and. verbose .ge. 1) then
            print*,'Setting dt to fixed_dt =',fixed_dt
         end if

      else

         dtold = dt
         dt = 1.d20

         do n = 1,nlevs
            call compute_dt(n,umac(n,:),dx(n,:),dtold,dt_lev)
            dt = min(dt,dt_lev)
         end do

         if (parallel_IOProcessor() .and. verbose .ge. 1) then
            print*,"Estimating new dt for step",istep
            print*,"Old dt=", dtold
            print*,"New dt=", dt
         end if

      end if

      if (stop_time >= 0.d0 .and. time+dt > stop_time) then
         dt = stop_time - time
         if (parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, "Stop time limits dt =",dt
         end if
      end if

      ! For uniformity with other interval counters perhaphs we should skip n_steps_skip here?
      if ( (print_int > 0) .and. (mod(istep,print_int) .eq. 0) ) then

         if ( .false..and.parallel_IOProcessor() ) then
            write(*, "('BEGINNING OF STEP = ',i6,1x,' TIME = ',f16.10,1x,'DT = ',f14.9)") &
                 istep,time,dt
            write(*,*) "-------------------------------------"     
         end if

         ! compute error norms
         if (print_error_norms) then
            call print_errors(mla,sold,snew,mold,mnew,dx,time)
         end if

         ! compute drift from eos
         if (print_eos_error) then
            call eos_check(mla,sold)
         end if

         if (print_conserved) then
            call sum_mass_momentum(mla,sold,mold)
         end if

      end if
      
      if ( (istep > n_steps_skip) ) then

         ! write checkfile
         if ( (chk_int > 0) .and. &
              ( (mod(istep-n_steps_skip-1,chk_int) .eq. 0) .or. &
                 exit_timeloop_after_analysis) ) then
            call write_checkfile(mla,sold,mold,time,dt,istep-n_steps_skip-1)
            last_chk_written = istep-n_steps_skip-1
         end if

         ! Extract just the solenoidal piece of umac if desired
         ! No need to mess with ghost/boundary values here as they are not used
         if ( (proj_int > 0) .and. &
              (mod(istep-n_steps_skip-1,proj_int) .eq. 0) ) then
              call macproject(mla,umac,sold,dx,the_bc_tower)
         end if

         ! write plotfile
         if ( (plot_int > 0) .and. &
              ( (mod(istep-n_steps_skip-1,plot_int) .eq. 0) .or. &
                 exit_timeloop_after_analysis) ) then
            call write_plotfile(mla,mold,umac,sold,dx,time,istep-n_steps_skip-1)
            last_plt_written = istep-n_steps_skip-1
         end if

         ! print out projection (average) and variance
         if ( (stats_int > 0) .and. &
              ( (mod(istep-n_steps_skip-1,stats_int) .eq. 0) .or. &
                 exit_timeloop_after_analysis) ) then
            call print_stats(mla,sold,mold,umac,prim,dx,istep-n_steps_skip-1,time)
            if (hydro_grid_int<0) then
               call analyze_hydro_grid(mla,sold,mold,umac,prim,dt,dx,istep-n_steps_skip-1,custom_analysis=.true.)
            end if   
         end if

         ! Add this snapshot to the average in HydroGrid
         if ( (hydro_grid_int > 0) .and. &
              ( mod(istep-n_steps_skip,hydro_grid_int) .eq. 0 ) ) then
            call analyze_hydro_grid(mla,sold,mold,umac,prim,dt,dx,istep-n_steps_skip-1,custom_analysis=.false.)
         end if

         if ( (hydro_grid_int > 0) .and. &
              (n_steps_save_stats > 0) .and. &
              ( mod(istep-n_steps_skip,n_steps_save_stats) .eq. 0 ) ) then
            call save_hydro_grid(id=(istep-n_steps_skip)/n_steps_save_stats, step=istep)
         end if

      end if

   end subroutine analyze_snapshot

end subroutine advance_timestep

end module advance_module 
