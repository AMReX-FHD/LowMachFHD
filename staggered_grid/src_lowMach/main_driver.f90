subroutine main_driver()

  use bl_IO_module
  use bl_error_module
  use init_lowmach_module
  use compute_mixture_properties_module
  use initial_projection_module
  use write_plotfile_module
  use advance_timestep_bousq_module
  use advance_timestep_inertial_module
  use advance_timestep_overdamped_module
  use advance_timestep_iterative_module
  use advance_timestep_imp_bousq_module
  use advance_timestep_inertial_midpoint_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_physbc_stag_module
  use analysis_module
  use analyze_spectra_module
  use div_and_grad_module
  use eos_check_module
  use estdt_module
  use stag_mg_layout_module
  use macproject_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use fill_rho_ghost_cells_module
  use ParallelRNGs 
  use bl_rng_module
  use bl_random_module
  use mass_flux_utilities_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use sum_momenta_module
  use restart_module
  use checkpoint_module
  use reservoir_bc_fill_module
  use utility_module
  use user_analysis
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, hydro_grid_int, &
                                  max_grid_size, n_steps_save_stats, n_steps_skip, &
                                  plot_int, chk_int, seed, stats_int, bc_lo, bc_hi, restart, &
                                  probin_common_init, print_int, nspecies, dx_saved, &
                                  advection_type, fixed_dt, dt_saved, max_step, cfl, &
                                  algorithm_type, variance_coef_mom, initial_variance_mom, &
                                  variance_coef_mass, barodiffusion_type, use_bl_rng, &
                                  density_weights, rhobar, rho0, analyze_conserved, rho_eos_form, &
                                  molmass, k_B, reset_tavg_vals, reset_tavg_step, stat_save_type, analyze_cuts

  use probin_multispecies_module, only: Dbar, start_time, probin_multispecies_init, c_init, T_init, c_bc

  use probin_gmres_module, only: probin_gmres_init

  use probin_charged_module, only: probin_charged_init, use_charged_fluid, dielectric_const, &
                                   dielectric_type, electroneutral, epot_mg_rel_tol, charge_per_mass, & 
                                   print_debye_len, bc_function_type, L_pos, L_trans, L_zero

  use probin_chemistry_module, only: probin_chemistry_init, nreactions

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2,Dbar_max,dt_diffusive
  integer                      :: n,nlevs,i,dm,istep,ng_s,init_step,n_Dbar,m
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: rho_old(:)
  type(multifab), allocatable  :: rhotot_old(:)
  type(multifab), allocatable  :: rho_new(:)
  type(multifab), allocatable  :: rhotot_new(:)
  type(multifab), allocatable  :: Temp(:)
  type(multifab), allocatable  :: Temp_ed(:,:)
  type(multifab), allocatable  :: diff_mass_fluxdiv(:)
  type(multifab), allocatable  :: stoch_mass_fluxdiv(:)
  type(multifab), allocatable  :: stoch_mass_flux(:,:)
  type(multifab), allocatable  :: umac(:,:)
  type(multifab), allocatable  :: mtemp(:,:)
  type(multifab), allocatable  :: rhotot_fc(:,:)
  type(multifab), allocatable  :: gradp_baro(:,:)
  type(multifab), allocatable  :: pi(:)
  type(multifab), allocatable  :: eta(:)
  type(multifab), allocatable  :: eta_ed(:,:)
  type(multifab), allocatable  :: kappa(:)

  real(kind=dp_t)              :: total_charge
  ! this is only used for implicit potential algorithms where we can't add this to
  ! the diffusive mass fluxes.
  type(multifab), allocatable  :: Epot_mass_fluxdiv(:)
  type(multifab), allocatable  :: charge_old(:)
  type(multifab), allocatable  :: charge_new(:)
  type(multifab), allocatable  :: permittivity(:)
  type(multifab), allocatable  :: grad_Epot_old(:,:)
  type(multifab), allocatable  :: grad_Epot_new(:,:)
  type(multifab), allocatable  :: Epot(:)
  type(multifab), allocatable  :: gradPhiApprox(:,:)

  type(multifab), allocatable  :: chem_rate(:)

  ! for writing time-averaged umac to plotfile
  type(multifab), allocatable :: umac_sum(:,:)
  type(multifab), allocatable :: umac_avg(:,:)

  ! for writing time-averaged species' densities to plotfile
  type(multifab), allocatable :: rho_sum(:) 
  type(multifab), allocatable :: rho_avg(:) 

  ! for writing time-averaged Electric_field to plotfile
  type(multifab), allocatable :: Epot_sum(:) 
  type(multifab), allocatable :: Epot_avg(:) 

  ! for algorithm_type=6, return total mass fluxes (diff + stoch + adv)
  type(multifab), allocatable :: total_mass_flux(:,:)
  
  ! For HydroGrid
  integer :: narg, farg, un, n_rngs_mass, n_rngs_mom, nscal
  character(len=128) :: fname
  logical :: lexist
  logical :: nodal_temp(3)

  ! misc
  real(kind=dp_t) :: max_charge, max_charge_abs, debye_len, sum_of_boundary_lens, delta_x

  ! DONEV FIXME
  real(kind=dp_t) :: rho_temp, w_temp(1:5), w_mol(1:3)
  real(kind=dp_t), parameter :: m_Na=3.817540700000000d-023, &	
                                m_Cl=5.887108600000000d-023, &
                                m_H =1.673723600000000d-024, &
                                m_OH=2.824068560000000d-023
  
  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  call probin_charged_init() 
  call probin_chemistry_init()

  ! first check that if we prescribe a nonconstant, piecewise cubic BC for Epot, 
  ! then the length of the various regions must equal Lx ie prob_hi(1)
  ! Additionally, the point at which L_pos/trans_zero regions change from one to another
  ! must align with cell boundaries. ie the x-loc of the | markers below must = M*dx, M = integer
  ! |---L_pos---|--L_trans--|---L_zero---|--L_trans--|---L_pos---| 
  if (bc_function_type.eq.1) then  
     sum_of_boundary_lens = 2.d0*L_pos + 2.d0*L_trans + L_zero
     if (sum_of_boundary_lens.ne.prob_hi(1)) then 
        call bl_error("If you impose a non-constant cubic BC for Epot, the regions summed length must match total domain length.")
     end if 
     delta_x = prob_hi(1)/n_cells(1)
     if (L_pos/delta_x.ne.floor(L_pos/delta_x)) then 
        call bl_error("L_pos must be an integer multiple of dx")
     end if 
     if (L_trans/delta_x.ne.floor(L_trans/delta_x)) then 
        call bl_error("L_trans must be an integer multiple of dx")
     end if 
     if (L_zero/delta_x.ne.floor(L_zero/delta_x)) then 
        call bl_error("L_zero must be an integer multiple of dx")
     end if 
  end if 

  ! for reservoirs, make sure the Dirichlet conditions for concentration sum to 1
  do i=1,dim_in
     if (bc_lo(i) .eq. NO_SLIP_RESERVOIR .or. bc_lo(i) .eq. SLIP_RESERVOIR) then
        c_bc(i,1,nspecies) = 1.d0 - sum(c_bc(i,1,1:nspecies-1))
     end if
     if (bc_hi(i) .eq. NO_SLIP_RESERVOIR .or. bc_hi(i) .eq. SLIP_RESERVOIR) then
        c_bc(i,2,nspecies) = 1.d0 - sum(c_bc(i,2,1:nspecies-1))
     end if
  end do

  if (use_bl_rng) then
     ! Build the random number engine and give initial distributions for the
     ! F_BaseLib/bl_random RNG module
     call rng_init()
  else
     ! Initialize random numbers *after* the global (root) seed has been set:
     ! This is for the RNG module that sits in Hydrogrid
     call SeedParallelRNG(seed)
  end if

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm),pmask(dm))
  allocate(rho_old(nlevs),rhotot_old(nlevs),pi(nlevs))
  allocate(rho_new(nlevs),rhotot_new(nlevs))
  allocate(rho_sum(nlevs)) 
  allocate(rho_avg(nlevs)) 
  allocate(Temp(nlevs))
  allocate(diff_mass_fluxdiv(nlevs),stoch_mass_fluxdiv(nlevs))
  allocate(stoch_mass_flux(nlevs,dm))
  allocate(umac(nlevs,dm),mtemp(nlevs,dm),rhotot_fc(nlevs,dm),gradp_baro(nlevs,dm))
  allocate(umac_sum(nlevs,dm))
  allocate(umac_avg(nlevs,dm))
  allocate(eta(nlevs),kappa(nlevs))

  ! 1 component in 2D, 3 components in 3D
  allocate(eta_ed(nlevs,2*dm-3))
  allocate(Temp_ed(nlevs,2*dm-3))

  allocate(Epot_mass_fluxdiv(nlevs))
  allocate(charge_old(nlevs),charge_new(nlevs))
  allocate(permittivity(nlevs))
  allocate(grad_Epot_old(nlevs,dm),grad_Epot_new(nlevs,dm))
  allocate(Epot(nlevs))
  allocate(Epot_sum(nlevs)) 
  allocate(Epot_avg(nlevs)) 
  allocate(gradPhiApprox(nlevs,dm))

  if (analyze_cuts>0) then
     if(algorithm_type /= 6) call bl_error("To analyze fluxes on a cut you must use Boussinesq algorithm_type=6")
     allocate(total_mass_flux(nlevs,dm))
  end if   

  allocate(chem_rate(nlevs))

  ! build pmask
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  if (advection_type .eq. 0) then
     ng_s = 2 ! centered advection
  else if (advection_type .le. 3) then
     ng_s = 3 ! bilinear bds or unlimited quadratic bds
  else if (advection_type .eq. 4) then
     ng_s = 4 ! limited quadratic bds
  end if

  if (restart .ge. 0) then

     init_step = restart + 1


     ! build the ml_layout
     ! read in time and dt from checkpoint
     ! build and fill rho, rhotot, pi, umac, and umac_sum
     call initialize_from_restart(mla,time,dt,rho_old,rho_sum,rhotot_old,pi,umac,umac_sum,pmask, &
                                  Epot,Epot_sum,grad_Epot_old) 

     ! enabled a fixed_dt that is different from the previous run
     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     end if

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
        call multifab_build(rho_old(n)   ,mla%la(n),nspecies,ng_s)
        call multifab_build(rho_sum(n)   ,mla%la(n),nspecies,ng_s) 
        call setval(rho_sum(n),0.d0) 
        call multifab_build(rhotot_old(n),mla%la(n),1       ,ng_s)
        call multifab_build(pi(n)        ,mla%la(n),1       ,1)
        do i=1,dm
           call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
           call multifab_build_edge(umac_sum(n,i),mla%la(n),1,0,i)
           call setval(umac_sum(n,i),0.d0)
        end do
     end do

  end if

  ! data structures to help with reservoirs
  call build_bc_multifabs(mla)

  ! set grid spacing at each level
  allocate(dx(nlevs,MAX_SPACEDIM))
  dx(1,1:MAX_SPACEDIM) = (prob_hi(1:MAX_SPACEDIM)-prob_lo(1:MAX_SPACEDIM)) &
       / n_cells(1:MAX_SPACEDIM)
  ! check that the grid spacing is the same in each direction for the first dm dimensions
  if (dx(1,1) .ge. (1.d0 + 1.d-12)*dx(1,2) .or. &
      dx(1,1) .le. (1.d0 - 1.d-12)*dx(1,2)) then
     call bl_error('ERROR: main_driver.f90, we only support dx=dy')
  end if

  if (dm .eq. 3) then
     if (dx(1,1) .ge. (1.d0 + 1.d-12)*dx(1,3) .or. &
         dx(1,1) .le. (1.d0 - 1.d-12)*dx(1,3)) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dz')
     end if
  end if

  if (dm .ne. 2 .and. dm .ne. 3) then
     call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end if

  ! use refined dx for next level
  ! assume refinement ratio is the same in each direction
  ! we do this because dx is allocated over MAX_SPACEDIM dimensions,
  ! whereas mba is only allocated over dm dimensions
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,1)
  end do
  
  dx_saved(1:dm) = dx(1,1:dm) ! Store this in the module so anyone can access it

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  !
  ! next, for scalars, scal_bc_comp=dm+2
  ! there are 2*nspecies+3 "scalars"
  ! scal_bc_comp = rhotot
  ! scal_bc_comp+1 = c_i
  ! scal_bc_comp+nspecies+1 = molfrac or massfrac (dimensionless fractions)
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = Epot_bc_comp = electric potential
  !
  ! next, for transport coefficients, tran_bc_comp = scal_bc_comp+2*nspecies+3
  ! we say there is "one" transport coefficient
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  ! tran_bc_comp = diffusion coefficients (eta,kappa,chi)
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies+3, &
                     num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! build layouts for staggered multigrid solver and macproject within preconditioner
  call stag_mg_layout_build(mla)
  call mgt_macproj_precon_build(mla,dx,the_bc_tower)

  if (restart .lt. 0) then

     ! initialize rho and umac in valid region only
     call init_rho_and_umac(mla,rho_old,umac,dx,time,the_bc_tower%bc_tower_array)

     ! initialize pi, including ghost cells
     do n=1,nlevs
        call multifab_setval(pi(n),0.d0,all=.true.)
     end do

  end if

  ! compute rhotot from rho in VALID REGION
  call compute_rhotot(mla,rho_old,rhotot_old)

  ! fill rho and rhotot ghost cells
  call fill_rho_rhotot_ghost(mla,rho_old,rhotot_old,dx,the_bc_tower)

  do n=1,nlevs
     ! pressure ghost cells
     call multifab_fill_boundary(pi(n))
     call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                          dx_in=dx(n,:))
  end do

if(.false.) then  
  ! DONEV FIXME temporary debugging:
  w_mol=(/0.0014d0,0.0053d0,0.0037d0/) ! Random values
  if(nspecies==5) then
      ! 1=HCl, 2=NaOH, 3=NaCl,        4=H2O
      ! 1=Na+, 2=Cl-,  3=H+,   4=OH-, 5=H2O 
      ! w_Na = m_Na / (m_Na+m_OH) * w_NaOH + m_Na / (m_Na+m_Cl) * w_NaCl
      w_temp(1) = m_Na / (m_Na+m_OH) * w_mol(2) + m_Na / (m_Na+m_Cl) * w_mol(3)
      ! w_Cl = m_Cl / (m_H+m_Cl) * w_HCl + m_Cl / (m_Na+m_Cl) * w_NaCl
      w_temp(2) = m_Cl / (m_H+m_Cl) * w_mol(1) + m_Cl / (m_Na+m_Cl) * w_mol(3)
      ! w_H = m_H / (m_H+m_Cl) * w_HCl
      w_temp(3) = m_H / (m_H+m_Cl) * w_mol(1)
      ! w_OH = m_OH / (m_Na+m_OH) * w_NaOH
      w_temp(4) = m_OH / (m_Na+m_OH) * w_mol(2)
  else
      w_temp(1:3)=w_mol    
  end if    
  w_temp(nspecies) = 1.0d0 - sum(w_temp(1:nspecies))
  call compute_rho_eos(w_temp(1:nspecies)*rho0, rho_temp)
  write(*,*) " w_H2O=", w_temp(nspecies), " rho_eos=", rho_temp, " w=", w_temp(1:nspecies-1)
  stop
end if

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  do n=1,nlevs 
     call multifab_build(rho_new(n)          ,mla%la(n),nspecies,ng_s)
     call multifab_build(rho_avg(n)          ,mla%la(n),nspecies,ng_s) 
     call setval(rho_avg(n),0.d0) 
     call multifab_build(rhotot_new(n)       ,mla%la(n),1       ,ng_s) 
     call multifab_build(Temp(n)             ,mla%la(n),1       ,ng_s)
     call multifab_build(diff_mass_fluxdiv(n),mla%la(n),nspecies,0) 
     call multifab_build(eta(n)              ,mla%la(n),1       ,1)
     call multifab_build(kappa(n)            ,mla%la(n),1       ,1)
     do i=1,dm
        call multifab_build_edge(     mtemp(n,i),mla%la(n),1,0,i)
        call multifab_build_edge( rhotot_fc(n,i),mla%la(n),1,0,i)
        call multifab_build_edge(gradp_baro(n,i),mla%la(n),1,0,i)
        call multifab_build_edge(umac_avg(n,i),mla%la(n),1,0,i)
     end do
  end do

  do n=1,nlevs
     ! eta and Temp on nodes (2d) or edges (3d)
     if (dm .eq. 2) then
        call multifab_build_nodal(eta_ed(n,1),mla%la(n),1,0)
        call multifab_build_nodal(Temp_ed(n,1),mla%la(n),1,0)
     else
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(eta_ed(n,1),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,1),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,2),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,2),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,3),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,3),mla%la(n),1,0,nodal_temp)
     end if
  end do

  if (variance_coef_mass .ne. 0.d0) then
     do n=1,nlevs
        call multifab_build(stoch_mass_fluxdiv(n),mla%la(n),nspecies,0) 
        do i=1,dm
           call multifab_build_edge(stoch_mass_flux(n,i),mla%la(n),nspecies,0,i)
        end do
     end do
  end if

  if (use_charged_fluid) then
     do n=1,nlevs
        call multifab_build(charge_old(n)       ,mla%la(n),1,1)
        call multifab_build(charge_new(n)       ,mla%la(n),1,1)
        call multifab_build(permittivity(n)     ,mla%la(n),1,1)
        call multifab_build(Epot_mass_fluxdiv(n),mla%la(n),nspecies,0) 
        call multifab_build(Epot_avg(n)      ,mla%la(n),1       ,1) 

        if (restart .lt. 0) then ! only build Epot, Epot_sum, and grad_Epot_old below if a checkpoint file wasn't provided
           call multifab_build(Epot(n)             ,mla%la(n),1,1)
           call multifab_build(Epot_sum(n)      ,mla%la(n),1       ,1) 
        endif
        do i=1,dm
           if (restart .lt. 0) then 
              call multifab_build_edge(grad_Epot_old(n,i),mla%la(n),1,1,i)
           endif
           call multifab_build_edge(grad_Epot_new(n,i),mla%la(n),1,1,i)
           call multifab_build_edge(gradPhiApprox(n,i),mla%la(n),1,0,i)
        end do
     end do
  end if

  if (analyze_cuts>0) then
     do n=1,nlevs
        do i=1,dm
           call multifab_build_edge( total_mass_flux(n,i),mla%la(n),nspecies,0,i)
        end do
     end do
  end if   

  if (nreactions .gt. 0) then
     do n=1,nlevs
        call multifab_build(chem_rate(n),mla%la(n),nspecies,0)
     end do
  end if

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 2 .or. algorithm_type .eq. 5) then
     n_rngs_mass = 2
     n_rngs_mom  = 2
  else if (algorithm_type .eq. 6) then
     n_rngs_mass = 2
     n_rngs_mom  = 1
  else
     n_rngs_mass = 1
     n_rngs_mom  = 1
  end if
  call init_mass_stochastic(mla,n_rngs_mass)
  call init_m_stochastic(mla,n_rngs_mom)

  if (use_bl_rng) then
     ! save random state for writing checkpoint
     call bl_rng_copy_engine(rng_eng_diffusion_chk,rng_eng_diffusion)
  end if

  !=====================================================================
  ! Initialize values
  !=====================================================================

  if (use_charged_fluid) then

     ! set these to zero
     do n=1,nlevs
        call multifab_setval(charge_old(n),0.d0,all=.true.)
        call multifab_setval(charge_new(n),0.d0,all=.true.)

        if (restart .lt. 0) then
           call multifab_setval(Epot(n),0.d0,all=.true.)
           call multifab_setval(Epot_sum(n),0.d0,all=.true.) 
           call multifab_setval(Epot_avg(n),0.d0,all=.true.) 
        endif

        call multifab_setval(Epot_mass_fluxdiv(n),0.d0,all=.true.)
        do i=1,dm
           if (restart .lt. 0) then
              call multifab_setval(grad_Epot_old(n,i),0.d0,all=.true.)
           endif
           call multifab_setval(grad_Epot_new(n,i),0.d0,all=.true.)
           call multifab_setval(gradPhiApprox(n,i),0.d0,all=.true.)
        end do
     end do

     if(electroneutral) then
        call dot_with_z(mla,rho_old,charge_old,abs_z=.true.)
        max_charge_abs = multifab_norm_inf(charge_old(1)) ! This will be saved for reuse
        if(max_charge_abs<=0.0d0) max_charge_abs=1.0d0 ! Avoid division by zero
        call dot_with_z(mla,rho_old,charge_old,abs_z=.false.)
        max_charge = multifab_norm_inf(charge_old(1))
     end if 
     ! compute total charge
     call dot_with_z(mla,rho_old,charge_old)
     ! multiply by total volume (all 3 dimensions, even for 2D problems)
    
     if (use_charged_fluid) then ! NOTE: we are using rho = 1 here, so the below is a close approximation to debye length
        debye_len =sqrt(dielectric_const*k_B*T_init(1)/ &
           (rho0*sum(c_init(1,1:nspecies)*molmass(1:nspecies)*charge_per_mass(1:nspecies)**2))) 
        if (parallel_IOprocessor()) then 
           print*, 'Debye length $\lambda_D$ is approx: ', debye_len
        endif 
     endif

     total_charge = multifab_sum_c(charge_old(1),1,1)*product(dx(1,1:3))    
     if (parallel_IOProcessor().and.electroneutral) then
        print*,'Initial total charge',total_charge
        print*," Rel max charge=", max_charge/max_charge_abs
     else if (parallel_IOProcessor()) then
        print*,'Initial total charge',total_charge          
     end if

     ! compute permittivity
     if (dielectric_type .eq. 0) then
        do n=1,nlevs
           call multifab_setval(permittivity(n),dielectric_const,all=.true.)
        end do
     else
        call compute_permittivity(mla,permittivity,rho_old,rhotot_old,the_bc_tower)
     end if

  end if

  ! initialize Temp
  call init_Temp(Temp,dx,time,the_bc_tower%bc_tower_array)
  if (dm .eq. 2) then
     call average_cc_to_node(nlevs,Temp,Temp_ed(:,1),1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  else if (dm .eq. 3) then
     call average_cc_to_edge(nlevs,Temp,Temp_ed,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  end if

  if (barodiffusion_type .gt. 0) then

     ! this computes an initial guess at p using HSE
     call compute_HSE_pres(mla,rhotot_old,pi,dx,the_bc_tower)

     ! compute grad p for barodiffusion
     call compute_grad(mla,pi,gradp_baro,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

  end if

  ! initialize eta and kappa
  call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_old,rhotot_old,Temp,dx, &
                         the_bc_tower%bc_tower_array)

  ! now that we have eta, we can initialize the inhomogeneous velocity bc's
  ! set inhomogeneous velocity bc's to values supplied in inhomogeneous_bc_val
  call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time, &
                                 the_bc_tower%bc_tower_array)

  do n=1,nlevs
     do i=1,dm
        ! set normal velocity on physical domain boundaries
        call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                       the_bc_tower%bc_tower_array(n), &
                                       dx(n,:),vel_bc_n(n,:))
        ! set transverse velocity behind physical boundaries
        call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                    the_bc_tower%bc_tower_array(n), &
                                    dx(n,:),vel_bc_t(n,:))
        ! fill periodic and interior ghost cells
        call multifab_fill_boundary(umac(n,i))
        ! protect against roundoff issues and sync up 
        ! faces with the same physical location
         call multifab_internal_sync(umac(n,i))
     end do
  end do

  if (restart .lt. 0) then

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "Initial state:"  
        call sum_mass(rho_old, 0) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        call sum_kinetic(mla,umac,rhotot_fc)
        call eos_check(mla,rho_old)
     end if

     ! add initial momentum fluctuations
     ! do not call for overdamped codes since the steady Stokes solver will 
     ! wipe out the initial condition to solver tolerance
     if ((algorithm_type .ne. 2) .and. &
         (initial_variance_mom .ne. 0.d0)) then
        call add_m_fluctuations(mla,dx,initial_variance_mom, &
                                umac,rhotot_old,Temp,the_bc_tower)

        do n=1,nlevs
           do i=1,dm
              ! set normal velocity on physical domain boundaries
              call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                             the_bc_tower%bc_tower_array(n), &
                                             dx(n,:),vel_bc_n(n,:))
              ! set transverse velocity behind physical boundaries
              call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                          the_bc_tower%bc_tower_array(n), &
                                          dx(n,:),vel_bc_t(n,:))
              ! fill periodic and interior ghost cells
              call multifab_fill_boundary(umac(n,i))
           end do
        end do

        if (print_int .gt. 0) then
           if (parallel_IOProcessor()) write(*,*) "After adding momentum fluctuations:"
           call sum_mass(rho_old, 0) ! print out the total mass to check conservation
           ! compute rhotot on faces
           call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                   the_bc_tower%bc_tower_array)
           ! compute momentum
           call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
           call sum_momenta(mla,mtemp)
           call sum_kinetic(mla,umac,rhotot_fc)
           call eos_check(mla,rho_old)
        end if

     end if

     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        if (parallel_IOProcessor()) write(*,*) &
           "Estimated dt_diff=", dt_diffusive, " dt_adv=", dt, " setting dt=", min(dt,dt_diffusive)
        dt = min(dt,dt_diffusive)
     end if
     
  end if
  
  dt_saved = dt ! Store this in a module so anyone can access it

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
           
           ! analyze electric potential (will pass in as a scalar)
           if (use_charged_fluid) then
              nscal = 1
           else
              nscal = 0
           end if

           ! We will also pass temperature
           call initialize_hydro_grid(mla,rho_old,dt,dx,namelist_file=un, & 
                                      nspecies_in=nspecies, &
                                      nscal_in=nscal, &
                                      exclude_last_species_in=.false., &
                                      analyze_velocity=.true., &
                                      analyze_density=.true., &
                                      analyze_temperature=.true.) 
           
           close(unit=un)
           
           if ((algorithm_type .eq. 6) .and. all(density_weights(0:nspecies)==0.0d0)) then
              ! Make rho be rho_eos for HydroGrid analysis purposes
              ! This can be overwridden by specifying density_weights in the input file
              if(rho_eos_form==2) then ! Linearized EOS
                 density_weights(1:nspecies) =  - (rho0/rhobar(1:nspecies)-1.0d0)
                 if(.not.analyze_conserved) density_weights(1:nspecies)=rho0*density_weights(1:nspecies)
                 density_weights(0) = rho0 ! No inversion here
              else ! Nonlinear EOS
                 density_weights(1:nspecies) = 1.0d0 / rhobar(1:nspecies)
                 if(analyze_conserved) then
                    density_weights(0) = -rho0
                 else
                    density_weights(0) = -1.0
                 end if
              end if   
           end if
           
        else
        
           call bl_error('HydroGrid initialization requires a namelist in an input file')
           
        end if
     end if
  end if
  
  ! this routine is only called for all inertial simulations (both restart and non-restart)
  ! it does the following:
  ! 1. fill mass random numbers
  ! 2. computes mass fluxes and flux divergences
  ! if restarting, the subroutine ends; otherwise
  ! 3. perform an initial projection
  !
  ! overdamped schemes need to do 1. and 2. within the advance_timestep routine
  ! in principle, performing an initial projection for overdamped will change
  ! the reference state for the GMRES solver
  ! For overdamped the first ever solve cannot have a good reference state
  ! so in general there is the danger it will be less accurate than subsequent solves
  ! but I do not see how one can avoid that
  ! From this perspective it may be useful to keep initial_projection even in overdamped
  ! because different gmres tolerances may be needed in the first step than in the rest
  if (algorithm_type .ne. 2 .and. algorithm_type .ne. 6) then
     call initial_projection(mla,umac,rho_old,rhotot_old,gradp_baro, &
                             diff_mass_fluxdiv, &
                             stoch_mass_fluxdiv,stoch_mass_flux,chem_rate, &
                             Temp,eta,eta_ed,dt,dx,the_bc_tower, &
                             charge_old,grad_Epot_old,Epot,permittivity)
  end if

  if (restart .lt. 0) then

     !=====================================================================
     ! Process initial conditions (non-restart runs)
     !=====================================================================

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "After initial projection:"  
        call sum_mass(rho_old,0) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        call sum_kinetic(mla,umac,rhotot_fc)
        call eos_check(mla,rho_old)
     end if

     !=====================================================================
     ! Hydrogrid analysis and output for initial data
     !=====================================================================

     do n=1,nlevs
        do i=1,dm
           call multifab_copy_c(umac_avg(n,i),1,umac(n,i),1,1,0)
        end do
     end do

     ! write initial plotfile
     if (plot_int .gt. 0) then
        if (parallel_IOProcessor()) then
           write(*,*) 'writing initial plotfile 0'
        end if
        call write_plotfile(mla,rho_old,rho_avg,rhotot_old,Temp,umac,umac_avg,pi,Epot, & 
                            Epot_avg,grad_Epot_old,gradPhiApprox,0,dx,time)
     end if

     ! write initial checkpoint
     if (chk_int .gt. 0) then
        if (parallel_IOProcessor()) then
           write(*,*) 'writing initial checkpoint 0'
        end if
        call checkpoint_write(mla,rho_old,rho_sum,rhotot_old,pi,umac,umac_sum,Epot,Epot_sum,grad_Epot_old,time,dt,0) 
     end if
     
     if (stats_int .gt. 0) then
        ! write initial vertical and horizontal averages (hstat and vstat files)   
        if (use_charged_fluid) then
           call print_stats(mla,dx,0,time,umac=umac,rho=rho_old,temperature=Temp,scalars=Epot)
        else
           call print_stats(mla,dx,0,time,umac=umac,rho=rho_old,temperature=Temp)
        end if
     end if

     ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
     if (n_steps_skip .eq. 0) then

        ! Add this snapshot to the average in HydroGrid
        if (hydro_grid_int > 0) then
           call analyze_hydro_grid(mla,dt,dx,istep,umac=umac,rho=rho_old,temperature=Temp)
        end if

        if (hydro_grid_int > 0 .and. n_steps_save_stats > 0) then
           call save_hydro_grid(id=0, step=0)
        end if

     end if

  end if
  
  if(analyze_cuts>0) then
     call initialize_planar_cut()
  end if
  

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     if (fixed_dt .le. 0.d0) then
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        dt = min(dt,dt_diffusive)
     end if

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
             print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     runtime1 = parallel_wtime()

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! advance the solution by dt
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! notes: eta, eta_ed, and kappa could be built and initialized within the advance routines
     ! but for now we pass them around (it does save a few flops)
     ! diff/stoch_mass_fluxdiv could be built locally within the overdamped
     ! routine, but since we have them around anyway for inertial we pass them in
     if (algorithm_type .eq. 0) then
        ! algorithm_type=0: inertial
        call advance_timestep_inertial(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                       gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                       diff_mass_fluxdiv, stoch_mass_fluxdiv,stoch_mass_flux, &
                                       dx,dt,time,the_bc_tower,istep, &
                                       grad_Epot_old,grad_Epot_new, &
                                       charge_old,charge_new,Epot, &
                                       permittivity)
     else if (algorithm_type .eq. 2) then
        ! algorithm_type=2: overdamped with 2 RNG
        call advance_timestep_overdamped(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                         gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                         diff_mass_fluxdiv, stoch_mass_fluxdiv, stoch_mass_flux, &
                                         chem_rate, &
                                         dx,dt,time,the_bc_tower,istep)
     else if (algorithm_type .eq. 3) then
        ! algorithm_type=3: iterative implicit electrodiffusion
        call advance_timestep_iterative(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, &
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux, &
                                        dx,dt,time,the_bc_tower,istep, &
                                        grad_Epot_old,grad_Epot_new, &
                                        charge_old,charge_new,Epot, &
                                        permittivity,gradPhiApprox)
     else if (algorithm_type .eq. 4) then
        ! algorithm_type=4: boussinesq implicit electrodiffusion
        call advance_timestep_imp_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        Epot_mass_fluxdiv, &
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, stoch_mass_flux, &
                                        dx,dt,time,the_bc_tower,istep, &
                                        grad_Epot_old,grad_Epot_new, &
                                        charge_old,charge_new,Epot, &
                                        permittivity,gradPhiApprox)
     else if (algorithm_type .eq. 5) then
        ! algorithm_type=5: inertial midpoint
        call advance_timestep_inertial_midpoint(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                                gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                                diff_mass_fluxdiv, stoch_mass_fluxdiv,stoch_mass_flux, &
                                                chem_rate, &
                                                dx,dt,time,the_bc_tower,istep, &
                                                grad_Epot_old,grad_Epot_new, &
                                                charge_old,charge_new,Epot, &
                                                permittivity)
     else if (algorithm_type .eq. 6) then
        ! algorithm_type=6: boussinesq
        if (analyze_cuts>0) then
           call advance_timestep_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                    gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux, &
                                    chem_rate, dx,dt,time,the_bc_tower,istep, &
                                    grad_Epot_old,grad_Epot_new, &
                                    charge_old,charge_new,Epot, &
                                    permittivity, total_mass_flux)
           call planar_cut(mla, mf=total_mass_flux(1,analyze_cuts), dir=analyze_cuts, id=istep)      
        else
           call advance_timestep_bousq(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                    gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                    diff_mass_fluxdiv,stoch_mass_fluxdiv,stoch_mass_flux, &
                                    chem_rate, dx,dt,time,the_bc_tower,istep, &
                                    grad_Epot_old,grad_Epot_new, &
                                    charge_old,charge_new,Epot,permittivity)        
        end if                            
     else
        call bl_error("Error: invalid algorithm_type")
     end if

     ! for electroneutral make sure charge does not build up
     if (use_charged_fluid .and. electroneutral) then
        max_charge = multifab_norm_inf(charge_new(1))        
        if (max_charge > 10*epot_mg_rel_tol*max_charge_abs) then
           if (parallel_IOProcessor()) then
              print*,''
              print*,'WARNING: Max charge density exceeds rel_tol, rel max=', max_charge/max_charge_abs
              print*,''
           end if
        end if
     end if        

     ! for writing time-averaged umac, rho, and Epot to plotfile
     if (istep .gt. n_steps_skip) then
        ! Note: reset time avg quantities is only possible if reset_tavg_step >= n_steps_skip
        ! Also note: if reset is turned ON, and istep is between n_steps_skip and reset_tavg_step,
        !            then the code below will not track the average (since it will get reset anyways)
        !
        ! TL;DR: don't turn reset_tavg_vals ON unless you want to use it. 
        if (reset_tavg_vals .and. (reset_tavg_step.ge.n_steps_skip)) then
           if (istep.eq.reset_tavg_step) then 
              do n=1,nlevs
                 ! reset rho_sum, epot_sum, and umac_sum to be 0
                 call setval(rho_sum(n),0.d0)
                 call setval(Epot_sum(n),0.d0, all=.true.) 
                 do i=1,dm
                    call setval(umac_sum(n,i),0.d0) 
                 end do 
              end do 
           else if (istep .gt. reset_tavg_step) then
              ! do the normal averaging, starting from the reset step
              do n=1,nlevs
                 ! first do rho
                 call multifab_plus_plus_c(rho_sum(n),1,rho_new(n),1,nspecies,0)
                 call multifab_copy_c(rho_avg(n),1,rho_sum(n),1,nspecies,0)
                 call multifab_mult_mult_s_c(rho_avg(n),1,(1.d0/(istep-reset_tavg_step)),nspecies,0)
 
                 ! next do Epot
                 if (use_charged_fluid) then
                    call multifab_plus_plus_c(Epot_sum(n),1,Epot(n),1,1,0)
                    call multifab_copy_c(Epot_avg(n),1,Epot_sum(n),1,1,0)
                    call multifab_mult_mult_s_c(Epot_avg(n),1,(1.d0/(istep-reset_tavg_step)),1,0)
                 end if
         
                 ! lastly do umac
                 do i=1,dm
                    call multifab_plus_plus_c(umac_sum(n,i),1,umac(n,i),1,1,0)
                    call multifab_copy_c(umac_avg(n,i),1,umac_sum(n,i),1,1,0)
                    call multifab_mult_mult_s_c(umac_avg(n,i),1,(1.d0/(istep-reset_tavg_step)),1,0)
                 end do
              end do
           end if
        else  
           do n=1,nlevs
              ! first do rho
              call multifab_plus_plus_c(rho_sum(n),1,rho_new(n),1,nspecies,0)
              call multifab_copy_c(rho_avg(n),1,rho_sum(n),1,nspecies,0)
              call multifab_mult_mult_s_c(rho_avg(n),1,(1.d0/(istep-n_steps_skip)),nspecies,0)

              ! next do Epot
              if (use_charged_fluid) then
                 call multifab_plus_plus_c(Epot_sum(n),1,Epot(n),1,1,0)
                 call multifab_copy_c(Epot_avg(n),1,Epot_sum(n),1,1,0)
                 call multifab_mult_mult_s_c(Epot_avg(n),1,(1.d0/(istep-n_steps_skip)),1,0)
              end if
      
              ! lastly do umac
              do i=1,dm
                 call multifab_plus_plus_c(umac_sum(n,i),1,umac(n,i),1,1,0)
                 call multifab_copy_c(umac_avg(n,i),1,umac_sum(n,i),1,1,0)
                 call multifab_mult_mult_s_c(umac_avg(n,i),1,(1.d0/(istep-n_steps_skip)),1,0)
              end do
           end do
        end if
     end if

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

     if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) .or. &
          (istep .eq. max_step) ) then
        if (parallel_IOProcessor()) write(*,*) "After time step ", istep, " t=", time
        call sum_mass(rho_new, istep) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        call sum_kinetic(mla,umac,rhotot_fc)
        call eos_check(mla,rho_new)

        if (use_charged_fluid) then
           ! multiply by total volume (all 3 dimensions, even for 2D problems)
           total_charge = multifab_sum_c(charge_new(1),1,1)*product(dx(1,1:3))
           if(electroneutral) then
              max_charge = multifab_norm_inf(charge_new(1))
           end if
           if (parallel_IOProcessor().and.electroneutral) then
              print*,'Total charge',total_charge
              print*," Rel max charge=", max_charge/max_charge_abs
           else if (parallel_IOProcessor()) then
              print*,'Total charge',total_charge          
              if (print_debye_len) then
                 print*, 'Debye length is: ', debye_len
              endif
           end if
        end if


     end if

     ! write plotfile at specific intervals
     if (plot_int.gt.0 .and. ( (mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) ) then
        if (parallel_IOProcessor()) then
           write(*,*) 'writing plotfiles after timestep =', istep 
        end if
        call write_plotfile(mla,rho_new,rho_avg,rhotot_new,Temp,umac,umac_avg,pi,Epot, & 
                            Epot_avg,grad_Epot_new,gradPhiApprox,istep,dx,time)
     end if

     ! write checkpoint at specific intervals
     if ((chk_int.gt.0 .and. mod(istep,chk_int).eq.0)) then

        if (parallel_IOProcessor()) then
           write(*,*) 'writing checkpoint after timestep =', istep 
        end if
        call checkpoint_write(mla,rho_new,rho_sum,rhotot_new,pi,umac,umac_sum,Epot,Epot_sum,grad_Epot_new,time,dt,istep)
     end if

     ! print out projection (average) and variance
     if ( (stats_int > 0) .and. &
          (mod(istep,stats_int) .eq. 0) ) then

        if (stat_save_type.eq.0) then  ! compute vertical/horizontal averages of instantaneous fields
           if (use_charged_fluid) then 
              call print_stats(mla,dx,istep,time,umac=umac,rho=rho_new,temperature=Temp,scalars=Epot) 
           else 
              call print_stats(mla,dx,istep,time,umac=umac,rho=rho_new,temperature=Temp) 
           endif 
        else                           ! compute vertical/horizontal averages of time averaged fields
           if (istep.lt.n_steps_skip) then 
              call bl_error("If stat_save_type is \ne 0, you cannot call print_stats until n_steps_skip has passed.")
           end if 
           if (use_charged_fluid) then 
              call print_stats(mla,dx,istep,time,umac=umac_avg,rho=rho_avg,temperature=Temp,scalars=Epot_avg) 
           else 
              call print_stats(mla,dx,istep,time,umac=umac_avg,rho=rho_avg,temperature=Temp) 
           endif 
        end if  
        
     end if

     if (istep .ge. n_steps_skip) then

        ! Add this snapshot to the average in HydroGrid
        if ( (hydro_grid_int > 0) .and. &
             ( mod(istep,hydro_grid_int) .eq. 0 ) ) then
           call analyze_hydro_grid(mla,dt,dx,istep,umac=umac,rho=rho_new,temperature=Temp)
        end if

        if ( (hydro_grid_int > 0) .and. &
             (n_steps_save_stats > 0) .and. &
             ( mod(istep,n_steps_save_stats) .eq. 0 ) ) then
           call save_hydro_grid(id=istep/n_steps_save_stats, step=istep)            
        end if

     end if

     ! set old state to new state
     do n=1,nlevs
        call multifab_copy_c(   rho_old(n),1,   rho_new(n),1,nspecies,   rho_old(n)%ng)
        call multifab_copy_c(rhotot_old(n),1,rhotot_new(n),1       ,1,rhotot_old(n)%ng)
     end do

     if (use_charged_fluid) then
        do n=1,nlevs
           call multifab_copy_c(charge_old(n),1,charge_new(n),1,1,charge_old(n)%ng)
           do i=1,dm
              call multifab_copy_c(grad_Epot_old(n,i),1,grad_Epot_new(n,i),1,1, &
                                   grad_Epot_old(n,i)%ng)
           end do
        end do
     end if

  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  if((hydro_grid_int>0) .or. (stats_int>0)) then
     call finalize_hydro_grid()
  end if

  call destroy_bc_multifabs(mla)
  call destroy_mass_stochastic(mla)
  call destroy_m_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(rho_old(n))
     call multifab_destroy(rhotot_old(n))
     call multifab_destroy(pi(n))
     do i=1,dm
        call multifab_destroy(umac(n,i))
     end do
  end do

  do n=1,nlevs
     call multifab_destroy(rho_new(n))
     call multifab_destroy(rhotot_new(n))
     call multifab_destroy(Temp(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     call multifab_destroy(diff_mass_fluxdiv(n))
     do i=1,dm
        call multifab_destroy(mtemp(n,i))
        call multifab_destroy(rhotot_fc(n,i))
        call multifab_destroy(gradp_baro(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
        call multifab_destroy(Temp_ed(n,i))
     end do
  end do

  if (variance_coef_mass .ne. 0.d0) then
     do n=1,nlevs
        call multifab_destroy(stoch_mass_fluxdiv(n))
        do i=1,dm
           call multifab_destroy(stoch_mass_flux(n,i))
        end do
     end do
  end if

  if (use_charged_fluid) then
     do n=1,nlevs
        call multifab_destroy(charge_old(n))
        call multifab_destroy(charge_new(n))
        call multifab_destroy(permittivity(n))
        call multifab_destroy(Epot(n))
        call multifab_destroy(Epot_mass_fluxdiv(n))
        do i=1,dm
           call multifab_destroy(grad_Epot_old(n,i))
           call multifab_destroy(grad_Epot_new(n,i))
           call multifab_destroy(gradPhiApprox(n,i))
        end do
     end do
  end if

  if (analyze_cuts>0) then
     do n=1,nlevs
        do i=1,dm
           call multifab_destroy(total_mass_flux(n,i))
        end do
     end do
     call destroy_planar_cut()
  end if
  
  if (nreactions .gt. 0) then
     do n=1,nlevs
        call multifab_destroy(chem_rate(n))
     end do
  end if

  do n=1,nlevs 
     call multifab_destroy(rho_sum(n))      
     call multifab_destroy(rho_avg(n))
     if (use_charged_fluid) then
        call multifab_destroy(Epot_sum(n))
        call multifab_destroy(Epot_avg(n))
     endif
     do i=1,dm
        call multifab_destroy(umac_sum(n,i))
        call multifab_destroy(umac_avg(n,i))
     end do
  end do

  deallocate(lo,hi,pmask)
  deallocate(rho_old,rhotot_old,pi)
  deallocate(rho_sum, rho_avg) 
  deallocate(rho_new,rhotot_new)
  deallocate(Temp)
  deallocate(diff_mass_fluxdiv,stoch_mass_fluxdiv)
  deallocate(stoch_mass_flux)
  deallocate(umac,mtemp,rhotot_fc,gradp_baro)
  deallocate(eta,kappa)
  deallocate(eta_ed)
  deallocate(Temp_ed)

  deallocate(Epot_mass_fluxdiv)
  deallocate(charge_old,charge_new)
  deallocate(permittivity)
  deallocate(grad_Epot_old,grad_Epot_new)
  deallocate(Epot)
  deallocate(Epot_sum, Epot_avg) 
  deallocate(gradPhiApprox)

  if (analyze_cuts>0) deallocate(total_mass_flux)

  deallocate(chem_rate)

  deallocate(dx)

  if (use_bl_rng) then
     call rng_destroy()
  end if

  call stag_mg_layout_destroy()
  call mgt_macproj_precon_destroy()
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
