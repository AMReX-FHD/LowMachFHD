!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HydroGrid analysis for spectra/rendering and other serial routines
  ! Make a copy of all the data on one of the processors and then call the serial routines
  ! Important note: RESTARTS NOT IMPLEMENTED in analysis code!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module analyze_spectra_module

  use fabio_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use convert_module
  !use mk_stoch_force_module

  use HydroGridModule
  use HydroGridCInterface 

  use probin_module, only : hydro_grid_int, nscal, pmask, n_cells, &
      variance_coeff, analyze_conserved, project_dir, center_snapshots, &
      max_grid_projection, stats_int, prob_lo, prob_hi, n_steps_save_stats

  implicit none

  private
  public :: initialize_hydro_grid, analyze_hydro_grid, save_hydro_grid, finalize_hydro_grid, &
       print_stats

  ! Molecular parameters (not used at present since there is no internal energy variable)
  real(dp_t), save :: k_B_over_m=1.0_dp_t

  ! For statistical analysis of fluctuating fields
  type (HydroGrid), target, save :: grid, grid_2D, grid_1D
  
  ! We collect all the data on processor 1 for the analysis stuff due to the FFTs etc:
  integer, save :: nvar, ncells(3), ncells2D(3), ncells1D(3)
  type(multifab), save ::  s_serial, s_dir, s_projected, s_var ! Only one level here
  ! Must keep the layout around until the end!
  type(layout), save :: la_serial, la_dir, la_projected

contains   

  subroutine initialize_hydro_grid(mla,s_in,m_in,dt,dx,namelist_file)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_in(:)
    type(multifab) , intent(inout) :: m_in(:,:)
    real(dp_t)     , intent(inout) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: namelist_file ! Where to read the namelists from

    ! local
    type(box)  :: bx_serial, bx_dir, bx_projected
    type(boxarray)  :: bxa_serial, bxa_dir, bxa_projected
    integer, dimension(mla%dim) :: lo, hi
    integer :: nlevs, dm, pdim, max_grid_dir(3), max_grid_projected(3)
    ! The analysis codes always works in 3D
    real(dp_t) :: grid_dx(3)

    nlevs = mla%nlevel
    dm = mla%dim

    if(nlevs>1) call parallel_abort("HydroGrid analysis only implemented for a single level!")

    nvar = dm + nscal ! mx, my, [mz,] rho, rho*c (or ux, uy, [uz,] rho, c if primitive)

    ncells(1) = n_cells(1)
    ncells(2) = n_cells(2)
    ncells(3) = 1
    if(dm>2) ncells(3) = n_cells(3)

    pdim=abs(project_dir) ! Axes to project along
    if(pdim>0) then
       ncells2D=ncells
       ncells2D(pdim)=1   
       ncells1D=1
       ncells1D(pdim)=ncells(pdim)
    end if

    if ((stats_int > 0) .or. (project_dir > 0)) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! build s_dir (multifab with "tall skinny boxes")
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! It is hard to write a general algorithm for figuring out the best way to split into boxes
       ! So we make the user specify that via max_grid_projection
       select case(pdim)
       case(1)
          max_grid_dir=(/n_cells(1),max_grid_projection(1),max_grid_projection(2)/)
       case(2)
          max_grid_dir=(/max_grid_projection(1),n_cells(2),max_grid_projection(2)/)
       case(3)
          max_grid_dir=(/max_grid_projection(1),max_grid_projection(2),n_cells(3)/)
       case default
          call bl_error("project_dir must be between 1 and 3 if stats_int>0")
       end select

       bx_dir = s_in(1)%la%lap%pd                   ! set bx_dir to problem domain
       call boxarray_build_bx(bxa_dir, bx_dir)      ! build a boxarray containing only one box
       call boxarray_maxsize(bxa_dir, max_grid_dir) ! chop into tall skinny boxes
       if ( parallel_IOprocessor() ) then
          call print(bxa_dir, 'Analysis tall/skinny BoxArray')
          print*,''
       end if
       call layout_build_ba (la_dir, bxa_dir, bx_dir, pmask)
       call destroy(bxa_dir)
       call multifab_build(s_dir, la_dir, nvar, 0)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! build s_projected (multifab with reduced dimension that holds the projection of s_dir)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! first set lo and hi to the problem domain
       lo = 0
       hi(1) = n_cells(1)-1
       if (dm > 1) then   
          hi(2) = n_cells(2) - 1        
          if (dm > 2)  then
             hi(3) = n_cells(3) -1
          endif
       endif
       ! then reduce the dimensionality by 1
       hi(pdim)=lo(pdim)
       select case(pdim)
       case(1)
          max_grid_projected=(/1,max_grid_projection(1),max_grid_projection(2)/)
       case(2)
          max_grid_projected=(/max_grid_projection(1),1,max_grid_projection(2)/)
       case(3)
          max_grid_projected=(/max_grid_projection(1),max_grid_projection(2),1/)
       case default
          call bl_error("project_dir must be between 1 and 3")
       end select
       
       call box_build_2(bx_projected,lo,hi)                     ! set bx_projected to reduced dimension problem domain
       call boxarray_build_bx(bxa_projected,bx_projected)       ! build a boxarray containing only one box
       call boxarray_maxsize(bxa_projected, max_grid_projected) ! chop up boxes
       if ( parallel_IOprocessor() ) then
          call print(bxa_projected, 'Analysis projected BoxArray')
          print*,''
       end if
       ! force same mapping as la_dir
       call layout_build_ba(la_projected, bxa_projected, bx_projected, pmask, &
                            explicit_mapping=get_proc(la_dir))
       call destroy(bxa_projected)
       call multifab_build(s_projected, la_projected, nvar, 0)

       if (stats_int > 0) then

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_var (multifab with reduced dimension that holds the variance of s_dir)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! has same layout as s_projected
          call multifab_build(s_var, la_projected, nvar, 0)

       end if

    end if

    if (abs(hydro_grid_int)>0) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Build serialized MultiFabs for HydroGrid analysis
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (project_dir > 0) then ! Serialize projection analysis only

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_serial (serial version of s_projected)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! set bx_serial to reduced dimension problem domain
          call box_build_2(bx_serial,lo,hi)
          ! build a boxarray containing only one box
          call boxarray_build_bx(bxa_serial,bx_serial)
          if ( parallel_IOprocessor() ) then
             call print(bxa_serial, 'Analysis serial BoxArray')
             print*,''
          end if
          call layout_build_ba(la_serial, bxa_serial, bx_serial, pmask, &
               explicit_mapping=(/parallel_IOProcessorNode()/) )
          call destroy(bxa_serial)
          call multifab_build(s_serial, la_serial, nvar, 0)

       else ! Serialize the whole analysis

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_serial (serial version of s_in)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! set bx_serial to problem domain
          bx_serial = s_in(1)%la%lap%pd
          ! build a boxarray containing only one box
          call boxarray_build_bx(bxa_serial, bx_serial) 
          call layout_build_ba(la_serial, bxa_serial, bx_serial, pmask, &
               explicit_mapping=(/parallel_IOProcessorNode()/) )
          if ( parallel_IOprocessor() ) then
             call print(bxa_serial, 'Analysis serial BoxArray')
             print*,''
          end if
          call destroy(bxa_serial)       
          call multifab_build(s_serial, la_serial, nvar, 0)

       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Now initialize the analysis code itself:   
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(parallel_IOProcessor()) then
       
          if( analyze_conserved .and. (.not.center_snapshots) ) then
            call bl_warn("If analyze_conserved is set then center_snapshots should also be set")
          end if

          grid_dx=1.0_dp_t ! Default grid spacing
          grid_dx(1:dm)=dx(1,1:dm) 

          if(project_dir <= 0) then ! Perform analysis on the full 3D grid
            call createHydroAnalysis (grid, ncells, nSpecies = 2, &
               isSingleFluid = .true., nVelocityDimensions = dm, nPassiveScalars=0, &
               systemLength = ncells*grid_dx, heatCapacity = (/1.5_dp_t*k_B_over_m, 1.5_dp_t*k_B_over_m /), &
               timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
               structFactMultiplier = 1.0_dp_t/max(variance_coeff, epsilon(1.0_dp_t)) )
          else ! Create a fake object that will not be used, with grid size 1x1x1
            ! This way the namelist input files do not have to change     
            call createHydroAnalysis (grid, nCells=(/1,1,1/), nSpecies = 2, &
               isSingleFluid = .true., nVelocityDimensions = dm, nPassiveScalars=0, &
               systemLength = ncells*grid_dx, heatCapacity = (/1.5_dp_t*k_B_over_m, 1.5_dp_t*k_B_over_m /), &
               timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
               structFactMultiplier = 1.0_dp_t/max(variance_coeff, epsilon(1.0_dp_t)) )
          end if     

          if(project_dir/=0) then
             ! Also perform analysis on a projected grid (averaged along project_dir axes)
             call createHydroAnalysis (grid_2D, nCells=ncells2D, nSpecies=2, &
                  isSingleFluid = .true., nVelocityDimensions = dm, &
                  systemLength = nCells*grid_dx, &
                  heatCapacity = (/1.5_dp_t*k_B_over_m, 1.5_dp_t*k_B_over_m /), &
                  timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
                  structFactMultiplier = 1.0_dp_t/max(variance_coeff, epsilon(1.0_dp_t)) )
          end if

          if(project_dir/=0) then ! Also perform analysis on a 1D grid (along project_dir only)
             call createHydroAnalysis (grid_1D, nCells=ncells1D, nSpecies=2, &
                  isSingleFluid = .true., nVelocityDimensions = dm, &
                  systemLength = nCells*grid_dx, &
                  heatCapacity = (/1.5_dp_t*k_B_over_m, 1.5_dp_t*k_B_over_m /), &
                  timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
                  structFactMultiplier = 1.0_dp_t/max(variance_coeff, epsilon(1.0_dp_t)) )
          end if

       end if

    end if
      
  end subroutine

  subroutine save_hydro_grid(id, step) ! This also *resets* all counters
      integer, intent(in) :: id, step ! We can use either one to number files here
      
      if((.not.parallel_IOProcessor()) .or. (hydro_grid_int<=0)) return
      
      if(project_dir<=0) call writeToFiles(grid, id=step)
      if(project_dir/=0) then
         call writeToFiles(grid_2D, id=step)
      end if
      if(project_dir/=0) then
         call writeToFiles(grid_1D, id=step)
      end if
  end subroutine
  
  subroutine finalize_hydro_grid()

    if (stats_int > 0 .or. project_dir .gt. 0) then
       call multifab_destroy(s_dir)
       call multifab_destroy(s_projected)
       call destroy(la_dir)
       call destroy(la_projected)
       if (stats_int > 0) then
          call multifab_destroy(s_var)
       end if
    end if
    if (abs(hydro_grid_int)>0) then
       call multifab_destroy(s_serial)
       call destroy(la_serial)
    end if
    
    if((.not.parallel_IOProcessor())) return
    
    if(hydro_grid_int<0) then ! Only do clean-up

       call destroyHydroAnalysis(grid)
       if(project_dir/=0) then
          call destroyHydroAnalysis(grid_2D) 
          call destroyHydroAnalysis(grid_1D)
       end if
    
    else if(hydro_grid_int>0) then ! Save HydroGrid data and clean-up

       if(project_dir<=0) then
         if(n_steps_save_stats <= 0) call writeToFiles(grid)
       end if  
       call destroyHydroAnalysis(grid)
       if(project_dir/=0) then
          if(n_steps_save_stats <= 0) call writeToFiles(grid_2D)
          call destroyHydroAnalysis(grid_2D) 
       end if
       if(project_dir/=0) then
          if(n_steps_save_stats <= 0) call writeToFiles(grid_1D)
          call destroyHydroAnalysis(grid_1D)
       end if
    
    end if   

  end subroutine
  
  !mcai----------start------------------------
  subroutine analyze_hydro_grid(mla,s_in,m_in,umac,prim,dt,dx,step,custom_analysis)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_in(:)
    type(multifab) , intent(inout) :: m_in(:,:)
    type(multifab) , intent(inout) ::    umac(:,:)
    type(multifab) , intent(inout) ::    prim(:)
    real(dp_t)     , intent(inout) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    logical, intent(in) :: custom_analysis
  
    if(project_dir>0) then ! Only do 2D analysis
      !write(*,*) "Calling analyze_hydro_grid_parallel"
      call analyze_hydro_grid_parallel(mla,s_in,m_in,umac,prim,dt,dx,step)
    else ! Analyze the whole grid
      !write(*,*) "Calling analyze_hydro_grid_serial"
      call analyze_hydro_grid_serial(mla,s_in,m_in,umac,prim,dt,dx,step,custom_analysis)
    end if

  !mcai----------end-----------------------

  end subroutine
  
   ! We need to make centered velocities to use FFTs directly
   ! These have a grid of size (nx,ny,nz) and not (nx+1,ny,nz) etc. like the staggered velocities do
   subroutine StaggeredToCentered(mac_cc, mla,s_in,m_in,umac,prim)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in) :: s_in(:)
    type(multifab) , intent(in) :: m_in(:,:)
    type(multifab) , intent(in) ::    umac(:,:)
    type(multifab) , intent(in) ::    prim(:)

    type(multifab), intent(inout) :: mac_cc(mla%nlevel) 

    integer :: i, n, nlevs, dm, pdim

    nlevs = mla%nlevel
    dm = mla%dim
    
    if(analyze_conserved) then    
       if(center_snapshots) then
          ! Momentum on a face here is split in half between the two adjacent cells
          ! This conserves total momentum even if boundaries are present
          ! But it will also smooth the spectrum of the fluctuations!
          do i=1,dm
             call average_face_to_cc(mla,m_in(:,i),1,mac_cc,i)
          end do
       else
          ! Pretend that velocities are cell-centered instead of staggered
          ! This is not right but for periodic should be OK as it preserves the spectrum of fluctuations
          do i=1,dm
             call shift_face_to_cc(mla,m_in(:,i),1,mac_cc,i)
          end do       
       end if  
   else ! Use primitive variables
       if(center_snapshots) then
          ! Momentum on a face here is split in half between the two adjacent cells
          ! This conserves total momentum even if boundaries are present
          ! But it will also smooth the spectrum of the fluctuations!
          do i=1,dm
             call average_face_to_cc(mla,m_in(:,i),1,mac_cc,i)
          end do
          ! Now convert momentum to velocity
          do n=1,nlevs
             do i=1,dm
                call multifab_div_div_c(mac_cc(n), i, s_in(n), 1, 1, 0)
             end do
          end do
       else
          ! Pretend that velocities are cell-centered instead of staggered
          ! This is not right but for periodic should be OK as it preserves the spectrum of fluctuations
          do i=1,dm
             call shift_face_to_cc(mla,umac(:,i),1,mac_cc,i)
          end do       
       end if
   end if        
  
  end subroutine  
  
!mcai----------start----------------------------
  ! This analysis routine is serialized:
  ! All the data from all the boxes is collected onto one processor and then analyzed
  subroutine analyze_hydro_grid_serial(mla,s_in,m_in,umac,prim,dt,dx, step, custom_analysis)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_in(:)
    type(multifab) , intent(inout) :: m_in(:,:)
    type(multifab) , intent(inout) ::    umac(:,:)
    type(multifab) , intent(inout) ::    prim(:)
    real(dp_t)     , intent(inout) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    logical, intent(in) :: custom_analysis    

    integer :: i, ii, jj, n, nlevs, dm, pdim
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: variables
    real(kind=dp_t), dimension(:,:), allocatable :: variables_1D ! Projection of data along project_dir

    !mcai------also give the file name for snapshot    
    character(len=25) :: cTemp
    character(len=1024) :: plot_filenameBase     !local variable, file name for snapshots

    type(multifab) :: mac_cc(mla%nlevel)
    
    nlevs = mla%nlevel
    dm = mla%dim
    pdim=abs(project_dir)

    ! Do a global gather to collect the full grid into a single box:
    ! -------------------------
    do n=1,nlevs
       call multifab_build(mac_cc(n), mla%la(n), dm, 0)
    end do    
    call StaggeredToCentered(mac_cc, mla,s_in,m_in,umac,prim)
    call copy(s_serial,1,mac_cc(1),1,dm)
    do n=1,nlevs
       call multifab_destroy(mac_cc(n))
    end do
    
    ! Now gather density and concentration
    if(analyze_conserved) then ! Use rho and rho1
       call copy(s_serial,dm+1,s_in(1),1,nscal)
    else ! Use rho and c
       call copy(s_serial,dm+1,prim(1),1,nscal)
    end if   

    if(parallel_IOProcessor()) then  
       !write(*,*) "Calling updateHydroAnalysis on 3D, 2D and 1D grids"

       ! Get to the actual data
       variables  => dataptr(s_serial,1) ! Gets all of the components
       
       if(pdim>0) then

          allocate ( variables_1D(lbound(variables,dim=pdim) : ubound(variables,dim=pdim), &
            lbound(variables,dim=4) : ubound(variables,dim=4)) )

          ! Average the data in the hyperplanes perpendicular to pdim:
          do jj =  lbound(variables,dim=4), ubound(variables,dim=4)
             do ii = lbound(variables,dim=pdim), ubound(variables,dim=pdim)     
                if (pdim .eq. 1) then
                   variables_1D(ii,jj) = sum(variables(ii,:,:,jj))
                else if (pdim .eq. 2) then
                   variables_1D(ii,jj) = sum(variables(:,ii,:,jj))
                else if (pdim .eq. 3) then
                   variables_1D(ii,jj) = sum(variables(:,:,ii,jj))
                end if
             end do
          end do
          
       end if   

       if(custom_analysis) then ! Do some serialized analysis for mixing experiments
       
         ! define the name of the statfile that will be written
         write( cTemp,'(i6.6)' ) step
         plot_filenameBase="projectedHydroGrid" // trim(adjustl(cTemp))
         if(project_dir/=0) then
            call resetHydroAnalysis(grid_2D)
         end if
         if(analyze_conserved) then
            if(nvar>=dm+2) then
               if(project_dir/=0) then
                  call projectHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,dm+2)/variables(:,:,:,dm+1), filename=plot_filenameBase, &
                    grid_2D=grid_2D )
               else     
                  call projectHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,dm+2)/variables(:,:,:,dm+1), filename=plot_filenameBase )
               end if     
               call writeHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,dm+2)/variables(:,:,:,dm+1), filename=plot_filenameBase )
            end if
         else
            if(project_dir/=0) then
               call projectHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,2+dm:nvar), filename=plot_filenameBase, &
                    grid_2D=grid_2D )
            else     
               call projectHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,2+dm:nvar), filename=plot_filenameBase)
            end if        
            call writeHydroGridMixture (grid, density=variables(:,:,:,dm+1), &
                    concentration=variables(:,:,:,2+dm:nvar), filename= plot_filenameBase)

         end if
         if(project_dir/=0) then
            call writeToFiles(grid_2D, id=step)
         end if
       
       else if(analyze_conserved) then
       
         ! This accepts the species densities as input and calculates 
         ! concentrations and velocity internally
	 ! It is not, however, quite right for a staggered grid since the analysis code
         ! does not use the staggered density for calculating velocities!
         call updateHydroAnalysisConserved (grid, density=variables(:,:,:,dm+1:nvar), &
                                           current=variables(:,:,:,1:dm)) 
         if(project_dir/=0) then ! Use momentum current
            call updateHydroAnalysisConserved (grid_2D, &
               density=SUM(variables(:,:,:,dm+1:nvar), DIM=pdim)/grid%nCells(pdim), &
               current=SUM(variables(:,:,:,1:dm), DIM=pdim)/grid%nCells(pdim))
            call updateHydroAnalysisConserved (grid_1D, density=variables_1D(:,dm+1:nvar), &
                                            current=variables_1D(:,1:dm))
         end if
       else ! Use concentration instead of density
         call updateHydroAnalysisPrimitive (grid, velocity=variables(:,:,:,1:dm), &
             density=variables(:,:,:,dm+1), concentration=variables(:,:,:,2+dm:nvar))

         if(project_dir/=0) then ! Use velocities
            call updateHydroAnalysisPrimitive (grid_2D, &
               density=SUM(variables(:,:,:,dm+1), DIM=pdim)/grid%nCells(pdim), &
               velocity=SUM(variables(:,:,:,1:dm), DIM=pdim)/grid%nCells(pdim), &
               concentration=SUM(variables(:,:,:,2+dm:nvar), DIM=pdim)/grid%nCells(pdim) )
            call updateHydroAnalysisPrimitive (grid_1D, velocity=variables_1D(:,1:dm), &
               density=variables_1D(:,dm+1), concentration=variables_1D(:,2+dm:nvar))
         end if
       end if
               
       if(pdim>0) deallocate(variables_1D)                                
    end if   

  end subroutine

  !mcai----------start------------------------
  ! This routine first projects onto the project_dir direction and the perpendicular 
  ! hyperpane in parallel.  Then it serializes the analysis of those projections
  subroutine analyze_hydro_grid_parallel(mla,s_in,m_in,umac,prim,dt,dx,step)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_in(:)
    type(multifab) , intent(inout) :: m_in(:,:)
    type(multifab) , intent(inout) ::    umac(:,:)
    type(multifab) , intent(inout) ::    prim(:)
    real(dp_t)     , intent(inout) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step     ! at what step?, optional
   !mcai----------end------------------------
   
    integer nlevs,dm,pdim,i,n,qdim,qdim1,qdim2
    integer lo(mla%dim),hi(mla%dim)
    integer ii,jj,kk

    real(kind=dp_t), pointer, dimension(:,:,:,:) :: variables, sdp, spp, svp
    real(kind=dp_t), dimension(:,:), allocatable :: variables_1D, variables_1D_proc
      ! Projection of data along project_dir

    type(multifab) :: mac_cc(mla%nlevel)

    nlevs = mla%nlevel
    dm = mla%dim
    pdim  = abs(project_dir)

    ! qdim is orthogonal to pdim
    if (pdim .eq. 1) then
       qdim  = 2
       qdim2 = 3
    else if (pdim .eq. 2) then
       qdim  = 1
       qdim2 = 3
    else if (pdim .eq. 3) then
       qdim  = 1
       qdim2 = 2
    end if

    ! Re-distribute the full grid into a grid of "tall skinny boxes"
    ! These boxes are not distributed along project_dim so we can do local analysis easily
    ! -------------------------
    do n=1,nlevs
       call multifab_build(mac_cc(n), mla%la(n), dm, 0)
    end do    
    call StaggeredToCentered(mac_cc, mla,s_in,m_in,umac,prim)
    call copy(s_dir,1,mac_cc(1),1,dm)
    do n=1,nlevs
       call multifab_destroy(mac_cc(n))
    end do
    
    ! Now gather density and concentration
    if(analyze_conserved) then ! Use rho and rho1
       call copy(s_dir,dm+1,s_in(1),1,nscal)
    else ! Use rho and c
       call copy(s_dir,dm+1,prim(1),1,nscal)
    end if   

    ! Compute s_projected as the average along project_dim
    ! -------------------------
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       spp => dataptr(s_projected, i)
       lo = lwb(get_box(s_dir, i))
       hi = upb(get_box(s_dir, i))
       ! put sum ( s_dir / ncell_pdim ) into s_projected
       if (pdim .eq. 1) then
          spp(0,:,:,:)=SUM( sdp, DIM=1 ) / dble(hi(1)-lo(1)+1)
       else if (pdim .eq. 2) then
          spp(:,0,:,:)=SUM( sdp, DIM=2 ) / dble(hi(2)-lo(2)+1)
       else if (pdim .eq. 3) then
          spp(:,:,0,:)=SUM( sdp, DIM=3 ) / dble(hi(3)-lo(3)+1)
       end if
    end do

    ! Now collect the projected data from s_projected into a single box in s_serial 
    call copy(s_serial,1,s_projected,1,nvar)
    
    ! We also want to collect 1D data for the projection onto the direction pdim
    ! We need to do the averaging for each of the boxes in s_dir and then do an mpi_reduction
    allocate ( variables_1D(1:ncells(pdim), nvar), variables_1D_proc(1:ncells(pdim), nvar) )    
    call average_1D(variables_1D_proc) ! Sum over cells in each box    
    ! sum reduction: Note that dividing by the number of cells is already done in average_1D
    do n=1,nvar
       call parallel_reduce(variables_1D(:,n), variables_1D_proc(:,n), MPI_SUM, &
                            proc=parallel_IOProcessorNode())
    end do

    if(parallel_IOProcessor()) then  
       !write(*,*) "Calling updateHydroAnalysis on 2D and 1D grids"

       ! Get to the actual data
       variables  => dataptr(s_serial,1) ! Gets all of the components
       
       if(analyze_conserved) then ! Use conserved variables
         call updateHydroAnalysisConserved (grid_2D, density=variables(:,:,:,dm+1:nvar), &
                                            current=variables(:,:,:,1:dm))
         call updateHydroAnalysisConserved (grid_1D, density=variables_1D(:,dm+1:nvar), &
                                            current=variables_1D(:,1:dm))
       else ! Use primitive variables
         call updateHydroAnalysisPrimitive (grid_2D, velocity=variables(:,:,:,1:dm), &
             density=variables(:,:,:,dm+1), concentration=variables(:,:,:,2+dm:nvar))
         call updateHydroAnalysisPrimitive (grid_1D, velocity=variables_1D(:,1:dm), &
             density=variables_1D(:,dm+1), concentration=variables_1D(:,2+dm:nvar))
       end if
                                      
    end if   

    deallocate(variables_1D)

  contains
  
    subroutine average_1D(variables_1D)
    
      real(dp_t), dimension(:,:) :: variables_1D
    
      integer :: ii, jj, kk
      integer :: lo(3), hi(3) ! It is better to make this dimension independent this way

      variables_1D = 0.0_dp_t
      
      do i=1,nfabs(s_dir)
         sdp => dataptr(s_dir, i)
         lo=1; lo(1:dm) = lwb(get_box(s_dir, i))
         hi=1; hi(1:dm) = upb(get_box(s_dir, i))
         select case(pdim)
         case (1)
            do kk=lo(3),hi(3)
               do jj=lo(2),hi(2)
                  variables_1D = variables_1D + sdp(:,jj,kk,1:nvar)
               end do
            end do
         case (2)
            do kk=lo(3),hi(3)
               do ii=lo(1),hi(1)
                  variables_1D = variables_1D + sdp(ii,:,kk,1:nvar)
               end do
            end do
         case (3)
            do jj=lo(2),hi(2)
               do ii=lo(1),hi(1)
                  variables_1D = variables_1D + sdp(ii,jj,:,1:nvar)
               end do
            end do
         end select
      end do

       ! Divide by number of cells to get average
       ! Note that there is no need to divide again after the MPI reduction
       lo=1; lo(1:dm) = lwb(s_dir%la%lap%pd)
       hi=1; hi(1:dm) = upb(s_dir%la%lap%pd)
       variables_1D = variables_1D / ( (hi(qdim)-lo(qdim)+1)*(hi(qdim2)-lo(qdim2)+1) )
    
    end subroutine 

  end subroutine

  subroutine print_stats(mla,s_in,m_in,umac,prim,dx,step,time)

    ! This subroutine writes out the mean and variance of all variables using
    ! different averaging techniques.
    !
    ! We work with m_in and s_in (if analyze_conserved=T) or umac and prim (=F).
    !
    ! First, it writes out "vertical" averages, as defined by the
    ! pdim=abs(project_dir) direction.
    ! Thus, in 3D, it writes out a "2D" plotfile (a 3D plotfile with ncell=1
    ! in the pdim direction) called "vstatXXXXXX".
    ! In 2D, it writes out a 1D text file, also called "vstatXXXXX"
    !
    ! Next, it writes out "horizontal" averages.
    ! Thus, in both 2D and 3D, it writes out a 1D text file called "hstatXXXXXX"

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s_in(:)
    type(multifab) , intent(in   ) :: m_in(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    real(kind=dp_t), intent(in   ) :: time

    ! local
    type(multifab) :: mac_cc(mla%nlevel)

    integer nlevs,dm,pdim,qdim,qdim2,i,n
    integer lo(3),hi(3)
    integer ii,jj,kk

    ! pointers to access s_dir, s_projected, and s_var multifabs
    real(kind=dp_t), pointer :: sdp(:,:,:,:)
    real(kind=dp_t), pointer :: spp(:,:,:,:)
    real(kind=dp_t), pointer :: svp(:,:,:,:)

    real(kind=dp_t), allocatable :: stats_1d(:,:), stats_1d_proc(:,:)

    character(len=20) :: plotfile_name
    character(len=20) :: variable_names(2*nvar)
    type(multifab) :: plotdata(1)
    integer :: rr(0)
    type(box) :: bx

    nlevs = mla%nlevel
    dm    = mla%dim
    pdim  = abs(project_dir)

    ! qdim is orthogonal to pdim
    if (pdim .eq. 1) then
       qdim  = 2
       qdim2 = 3
    else if (pdim .eq. 2) then
       qdim  = 1
       qdim2 = 3
    else if (pdim .eq. 3) then
       qdim  = 1
       qdim2 = 2
    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    ! COMPUTE AND WRITE OUT VERTICAL AVERAGE/VARIANCE
    !!!!!!!!!!!!!!!!!!!!!!!

    ! copy m_in and s_in (if analyze_conserved=T) or umac and prim (=F)
    ! into the tall skinny s_dir
    do n=1,nlevs
       call multifab_build(mac_cc(n), mla%la(n), dm, 0)
    end do    
    call StaggeredToCentered(mac_cc, mla,s_in,m_in,umac,prim)
    call copy(s_dir,1,mac_cc(1),1,dm)
    do n=1,nlevs
       call multifab_destroy(mac_cc(n))
    end do
    if(analyze_conserved) then
       call copy(s_dir,dm+1,s_in(1),1,nscal)
    else
       call copy(s_dir,dm+1,prim(1),1,nscal)
    end if

    ! Compute s_projected (average) and s_var (variance)
    ! -------------------------
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       spp => dataptr(s_projected, i)
       svp => dataptr(s_var, i)
       lo(1:dm) = lwb(get_box(s_dir, i))
       hi(1:dm) = upb(get_box(s_dir, i))
       ! first put average, <x>, into s_projected
       ! then put variance, < (x - <x> )^2 > = <x^2> - <x>^2, into s_var
       select case (pdim)
       case (1)
          spp(0,:,:,:)=SUM( sdp, DIM=1 ) / dble(hi(1)-lo(1)+1)
          svp(0,:,:,:)=SUM( sdp**2, DIM=1 ) / dble(hi(1)-lo(1)+1) - spp(0,:,:,:)**2
       case (2)
          spp(:,0,:,:)=SUM( sdp, DIM=2 ) / dble(hi(2)-lo(2)+1)
          svp(:,0,:,:)=SUM( sdp**2, DIM=2 ) / dble(hi(2)-lo(2)+1) - spp(:,0,:,:)**2
       case (3)
          spp(:,:,0,:)=SUM( sdp, DIM=3 ) / dble(hi(3)-lo(3)+1)
          svp(:,:,0,:)=SUM( sdp**2, DIM=3 ) / dble(hi(3)-lo(3)+1) - spp(:,:,0,:)**2
       end select
    end do

    ! For 2D simulations, write the vertical average/variance to a text file
    if (dm .eq. 2) then
       
       ! collect the vertical average/variance in a 1D
       ! array so we can write it out in proper order

       bx = s_projected%la%lap%pd
       lo(1:dm) = lwb(bx)
       hi(1:dm) = upb(bx)

       ! components 1:nvar        will hold the average
       ! components nvar+1:2*nvar will hold the variance
       allocate(stats_1d_proc(lo(qdim):hi(qdim),2*nvar))
       allocate(stats_1d     (lo(qdim):hi(qdim),2*nvar))

       stats_1d_proc = 0.d0
       stats_1d      = 0.d0

       ! collect average and variance
       do i=1,nfabs(s_projected)
          spp => dataptr(s_projected, i)
          svp => dataptr(s_var, i)
          lo(1:dm) = lwb(get_box(s_projected, i))
          hi(1:dm) = upb(get_box(s_projected, i))
          select case(pdim)
          case (1)
             stats_1d_proc(lo(2):hi(2),1:nvar)        = spp(0,lo(2):hi(2),1,1:nvar)
             stats_1d_proc(lo(2):hi(2),nvar+1:2*nvar) = svp(0,lo(2):hi(2),1,1:nvar)
          case (2)
             stats_1d_proc(lo(1):hi(1),1:nvar)        = spp(lo(1):hi(1),0,1,1:nvar)
             stats_1d_proc(lo(1):hi(1),nvar+1:2*nvar) = svp(lo(1):hi(1),0,1,1:nvar)
          end select
       end do

       ! sum reduction
       do n=1,2*nvar
          call parallel_reduce(stats_1d(:,n), stats_1d_proc(:,n), MPI_SUM, &
                               proc=parallel_IOProcessorNode())
       end do

       if ( parallel_IOProcessor() ) then

          ! define the name of the statfile that will be written
          write(unit=plotfile_name,fmt='("vstat",i6.6)') step
          write(*,'(2A)') "Saving vSTAT FILEs to file ", trim(plotfile_name)
          write(*,*)

          lo(1:dm) = lwb(bx)
          hi(1:dm) = upb(bx)
          open(1000, file=trim(plotfile_name), status = "unknown", action = "write")

          if (analyze_conserved) then
             write(1000,'(A)') "# mx_avg my_avg rho_avg rho*c_avg mx_var my_var rho_var rho*c_var"
          else
             write(1000,'(A)') "# umac_avg vmac_avg rho_avg c_avg umac_var vmac_var rho_var c_var"
          end if
          do i=lo(qdim),hi(qdim)
             write(1000,'(1000(g17.9))') prob_lo(qdim) + (i+0.5d0)*dx(1,qdim), &
                  stats_1d(i,:)
          end do
          
          close(1000)
       end if

       deallocate(stats_1d,stats_1d_proc)
   
    ! For 3D simulations, write the vertical average/variance to a plotfile
    else if (dm .eq. 3) then

       ! components 1:nvar        will hold the average
       ! components nvar+1:2*nvar will hold the variance
       call multifab_build(plotdata(1),s_projected%la,2*nvar,0)

       ! copy s_projected and s_var into plotdata
       call multifab_copy_c(plotdata(1),1     ,s_projected,1,nvar)
       call multifab_copy_c(plotdata(1),nvar+1,s_var      ,1,nvar)

       if (analyze_conserved) then
          variable_names(1) = "mx_avg"
          variable_names(2) = "my_avg"
          variable_names(3) = "mz_avg"
          variable_names(4) = "rho_avg"
          variable_names(5) = "rho*c_avg"
          variable_names(6) = "mx_var"
          variable_names(7) = "my_var"
          variable_names(8) = "mz_var"
          variable_names(9) = "rho_var"
          variable_names(10) = "rho*c_var"
       else
          variable_names(1) = "umac_avg"
          variable_names(2) = "vmac_avg"
          variable_names(3) = "wmac_avg"
          variable_names(4) = "rho_avg"
          variable_names(5) = "c_avg"
          variable_names(6) = "umac_var"
          variable_names(7) = "vmac_var"
          variable_names(8) = "wmac_var"
          variable_names(9) = "rho_var"
          variable_names(10) = "c_var"
       end if

       ! define the name of the plotfile that will be written
       write(unit=plotfile_name,fmt='("vstat",i6.6)') step
       if ( parallel_IOProcessor() ) then
          write(*,'(2A)') "Saving vSTAT FILEs to directory ", trim(plotfile_name)
          write(*,*)
       end if

       ! write the plotfile
       call fabio_ml_multifab_write_d(plotdata, rr, plotfile_name, variable_names, &
                                      plotdata(1)%la%lap%pd, prob_lo, prob_hi, &
                                      time, dx(1,:))

       call multifab_destroy(plotdata(1))

    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    ! COMPUTE AND WRITE OUT HORIZONTAL AVERAGE/VARIANCE
    !!!!!!!!!!!!!!!!!!!!!!!

    bx = s_dir%la%lap%pd
    lo=1; lo(1:dm) = lwb(bx)
    hi=1; hi(1:dm) = upb(bx)

    ! will hold average (components 1:nvar) and variance (nvar+1:2*nvar)
    ! as a function of pdim
    allocate(stats_1d_proc(lo(pdim):hi(pdim),2*nvar))
    allocate(stats_1d     (lo(pdim):hi(pdim),2*nvar))

    stats_1d_proc = 0.d0
    stats_1d      = 0.d0

    ! put sum of x into average 
    ! put sum of x^2 and variance
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       lo=1; lo(1:dm) = lwb(get_box(s_dir, i))
       hi=1; hi(1:dm) = upb(get_box(s_dir, i))
       select case(pdim)
       case (1)
          do kk=lo(3),hi(3)
             do jj=lo(2),hi(2)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(:,jj,kk,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(:,jj,kk,1:nvar)**2
             end do
          end do
       case (2)
          do kk=lo(3),hi(3)
             do ii=lo(1),hi(1)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(ii,:,kk,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(ii,:,kk,1:nvar)**2
             end do
          end do
       case (3)
          do jj=lo(2),hi(2)
             do ii=lo(1),hi(1)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(ii,jj,:,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(ii,jj,:,1:nvar)**2
             end do
          end do
       end select
    end do

    ! sum reduction
    do n=1,2*nvar
       call parallel_reduce(stats_1d(:,n), stats_1d_proc(:,n), MPI_SUM, &
                            proc=parallel_IOProcessorNode())
    end do

    if ( parallel_IOProcessor() ) then

       ! divide by number of cells, so now average contains <x> and variance contains <x^2>
       lo=1; lo(1:dm) = lwb(bx)
       hi=1; hi(1:dm) = upb(bx)
       stats_1d = stats_1d / ( (hi(qdim)-lo(qdim)+1)*(hi(qdim2)-lo(qdim2)+1) )

       ! calculate variance = <x^2> - <x>^2
       stats_1d(:,nvar+1:2*nvar) = stats_1d(:,nvar+1:2*nvar) - stats_1d(:,1:nvar)**2

       ! define the name of the statfile that will be written
       write(unit=plotfile_name,fmt='("hstat",i6.6)') step
       write(*,'(2A)') "Saving hSTAT FILEs to file ", trim(plotfile_name)
       write(*,*)

       lo(1:dm) = lwb(bx)
       hi(1:dm) = upb(bx)
       open(1000, file=trim(plotfile_name), status = "unknown", action = "write")

       if (analyze_conserved) then
          if (dm .eq. 2) then
             write(1000,'(A)') "# mx my rho rho*c mx_var my_var rho_var rho*c_var"
          else if (dm .eq. 3) then
             write(1000,'(A)') "# mx my mz rho rho*c mx_var my_var mz_var rho_var rho*c_var"
          end if
       else
          if (dm .eq. 2) then
             write(1000,'(A)') &
             "# umac vmac rho c umac_var vmac_var rho_var c_var"
          else if (dm .eq. 3) then
             write(1000,'(A)') &
             "# umac vmac wmac rho c umac_var vmac_var wmac_var rho_var c_var"
          end if
       end if
       do i=lo(pdim),hi(pdim)
          write(1000,'(1000(g17.9))') prob_lo(pdim) + (i+0.5d0)*dx(1,pdim), stats_1d(i,:)
       end do

       close(1000)

    end if

  end subroutine print_stats


end module
