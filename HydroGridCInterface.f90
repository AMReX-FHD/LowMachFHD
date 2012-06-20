MODULE HydroGridCInterface ! Interface to my HydroGrid module
   USE, INTRINSIC :: ISO_C_BINDING
   USE HydroGridModule
   IMPLICIT NONE
   PUBLIC

   INTEGER, PARAMETER :: wp = c_double
   type(HydroGrid), save, private :: grid, grid_2D 
!   type(HydroGrid), save :: grid, grid_2D
      ! I use a global variable here but one could use C_F_POINTER
   logical, public, save :: project_2D=.false.
   
   character(kind=c_char, len=1024) :: input_file="hydroGridOptions.nml"

CONTAINS   

! The C interface to this routine is
! void setHydroInputFile_C(char * filename)
subroutine setHydroInputFile_C(filename) BIND(C,NAME="setHydroInputFile_C")
   character(kind=c_char), dimension(*), intent(in) :: filename

   integer :: i
   
   input_file=""
   do i=1, len(input_file)
      if(filename(i)==C_NULL_CHAR) exit
      input_file(i:i)=filename(i)
   end do
   write(*,*) "Will read Fortran namelist from file ", trim(input_file)
   
end subroutine

! void createHydroAnalysis_C (int nCells[3], int nSpecies, int nVelocityDimensions, \
!     int isSingleFluid, double systemLength[3], double heatCapacity[], double timestep, \
!     int nPassiveScalars, double structFactMultiplier, int project2D)
! 
! This will read input from the file hydroGridOptions.nml
!
subroutine createHydroAnalysis_C (nCells, nSpecies, nVelocityDimensions, isSingleFluid, &
      systemLength, heatCapacity, timestep, &
      nPassiveScalars, structFactMultiplier, project2D) &
      BIND(C,NAME="createHydroAnalysis_C")
   integer(c_int), intent(in) :: nCells(nMaxDims)
   integer(c_int), value :: nSpecies, nVelocityDimensions
   real (wp), intent(in) :: systemLength(nMaxDims), heatCapacity(nSpecies)
   real (wp), value :: timestep
   integer(c_int), value :: isSingleFluid, nPassiveScalars, project2D
   real (wp), value :: structFactMultiplier  
   
   integer :: nameListFile
   
   nameListFile = 114
   open (nameListFile, file = trim(input_file), status="old", action="read")
   
   !write(*,*) "PROJECT_Y=", project2D   
   call createHydroAnalysis (grid, nCells, nSpecies, nVelocityDimensions, &
      merge(.false.,.true.,isSingleFluid==0), systemLength, heatCapacity, timestep, &
      fileUnit=nameListFile, nPassiveScalars=nPassiveScalars, structFactMultiplier=structFactMultiplier)
   
   if(project2D>0) then ! Also do analysis for vertically-averaged fields
      project_2D=.true.
      call createHydroAnalysis (grid_2D, (/nCells(1), nCells(3), 1/), nSpecies, nVelocityDimensions, &
         merge(.false.,.true.,isSingleFluid==0), &
         (/systemLength(1),systemLength(3),systemLength(2)/), heatCapacity, timestep, &
         fileUnit=nameListFile, nPassiveScalars=nPassiveScalars, structFactMultiplier=structFactMultiplier)
   end if
   
   if(nameListFile>0) then
      close(nameListFile)
   end if
      
end subroutine

! void destroyHydroAnalysis_C ()
subroutine destroyHydroAnalysis_C () BIND(C,NAME="destroyHydroAnalysis_C")
   call destroyHydroAnalysis (grid)
   if(project_2D) call destroyHydroAnalysis (grid_2D)
end subroutine

! void resetHydroAnalysis_C ()
subroutine resetHydroAnalysis_C () BIND(C,NAME="resetHydroAnalysis_C")
   call resetHydroAnalysis (grid)
   if(project_2D) call resetHydroAnalysis (grid_2D)
end subroutine

! void updateHydroAnalysisIsothermal_C (double velocity[], double density[])
subroutine updateHydroAnalysisIsothermal_C (velocity, density) &
              BIND(C, NAME="updateHydroAnalysisIsothermal_C")
   real (wp), intent(in) :: velocity(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nDimensions, 0:grid%nFluids)
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)

   call updateHydroAnalysisPrimitive (grid, velocity, density)
   
   if(project_2D) then
      call updateHydroAnalysisPrimitive (grid_2D, velocity=SUM(velocity,DIM=2)/grid%nCells(2), &
         density=SUM(density,DIM=2)/grid%nCells(2))
   end if
   
end subroutine

! void updateHydroAnalysisMixture_C (double velocity[], double density[], double concentration[])
subroutine updateHydroAnalysisMixture_C (velocity, density, concentration) &
              BIND(C, NAME="updateHydroAnalysisMixture_C")
   real (wp), intent(in) :: velocity(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nDimensions, 0:grid%nFluids)
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in) :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)

   call updateHydroAnalysisPrimitive (grid, velocity=velocity, density=density, concentration=concentration)

   if(project_2D) then
      call updateHydroAnalysisPrimitive (grid_2D, velocity=SUM(velocity,DIM=2)/grid%nCells(2), &
         density=SUM(density,DIM=2)/grid%nCells(2), concentration=SUM(concentration,DIM=2)/grid%nCells(2))
   end if
   
end subroutine

! void updateHydroAnalysisStaggered_C (int nGhost, int nPassiveScalars, double vx[], double vy[], double vz[], int nGhostScalar, double scalar[])
! For simplicity, we alias density as an advected scalar here, if grid%nPassiveScalars>0, otherwise ignore it
! Note: You can pass NULL or a zero-sized array for vz or density if there is no data (it simply won't be touched)
subroutine updateHydroAnalysisStaggered_C (nGhost, nPassiveScalars, vx, vy, vz, nGhostScalar, scalar) &
              BIND(C, NAME="updateHydroAnalysisStaggered_C")
   integer(c_int), value :: nGhost, nPassiveScalars, nGhostScalar
   real (wp), intent(in), dimension(grid%nCells(1)+nGhost, grid%nCells(2), grid%nCells(3), 0:grid%nFluids) :: vx
   real (wp), intent(in), dimension(grid%nCells(1), grid%nCells(2)+nGhost, grid%nCells(3), 0:grid%nFluids) :: vy
   real (wp), intent(in), dimension(grid%nCells(1), grid%nCells(2), grid%nCells(3)+nGhost, 0:grid%nFluids) :: vz
   !real (wp), intent(in) :: scalar(grid%nCells(1), grid%nCells(2), grid%nCells(3), nPassiveScalars) ! Without ghost cells
   ! Allowing for ghost cells:
   real (wp), intent(in) :: scalar(1-nGhostScalar:grid%nCells(1)+nGhostScalar, &
                                   1-nGhostScalar:grid%nCells(2)+nGhostScalar, &
       1-merge(nGhostScalar,0,grid%nDimensions>2):grid%nCells(3)+merge(nGhostScalar,0,grid%nDimensions>2), &
                                   nPassiveScalars)
   
   !write(*,*) "HYDRO: Called with nPassiveScalars=", nPassiveScalars
      
   if(grid%nDimensions>2) then ! vz must be present if 3D      
      if(nPassiveScalars>0) then 
         call updateHydroAnalysisStaggered (grid, nGhost, vx, vy, vz, nGhostScalar=nGhostScalar, density=scalar)
      else
         call updateHydroAnalysisStaggered (grid, nGhost, vx, vy, vz, nGhostScalar=nGhostScalar)
      end if
   else
      if(nPassiveScalars>0) then
         call updateHydroAnalysisStaggered (grid, nGhost, vx, vy, nGhostScalar=nGhostScalar, density=scalar)
      else
         call updateHydroAnalysisStaggered (grid, nGhost, vx, vy, nGhostScalar=nGhostScalar)
      end if   
   end if      

   if(project_2D) then
      if(grid%nDimensions>2) then ! vz must be present if 3D      
         if(nPassiveScalars>0) then
            ! Observe that we average over all of vy here since there may be hard walls along the projected dimension
            call updateHydroAnalysisStaggered (grid_2D, nGhost=0, &
               vx=SUM(vx(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vy=SUM(vz(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vz=SUM(vy(1:grid%nCells(1), :, 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               nGhostScalar=0, density=&
                  SUM(scalar(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 1:nPassiveScalars),DIM=2)/grid%nCells(2))
         else
            call updateHydroAnalysisStaggered (grid_2D, nGhost=0, &
               vx=SUM(vx(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vy=SUM(vz(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vz=SUM(vy(1:grid%nCells(1), :, 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), nGhostScalar=0)
         end if
      else
         if(nPassiveScalars>0) then
            call updateHydroAnalysisStaggered (grid_2D, nGhost=0, &
               vx=SUM(vx(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vy=SUM(vy(1:grid%nCells(1), :, 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               nGhostScalar=0, density=&
                  SUM(scalar(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 1:nPassiveScalars),DIM=2)/grid%nCells(2))
         else
            call updateHydroAnalysisStaggered (grid_2D, nGhost=0, &
               vx=SUM(vx(1:grid%nCells(1), 1:grid%nCells(2), 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), &
               vy=SUM(vy(1:grid%nCells(1), :, 1:grid%nCells(3), 0:grid%nFluids), DIM=2)/grid%nCells(2), nGhostScalar=nGhostScalar)
         end if   
      end if      
   end if
   
end subroutine

! void writeToFiles_C(int id)
subroutine writeToFiles_C(id) BIND(C, NAME="writeToFiles_C")
   integer(c_int), value :: id ! An additional integer to append to file names, if positive 

   if(id>=0) then  
      call writeToFiles(grid, id)
   else
      call writeToFiles(grid)
   end if

   if(project_2D) then
      if(id>=0) then  
         call writeToFiles(grid_2D, id)
      else
         call writeToFiles(grid_2D)
      end if
   end if
       
end subroutine

!mcai-----begin-----------------------------

! void projectHydroGridMixture_C (double density[], double concentration[])
subroutine projectHydroGrid_C (density, concentration, filename) &
              BIND(C, NAME="projectHydroGrid_C")
   ! state the in/out args
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in) :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
   character(kind=c_char), target, dimension(*), intent(in) :: filename
   
   call projectHydroGridMixture (grid, density, concentration, filename)
   
end subroutine


! A. Donev: This routine is made to be callable from either Fortran or C codes:
subroutine projectHydroGridMixture (grid, density, concentration, filename)

   type(HydroGrid), target, intent(inout) :: grid ! A. Donev: This can be different from the module variable grid!
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in) :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
   character, dimension(*), intent(in) :: filename
   
   !local variables
   integer :: i, j, k
   ! Density temp, Density times Y coordinates, Concentratino times Y coordinates
   real (wp) :: density_tmp(grid%nCells(2)), concentration_tmp(grid%nCells(2)), Y_coor_tmp(grid%nCells(2))
   real (wp) :: density1_tmp(grid%nCells(2))

   ! size of DensityTimesY_coor and size of ConcentrTimesY_coor
   ! real (wp) :: DensityTimesY_coor(1, grid%nCells(1), grid%nCells(3)), ConcentrTimesY_coor(1, grid%nCells(1), grid%nCells(3))  
   real(wp), dimension(:,:), allocatable, target :: DensityTimesY_coor, Density1TimesY_coor, ConcentrTimesY_coor
   ! averaging values along Y direction
   real(wp), dimension(:,:), allocatable, target :: rho_avg_Y, rho1_avg_Y, c_avg_Y
   ! averaing values along XZ direction
   real(wp) :: rho_XZ_tmp, rho1_XZ_tmp, c_XZ_tmp
   real(wp), dimension(:), allocatable, target :: rho_avg_XZ, rho1_avg_XZ, c_avg_XZ
   
   !for writing vtk file
   character(len=30), target :: input_file_name     !local variable, file name for snapshots
   character(len=30), target :: hstat_file_name  ! file name for storing horizontal data as DAT file
   integer :: mesh_dims(3), dim, iVariance
   character(len=16), dimension(max(6,grid%nVariables)), target :: varnames
   

   input_file=""
   do i=1, len(input_file)
      if(filename(i)==C_NULL_CHAR) exit
      input_file(i:i)=filename(i)
   end do

  
   ! allocate memory for local variables
   allocate(DensityTimesY_coor(grid%nCells(1), grid%nCells(3)))
   allocate(Density1TimesY_coor(grid%nCells(1), grid%nCells(3)))
   allocate(ConcentrTimesY_coor(grid%nCells(1), grid%nCells(3)))

   allocate(rho_avg_Y(grid%nCells(1), grid%nCells(3)))
   allocate(rho1_avg_Y(grid%nCells(1), grid%nCells(3)))
   allocate(c_avg_Y(grid%nCells(1), grid%nCells(3)))
   
   allocate(rho_avg_XZ(grid%nCells(2)))
   allocate(rho1_avg_XZ(grid%nCells(2)))
   allocate(c_avg_XZ(grid%nCells(2)))


   ! sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i), it works for both 2D(grid%nCells(3)=1) and 3D
   do i=1, grid%nCells(1)
      do j=1, grid%nCells(3)
         do k=1, grid%nCells(2)
            density_tmp(k)=density(i, k, j, 0)
            ! A. Donev: We also want the center of mass using rho1:
            density1_tmp(k)=density(i, k, j, 0)*concentration(i, k, j, 1)
            concentration_tmp(k)=concentration(i, k, j, 1)
            ! Y_i cooridinates
            Y_coor_tmp(k)=(k-0.5_wp)*grid%systemLength(2)/grid%nCells(2)    ! it should be (k-1)*dy + dy/2+prob_lo (we assume prob_lo=0)
         enddo 
   
         DensityTimesY_coor(i, j)=DOT_PRODUCT(density_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         Density1TimesY_coor(i, j)=DOT_PRODUCT(density1_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         ConcentrTimesY_coor(i, j)=DOT_PRODUCT(concentration_tmp, Y_coor_tmp)/(sum(concentration_tmp, dim=1))
        
         ! Donev: Calculate also here the average rho, rho1 and c along the y direction
         rho_avg_Y(i, j)=sum(density_tmp, dim=1)/grid%nCells(2)
         rho1_avg_Y(i, j)=sum(density1_tmp, dim=1)/grid%nCells(2)
         c_avg_Y(i, j)=sum(concentration_tmp, dim=1)/grid%nCells(2)
         
      enddo
   enddo

   ! calcuate horizontal average for rho, rho1, and c
   rho_XZ_tmp=0.0_wp
   rho1_XZ_tmp=0.0_wp
   c_XZ_tmp=0.0_wp
   do k=1, grid%nCells(2)
      do i= 1, grid%nCells(1)
         do j=1, grid%nCells(3)
            rho_XZ_tmp=rho_XZ_tmp+density(i, k, j, 0)
            rho1_XZ_tmp=rho_XZ_tmp+density(i, k, j, 0)*concentration(i, k, j, 1)
            c_XZ_tmp=c_XZ_tmp+concentration(i, k, j, 1)
         end do 
      end do
      rho_avg_XZ(k)=rho_XZ_tmp
      rho1_avg_XZ(k)=rho1_XZ_tmp
      c_avg_XZ(k)=c_XZ_tmp
   end do  
     

   ! grid%nCells(3)=1 if it is 2D case 
   if(grid%nCells(3)>1) then

      !To write the data into VTK file, call WriteRectilinearVTKMesh
      !we have x z coorinates and the sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i) in stride

      input_file_name=trim(input_file)// ".CofM.vtk"
      write(*,*) "Writing average rho*Y and c*Y variables to file ", trim(input_file_name) 
      input_file_name=trim(input_file_name) // C_NULL_CHAR

      ! swap the axes
      mesh_dims(1)=grid%nCells(1)
      mesh_dims(2)=grid%nCells(3)
      mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)
   
      !~~~ we need XZ plane and plot the quantity~~~~~~~~~~~~~~~~
      varnames(1) = "rho_CofM" // C_NULL_CHAR
      varnames(2) = "c_CofM" // C_NULL_CHAR
      varnames(3) = "rho1_CofM" // C_NULL_CHAR
      varnames(4) = "rho_Avg" // C_NULL_CHAR
      varnames(5) = "rho1_Avg" // C_NULL_CHAR
      varnames(6) = "c_Avg" // C_NULL_CHAR
      
      call WriteRectilinearVTKMesh(filename=C_LOC(input_file_name(1:1)), &
         ub=0_c_int, dims=mesh_dims+1, &
         x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
         y=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(2)) /), &            
         z=(/ (0.0_wp, dim=0,1) /), &
         nvars=2, vardim=(/(1, dim=1, 2)/), centering=(/(0, dim=1, 2)/), &
         varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,2) /), &
         vars=(/ C_LOC(DensityTimesY_coor), C_LOC(Density1TimesY_coor), &
                C_LOC(ConcentrTimesY_coor), C_LOC(rho_avg_Y), C_LOC(rho1_avg_Y), C_LOC(c_avg_Y) /))

      ! write horizontal stat data into dat file
      hstat_file_name=trim(input_file) // ".avgXZ.dat"
      write(*,*) "Writing fluid density and concentraion to file ", trim(hstat_file_name) 
      open(2000, file=trim(hstat_file_name), status = "unknown", action = "write")
      hstat_file_name=trim(hstat_file_name)
      do k = 1, grid%nCells(2)            
         write(2000, '(1000(g17.9))') ((k-0.5_wp)*grid%systemLength(2)/grid%nCells(2)), &   ! Y coord at cell center
                      (rho_avg_XZ(k)), (rho1_avg_XZ(k)), (c_avg_XZ(k))    
      enddo  
      close(2000)

   else  ! if the grid is a 2D grid, we directly write down X coordinates and sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i)  
      
      ! write the 2D data into input_file directly
      input_file_name=trim(input_file) // ".2d.dat"
      open(1000, file=trim(input_file_name), status = "unknown", action = "write")
      do i = 1, grid%nCells(1)            
         write(1000, '(1000(g17.9))') ((i-0.5_wp)*grid%systemLength(1)/grid%nCells(1)), &           ! x coord at cell center
            (DensityTimesY_coor(i, 1)), (Density1TimesY_coor(i, 1)), (ConcentrTimesY_coor(i, 1)), &
            (rho_avg_Y(i, 1)), (rho1_avg_Y(i, 1)), (c_avg_Y(i, 1))    
      enddo  
      close(1000)
      
   end if

   !release memory
   deallocate(DensityTimesY_coor)
   deallocate(Density1TimesY_coor)
   deallocate(ConcentrTimesY_coor)

   deallocate(rho_avg_Y)
   deallocate(rho1_avg_Y)
   deallocate(c_avg_Y)
   
   deallocate(rho_avg_XZ)
   deallocate(rho1_avg_XZ)
   deallocate(c_avg_XZ)

end subroutine


! also write the 3D data file into a different vtk file 
! void writeHydroGridMixture_C (double density[], double concentration[])
subroutine writeHydroGridMixture_C (density, concentration, filename) &
              BIND(C, NAME="writeHydroGridMixture_C")
   ! state the in/out args
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in) :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
   character(kind=c_char), target, dimension(*), intent(in) :: filename
   
   call writeHydroGridMixture (grid, density, concentration, filename)
   
end subroutine

! A. Donev: Mingchao, make this routine accept density and concentration as arguments and only write them to the VTK file
subroutine writeHydroGridMixture(grid, density, concentration, filename) 

   type(HydroGrid), target, intent(inout) :: grid ! A. Donev: This can be different from the module variable grid!
   real (wp), intent(in), target :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in), target :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
   character, dimension(*), intent(in) :: filename
    
   !local variables
   integer :: mesh_dims(3), dim, iVariance
   character(len=16), dimension(max(4,grid%nVariables)), target :: varnames
   
   real(wp), dimension(:,:,:,:), allocatable, target :: velocity
   character(len=30), target :: input_file_name     !local variable, file name for snapshots
   integer :: i

   ! define the name of the statfile that will be written
  
   input_file=""
   do i=1, len(input_file)
      if(filename(i)==C_NULL_CHAR) exit
         input_file(i:i)=filename(i)
   end do 
    
   input_file_name=trim(input_file) //".DenCon.vtk"
   write(*,*) "Writing fluid density and concentraion to file ", trim(input_file_name) 
   input_file_name = trim(input_file_name) // C_NULL_CHAR

   
   mesh_dims=grid%nCells(1:3)
   if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)

   varnames(1) = "rho" // C_NULL_CHAR
   varnames(2) = "c" // C_NULL_CHAR
   varnames(3) = "c1" // C_NULL_CHAR
   varnames(4) = "c2" // C_NULL_CHAR
   CALL WriteRectilinearVTKMesh(filename=C_LOC(input_file_name(1:1)), &
        ub=0_c_int, dims=mesh_dims+1, &
        x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
        y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
        z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
        nvars=grid%nSpecies, vardim=(/(1, dim=1,grid%nSpecies)/), centering=(/(0, dim=1,grid%nSpecies)/), &
        varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,grid%nSpecies) /), &
        vars=(/ C_LOC(density), (C_LOC(concentration(:,:,:,dim)),dim=1,grid%nSpecies-1) /) )

end subroutine
!mcai-----end-----------------------------


END MODULE
