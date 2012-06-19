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
   
   character(kind=c_char, len=1024) :: file_name="hydroGridOptions.nml"

CONTAINS   

! The C interface to this routine is
! void setHydroInputFile_C(char * filename)
subroutine setHydroInputFile_C(filename) BIND(C,NAME="setHydroInputFile_C")
   character(kind=c_char), dimension(*), intent(in) :: filename

   integer :: i
   
   file_name=""
   do i=1, len(file_name)
      if(filename(i)==C_NULL_CHAR) exit
      file_name(i:i)=filename(i)
   end do
   write(*,*) "Will read Fortran namelist from file ", trim(file_name)
   
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
   open (nameListFile, file = trim(file_name), status="old", action="read")
   
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

   type(HydroGrid), intent(inout) :: grid ! A. Donev: This can be different from the module variable grid!
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in) :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
   character, target, dimension(*), intent(in) :: filename

   !local variables
   integer :: i, j, k
   ! Density temp, Density times Y coordinates, Concentratino times Y coordinates
   real (wp) :: density_tmp(grid%nCells(2)), concentration_tmp(grid%nCells(2)), Y_coor_tmp(grid%nCells(2))
   real (wp) :: density1_tmp(grid%nCells(2))

   ! size of DensityTimesY_coor and size of ConcentrTimesY_coor
   ! real (wp) :: DensityTimesY_coor(1, grid%nCells(1), grid%nCells(3)), ConcentrTimesY_coor(1, grid%nCells(1), grid%nCells(3))  
   real(wp), dimension(:,:), allocatable, target :: DensityTimesY_coor, Density1TimesY_coor, ConcentrTimesY_coor
   ! averaging values along Y direction
   real(wp), dimension(:,:), allocatable, target :: rho_avg, rho1_avg, c_avg
  
   !for writing vtk file
   integer :: mesh_dims(3), dim, iVariance
   character(len=16), dimension(max(6,grid%nVariables)), target :: varnames
   

   ! write all data into the 3D VTK file
   !if(grid%writeMeansVTK) then
   !  call writeVelocityToVTK()
   !end if 

   ! allocate memory for local variables
   allocate(DensityTimesY_coor(grid%nCells(1), grid%nCells(3)))
   allocate(Density1TimesY_coor(grid%nCells(1), grid%nCells(3)))
   allocate(ConcentrTimesY_coor(grid%nCells(1), grid%nCells(3)))

   allocate(rho_avg(grid%nCells(1), grid%nCells(3)))
   allocate(rho1_avg(grid%nCells(1), grid%nCells(3)))
   allocate(c_avg(grid%nCells(1), grid%nCells(3)))


   ! sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i)
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
         ! Donev: Calculate also here the average rho, rho1 and c along the y direction

         DensityTimesY_coor(i, j)=DOT_PRODUCT(density_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         Density1TimesY_coor(i, j)=DOT_PRODUCT(density1_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         ConcentrTimesY_coor(i, j)=DOT_PRODUCT(concentration_tmp, Y_coor_tmp)/(sum(concentration_tmp, dim=1))
         
         rho_avg(i, j)=sum(density_tmp, dim=1)/grid%nCells(2)
         rho1_avg(i, j)=sum(density1_tmp, dim=1)/grid%nCells(2)
         c_avg(i, j)=sum(concentration_tmp, dim=1)/grid%nCells(2)
      enddo
   enddo

   ! grid%nCells(3)=0     ! For testing 2D case 
   if(grid%nCells(3)>1) then

   !To write the data into VTK file, call WriteRectilinearVTKMesh
   !we have x z coorinates and the sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i) in stride
 
      file_name=""
      do i=1, len(file_name)
         if(filename(i)==C_NULL_CHAR) exit
         file_name(i:i)=filename(i)
      end do
      write(*,*) "Writing average rho*Y and c*Y variables to file ", trim(file_name)
      
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
      
     
      call WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
         ub=0_c_int, dims=mesh_dims+1, &
         x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
         y=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(2)) /), &            
         z=(/ (0.0_wp, dim=0,1) /), &
         nvars=2, vardim=(/(1, dim=1, 2)/), centering=(/(0, dim=1, 2)/), &
         varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,2) /), &
         vars=(/ C_LOC(DensityTimesY_coor), C_LOC(Density1TimesY_coor), &
                C_LOC(ConcentrTimesY_coor), C_LOC(rho_avg), C_LOC(rho1_avg), C_LOC(c_avg) /))

   else  ! if the grid is a 2D grid, we directly write down X coordinates and sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i)  
      ! sum(rho1_i*Y_i)/sum(rho1_i) and sum(c_i*Y_i)/sum(c_i)
      do i=1, grid%nCells(1)
         do k=1, grid%nCells(2)
            density_tmp(k)=density(i, k, 1, 0)
            density1_tmp(k)=density(i, k, j, 0)*concentration(i, k, j, 1)
            concentration_tmp(k)=concentration(i, k, 1, 1)
            Y_coor_tmp(k)=(k-1/2)*grid%systemLength(2)/grid%nCells(2)    ! it should be (k-1)*dy + dy/2+prob_lo (we assume prob_lo=0)
         enddo 
         DensityTimesY_coor(i, 1)=DOT_PRODUCT(density_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         Density1TimesY_coor(i, j)=DOT_PRODUCT(density1_tmp, Y_coor_tmp)/(sum(density_tmp, dim=1))
         ConcentrTimesY_coor(i, 1)=DOT_PRODUCT(concentration_tmp, Y_coor_tmp)/(sum(concentration_tmp, dim=1))
         
         ! for 2D case, we actually just use 1D array
         rho_avg(i, 1)=sum(density_tmp, dim=1)/grid%nCells(2)
         rho1_avg(i, 1)=sum(density1_tmp, dim=1)/grid%nCells(2)
         c_avg(i, 1)=sum(concentration_tmp, dim=1)/grid%nCells(2)
      enddo

      file_name=""
      do i=1, len(file_name)
         if(filename(i)==C_NULL_CHAR) exit
         file_name(i:i)=filename(i)
      end do
      write(*,*) "Writing average rho*Y and c*Y variables to file ", trim(file_name)

      open(1000, file=trim(file_name), status = "unknown", action = "write")
      do i = 1, grid%nCells(1)            
         ! A. Donev: Do not use format * when writing data files
         write(11, '(1000(g17.9))') ((i-0.5_wp)*grid%systemLength(1)/grid%nCells(1)), &           ! x coord at cell center
            (DensityTimesY_coor(i, 1)), (Density1TimesY_coor(i, 1)), (ConcentrTimesY_coor(i, 1)), &
            (rho_avg(i, 1)), (rho1_avg(i, 1)), (c_avg(i, 1))    
      enddo  
      close(1000)
      
   end if

end subroutine



!subroutine writeSnapshotToVTK_C() &
!           BIND(C, NAME="writeSnapshotToVTK_C")
!      ! Write a snapshot of the instantaneous fields
!   call writeSnapshotToVTK()

!end subroutine

!subroutine writeSnapshotToVTK() 
   ! Write a snapshot of the instantaneous fields
!   integer :: mesh_dims(3), dim, iVariance
!   character(len=16), dimension(max(6,grid%nVariables)), target :: varnames

!   character(len=nMaxCharacters), target :: filename
      
!   real(wp), dimension(:,:,:,:), allocatable, target :: velocity

!   filename = trim(filenameBase) // ".snapshot.vtk"
!   write(*,*) "Writing instantaneous single-fluid variables to file ", trim(filename)
!   filename = trim(filename) // C_NULL_CHAR

!   mesh_dims=grid%nCells(1:3)
!   if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)

!   allocate(velocity(3, grid%nCells(1), grid%nCells(2), grid%nCells(3)))
!   do dim=1, grid%nDimensions
!      velocity(dim, :, :, :) = grid%primitive(:, :, :, grid%jIdx1 + dim - 1 , 0)
!   end do   
!   velocity(grid%nDimensions+1 : 3, :, :, :) = 0.0_wp

!   varnames(1) = "Density" // C_NULL_CHAR
!   varnames(2) = "Velocity" // C_NULL_CHAR
!   varnames(3) = "Scalar1" // C_NULL_CHAR
!   varnames(4) = "Scalar2" // C_NULL_CHAR
!   varnames(5) = "Scalar3" // C_NULL_CHAR
!   varnames(6:) = "Scalar" // C_NULL_CHAR
!   CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
!        ub=0_c_int, dims=mesh_dims+1, &
!        x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
!        y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
!        z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
!        nvars=grid%nVariables+1-grid%nDimensions, vardim=(/1, 3, (1, dim=3,grid%nVariables)/), &
!        centering=(/(0, dim=1,grid%nVariables)/), &
!        varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,grid%nVariables) /), &
!        vars=(/ C_LOC(grid%primitive(1,1,1,grid%mIdx,0)), C_LOC(velocity), &
!             (C_LOC(grid%primitive(1,1,1,grid%eIdx+dim,0)), dim=0,grid%nVariables-grid%eIdx) /))

!end subroutine

!mcai-----end-----------------------------


END MODULE
