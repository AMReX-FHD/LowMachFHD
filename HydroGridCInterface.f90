MODULE HydroGridCInterface ! Interface to my HydroGrid module
   USE, INTRINSIC :: ISO_C_BINDING
   USE HydroGridModule
   IMPLICIT NONE
   PUBLIC

   INTEGER, PARAMETER :: wp = c_double
   
   type(HydroGrid), save :: grid, grid_2D
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

END MODULE
