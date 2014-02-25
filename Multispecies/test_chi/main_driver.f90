subroutine main_driver()

  use convert_variables_module
  use fluid_model_module
  use matrix_utilities 
  use F95_LAPACK
  use probin_multispecies_module

  implicit none

  nspecies=2
  call test_chi(2)

contains

subroutine test_chi(nspecies)
  implicit none

  integer, intent(in) :: nspecies
  real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda,chi,D_MS,Gama
  real(kind=dp_t), dimension(nspecies)          :: W,rho,drho,molarconc,molmass,molmass_in,Dbar_in,chiw 
  real(kind=dp_t)                               :: rho_tot,molmtot,Sum_woverm,Sum_knoti
  integer                                       :: i,j,k,n,row,column,loop

  ! free up memory 
  D_MS       = 0.d0         
  Lambda     = 0.d0         
  chi        = 0.d0         
  rho_tot    = 0.d0
  molarconc  = 0.d0
  Sum_knoti  = 0.d0
  Sum_woverm = 0.d0

  ! initialize conserved and constant quantities
  rho(1)        = 0 
  rho(2)        = 1 
  !rho(3)        = 1.35d0
  molmass_in(1) = 1.0d0 
  molmass_in(2) = 2.0d0 
  !molmass_in(3) = 3.0d0 
  Dbar_in(1)    = 1.0d0 
  !Dbar_in(2)    = 0.5d0 
  !Dbar_in(3)    = 1.5d0 
  fraction_tolerance = 1e-6

  ! change 0 with tolerance to prevent division by zero in case species
  ! density, molar concentration or total density = 0.
  rho_tot = sum(rho)
  write(*,*) "TEST=", rho_tot, fraction_tolerance
  do row=1, nspecies
        if(rho(row) .lt. fraction_tolerance*rho_tot) then
           drho(row) = fraction_tolerance*rho_tot
       else
           drho(row) = 0.0d0
       endif
       write(*,*) row, rho(row), fraction_tolerance*rho_tot, drho(row)
 enddo
  rho = rho + drho ! Add a correction to make sure no mass or mole fraction is zero
  write(*,*) "new rho=", rho, " new w=", rho/sum(rho)
  
  molmass = molmass_in
  
  ! Compute quantities consistently now
  call compute_molconc_rhotot_local(rho,rho_tot,molarconc,molmass,molmtot)  
  write(*,*) "rho=",rho_tot, " x=", molarconc, " m=", molmtot
  
  ! populate D_MS, Gama and molar masses 
  call compute_D_MSGama_local(rho,rho_tot,molarconc,molmtot,D_MS,Gama)

  do loop=1,2
  
     ! compute chi either selecting inverse/pseudoinverse or iterative methods 
     if(loop==1) then
        print*, 'compute chi via inverse/p-inverse'
        use_lapack = .true.
     else
        print*, 'compute chi via iterative methods'
        use_lapack = .false.
     endif
     call compute_chi_local(rho,rho_tot,molarconc,molmass,chi,D_MS)

     chiw = matmul(chi,W)
     do row=1, nspecies
        do column=1, nspecies
           print*, chi(row, column)
        enddo
        print*, ''
     enddo
 
     print*, 'print chi*w' 
     do row=1, nspecies
        print*, chiw(row)
     enddo
     print*, ''
 
  end do
  
  rho = rho - drho ! UNDO the correction so we don't mess up conservation
  
end subroutine test_chi

end subroutine main_driver
