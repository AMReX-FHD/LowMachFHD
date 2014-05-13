subroutine main_driver()

  use convert_variables_module
  use fluid_model_module
  use matrix_utilities 
  use F95_LAPACK
  use probin_multispecies_module

  implicit none

  nspecies=3
  call test_chi(nspecies)

contains

subroutine test_chi(nspecies)
  implicit none

  integer, intent(in) :: nspecies
  real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda,chi,D_MS,Gama
  real(kind=dp_t), dimension(nspecies)          :: W,rho,drho,molarconc,molmass_in,molmass,chiw 
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
  rho(1)        = 0.0d0 ! 0.60
  rho(2)        = 0.0d0 ! 1.05d0 
  rho(3)        = 1.35d0
  molmass_in(1) = 1.0d0 
  molmass_in(2) = 2.0d0 
  molmass_in(3) = 3.0d0 
  Dbar_in(1)    = 1.0d0 
  Dbar_in(2)    = 0.5d0 
  Dbar_in(3)    = 1.5d0 
  inverse_type  = 1
  fraction_tolerance = 1e-14

  ! populate molar masses 
  molmass = molmass_in
  
  ! populate drho and add it to rho
  call correct_rho_with_drho_local(rho,drho)
  W = rho/sum(rho)

  ! Compute quantities consistently now
  call compute_molconc_rhotot_local(rho,rho_tot,molarconc,molmass,molmtot)  
  
  ! populate D_MS, Gama 
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

     print*, 'print chi' 
     print*, chi
     print*, 'print w' 
     print*, W
     print*, 'print chi*w' 
     print*, chiw
 
  enddo
 
  ! Undo the correction so we don't mess up conservation
  rho = rho - drho 
  
end subroutine test_chi

end subroutine main_driver
