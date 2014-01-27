program test_chi

  use multifab_module
  use init_module
  use populate_DbarGama_module
  use convert_variables_module
  use probin_common_module
  use probin_multispecies_module
  use multifab_physbc_module
  use matrix_utilities 
  use F95_LAPACK

  implicit none

  real(kind=dp_t), allocatable :: Dbar(:,:)
  real(kind=dp_t), allocatable :: Gama(:,:)
  real(kind=dp_t), allocatable :: mass(:) 

  call probin_common_init()
  call probin_multispecies_init()
  
  call test_chi_routine(nspecies)

contains


subroutine test_chi_routine(nspecies)
  implicit none

  integer, intent(in) :: nspecies

  real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda,chi
  real(kind=dp_t), dimension(nspecies)          :: W,rho,molarconc 
  real(kind=dp_t)                               :: rho_tot,molmtot,Sum_woverm,tolerance
  integer                                       :: i,j,k,n,row,column 

  allocate(Dbar(nspecies,nspecies))
  allocate(Gama(nspecies,nspecies))
  allocate(mass(nspecies))

  ! populate Dbar, Gama and molar masses 
  call populate_DbarGama(Dbar,Gama,mass) 

  Lambda    = 0.d0         
  chi       = 0.d0         
  W         = 0.d0
  rho       = 0.d0
  rho_tot   = 0.d0
  molarconc = 0.d0
  tolerance = 1e-13
  
  ! populate rho 
  rho(1)    = 0.6d0
  rho(2)    = 1.05d0
  rho(3)    = 1.35d0

  ! compute rho_tot
  do n=1, nspecies  
     rho_tot = rho_tot + rho(n)
  enddo         
 
  ! calculate mass fraction,total molar mass (1/m=Sum(w_i/m_i)), molar
  ! concentration (x_i=m*w_i/m_i) 
  Sum_woverm=0.d0
  do n=1, nspecies  
     W(n) = rho(n)/rho_tot
     Sum_woverm = Sum_woverm + W(n)/mass(n)
  enddo
  molmtot = 1.0d0/Sum_woverm 
  do n=1, nspecies 
     molarconc(n) = molmtot*W(n)/mass(n)
  enddo

  ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
  ! expressed in terms of molmtot,mi,rhotot etc. 
  do row=1, nspecies  
     do column=1, row-1
        Lambda(row, column) = -molarconc(i,j,row)*molarconc(i,j,column)/Dbar(row,column)
        Lambda(column, row) = Lambda(row, column) 
     enddo
     W(row) = rho(i,j,row)/rho_tot(i,j)
  enddo

  ! compute Lambda_ii
  do row=1,nspecies
     Sum_knoti = 0.d0
     do column=1,nspecies
        if(column.ne.row) then
           Sum_knoti = Sum_knoti - Lambda(row,column)
        endif
        Lambda(row,row) = Sum_knoti
     enddo
  enddo

  ! compute chi either selecting inverse/pseudoinverse or iterative methods 
  if(use_lapack) then
     call populate_coefficient(Lambda(:,:),chi(:,:),Gama(:,:),W(:),tolerance)
  else
     call Dbar2chi_iterative(nspecies,3,Dbar(:,:),W(:),molarconc(i,j,:),chi(:,:))
  endif

  ! print to match with previous code
  if(use_lapack) then 
     print*, 'print chi via inverse/p-inverse'
  else 
     print*, 'print chi via iterative methods'
  endif 
     do row=1, nspecies
        do column=1, nspecies
           print*, chi(row, column)
        enddo
        print*, ''
     enddo

  deallocate(Dbar)
  deallocate(Gama)
  deallocate(W)

end subroutine

end program

