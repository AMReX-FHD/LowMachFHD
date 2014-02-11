subroutine main_driver()

  use multifab_module
  use init_module
  use convert_variables_module
  use probin_multispecies_module
  use multifab_physbc_module
  use matrix_utilities 
  use F95_LAPACK

  implicit none

  call probin_multispecies_init()
  call test_chi(nspecies)

contains

subroutine test_chi(nspecies)
  implicit none

  integer, intent(in) :: nspecies
  real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda,chi,D_MS,Gama
  real(kind=dp_t), dimension(nspecies)          :: W,rho,molarconc,molmass,molmass_in,Dbar_in 
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
  rho(1)        = 0.6d0
  rho(2)        = 1.05d0
  rho(3)        = 1.35d0
  molmass_in(1) = 1.0d0 
  molmass_in(2) = 2.0d0 
  molmass_in(3) = 3.0d0 
  Dbar_in(1)    = 0.5d0 
  Dbar_in(2)    = 1.0d0 
  Dbar_in(3)    = 1.5d0 
 
  ! populate D_MS, Gama and molar masses 
  n=0; 
  do row=1, nspecies  
     do column=1, row-1
        n=n+1
        D_MS(row, column) = Dbar_in(n)
        D_MS(column, row) = D_MS(row, column) ! symmetric
        Gama(row, column) = 0.d0       
        Gama(column, row) = Gama(row, column) ! symmetric
     enddo
     D_MS(row, row) = 0.d0 ! self-diffusion is zero
     Gama(row, row) = 1.d0 ! set to unit matrix for time being
     molmass(row) = molmass_in(row)
  enddo
 
  ! compute rho_tot
  do n=1, nspecies  
     rho_tot = rho_tot + rho(n)
  enddo         
 
  ! calculate mass fraction,total molar mass (1/m=Sum(w_i/m_i)), molar
  ! concentration (x_i=m*w_i/m_i) 
  Sum_woverm=0.d0
  do n=1, nspecies  
     W(n) = rho(n)/rho_tot
     Sum_woverm = Sum_woverm + W(n)/molmass(n)
  enddo
  molmtot = 1.0d0/Sum_woverm 
  do n=1, nspecies 
     molarconc(n) = molmtot*W(n)/molmass(n)
  enddo

  ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
  ! expressed in terms of molmtot,mi,rhotot etc. 
  do row=1, nspecies  
     do column=1, row-1
        Lambda(row, column) = -molarconc(row)*molarconc(column)/D_MS(row,column)
        Lambda(column, row) = Lambda(row, column)
     enddo
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

  do loop=1,2
  
     ! compute chi either selecting inverse/pseudoinverse or iterative methods 
     if(loop==1) then
        call compute_chi_lapack(Lambda,chi,W)
        print*, 'compute chi via inverse/p-inverse'
     else
!       call Dbar2chi_iterative(nspecies,10,D_MS,W,molarconc,chi)
!       call Dbar2chi_iterative(nspecies,5,D_MS,molmass,molarconc,chi)
        call Dbar2chi_iterative(nspecies,3,D_MS,molmass,molarconc,chi)
        print*, 'compute chi via iterative methods'
     endif

     !if(.false.) then
     do row=1, nspecies
        do column=1, nspecies
           print*, chi(row, column)
        enddo
        print*, ''
     enddo
     !endif 
 
  end do
  
end subroutine test_chi

end subroutine main_driver
