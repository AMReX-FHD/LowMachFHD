module eos_model_module

  use bl_types
  use bl_IO_module
  use bl_constants_module
  use probin_common_module, only: k_B, Runiv, molmass
  use probin_multispecies_module, only: nspecies
  use probin_energy_module, only: use_fake_diff, fake_diff_coeff, fake_soret_factor, &
                                  heating_type, p0_in, dia_in, int_deg_free_in

  implicit none

  integer, parameter  :: MAXSPECIES = 20

  REAL*8 , SAVE :: AVOGADRO

  ! molmass from namelist is in g/molecule.  molecular_weight is in g/mole
  REAL*8, SAVE, allocatable :: molecular_weight(:) 
  REAL*8, SAVE, allocatable :: dia(:)
  REAL*8, SAVE, allocatable :: int_deg_free(:)
  REAL*8, save, allocatable :: cvgas(:)
  REAL*8, save, allocatable :: cpgas(:)
  REAL*8, save, allocatable :: e0ref(:)

!  matrix components for computation of diffusion that are independent of state
!  these are built duruing initialization and saved for efficiency
  REAL*8, save, allocatable :: Dbinbar(:,:),omega11bar(:,:),sigma11bar(:,:),diamat(:,:)
  REAL*8, save, allocatable :: amat1bar(:,:),amat2bar(:,:),alphabar(:,:)

contains

  subroutine eos_model_init()

    integer :: iwrk, nfit, i, ic, ii, j
    integer :: dochem, dostrang
    double precision :: rwrk
    integer ns
    real*8 :: mu,Fij,Fijstar,fact1
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs

    allocate(dia(nspecies))
    allocate(molecular_weight(nspecies))
    allocate(int_deg_free(nspecies))
    allocate(cvgas(nspecies))
    allocate(cpgas(nspecies))
    allocate(e0ref(nspecies))

    AVOGADRO = Runiv / k_B

    do ns=1,nspecies

       dia(ns) = dia_in(ns)

       ! molmass from namelist is in g/molecule.  Converting to g/mole
       molecular_weight(ns) = molmass(ns)*AVOGADRO

       int_deg_free(ns) = int_deg_free_in(ns)

    enddo

    do ns =  1,nspecies
       e0ref(ns) = 0.d0
       cvgas(ns) = 0.5d0*(3+int_deg_free(ns))*Runiv/molecular_weight(ns)
       cpgas(ns) = 0.5d0*(5+int_deg_free(ns))*Runiv/molecular_weight(ns)
    enddo

    write(6,*)"e0fref",e0ref
    write(6,*)"dia",dia
    write(6,*)"molecular_weight",molecular_weight
    write(6,*)"molecular_mass",molmass(1:nspecies)
    write(6,*)"p0_in",p0_in
    write(6,*)"heating_type",heating_type
    write(6,*)"int_deg_free",int_deg_free
    write(6,*)"cvgas",cvgas
    write(6,*)"cpgas",cpgas

    !    build matrices for computing diffusion
    allocate(diamat(1:nspecies,1:nspecies))
    allocate(Dbinbar(1:nspecies,1:nspecies))
    allocate(omega11bar(1:nspecies,1:nspecies))
    allocate(sigma11bar(1:nspecies,1:nspecies))
    allocate(amat1bar(1:nspecies,1:nspecies))
    allocate(amat2bar(1:nspecies,1:nspecies))
    allocate(alphabar(1:nspecies,1:nspecies))


    do i = 1, nspecies
       do j = 1, nspecies

          diamat(i,j) = 0.5d0*(dia(i) + dia(j))

          Dbinbar(i,j) = 3.0d0/16.0d0*sqrt(2.0d0*M_PI*k_B**3.d00    &   
               *(molmass(i)+molmass(j))/molmass(i)/molmass(j))/(M_PI*diamat(i,j)**2.0d0)

          mu = molmass(i)*molmass(j)/(molmass(i) + molmass(j))
          omega11bar(i,j) = sqrt(M_PI*k_B/(2.0d0*mu))*diamat(i,j)**2.0d0

          sigma11bar(i,j) = sqrt( k_B/(2.0d0*M_PI*molmass(i)*molmass(j)/  &
               (molmass(i)+molmass(j))) )*M_PI*(diamat(i,j)**2.0d0)

          Fij = (6.0d0*molmass(i)*molmass(i) + 13.0d0/5.0d0*molmass(j)*molmass(j) +   &
               16.0d0/5.0d0*molmass(i)*molmass(j))/((molmass(i)+molmass(j))**2.0d0)
          Fijstar = -27.0d0/5.0d0

          amat1bar(i,j) = 5.0d0/(k_B)*molmass(i)*molmass(j)/    &
               (molmass(i)+molmass(j))*Fij*sigma11bar(i,j)
          amat2bar(i,j) = 5.0d0/(k_B)*molmass(i)*molmass(j)*molmass(i)*molmass(j)/    &
               ((molmass(i)+molmass(j))**3.0d0)*Fijstar*sigma11bar(i,j)
          fact1 = molmass(i)*molmass(j)/(molmass(i)+molmass(j))

          alphabar(i,j) = 8.0d0/(3.0d0*k_B)*fact1*fact1*   &
               (-.5d0*sigma11bar(i,j))

       enddo
    enddo

  end subroutine eos_model_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes $c_p(\wb,T)$.  In this EOS this is not dependent on temperature
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKCPBS(temp,Yk,IWRK,RWRK,cpmix)
    
    real(kind=8) :: temp,Yk(1:nspecies),RWRK,cpmix
    integer :: IWRK, ns

    cpmix = 0.0d0
    do ns = 1, nspecies
       cpmix = cpmix + Yk(ns)*cpgas(ns)
    enddo
 
    return
  end subroutine CKCPBS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes $c_v(\wb,T)$.  In this EOS this is not dependent on temperature.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKCVBS(temp,Yk,IWRK,RWRK,cvmix)

    real(kind=8) :: temp,Yk(1:nspecies),RWRK,cvmix
    integer :: IWRK, ns
    
    cvmix = 0.0d0
    do ns = 1, nspecies
       cvmix = cvmix + Yk(ns)*cvgas(ns)
    enddo
 
    return
  end subroutine CKCVBS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes $P(\rho,\wb,T)$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKPY(rho, temp, Yk, iwrk, rwrk, pt) 

    real(kind=8) :: rho, temp, Yk(1:nspecies), rwrk, pt
    integer :: iwrk, ns
    real(kind=8) :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Yk(ns)/molecular_weight(ns)
    enddo
    molmix = 1.0d0/molmix

    pt = rho*(Runiv/molmix)*temp  

    if(temp.ge.6000.0d0) then
       print*, 'BUG IN CKPY ', rho, temp,Yk, molmix  
       stop
    endif
 
    return
  end subroutine CKPY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  computes $h_m(T)$.  Note $h = \sum_k w_k h_k$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ckhms(T, iwrk, rwrk, hk) 
    
    real(kind=8) :: T, rwrk, hk(1:nspecies)
    integer :: iwrk, ns   
    
    do ns = 1, nspecies 
       hk(ns) = e0ref(ns) + cpgas(ns)*T
    enddo

    return
  end subroutine ckhms

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes mole fractions, $\xb$, from mass fractions, $\wb$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKYTX(Yk, iwrk, rwrk, Xk)

    real(kind=8) :: Yk(1:nspecies), Xk(1:nspecies), rwrk
    integer :: iwrk, ns
    real*8 :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Yk(ns)/molecular_weight(ns)
    enddo
    do ns = 1, nspecies
       Xk(ns) = Yk(ns)/(molmix*molecular_weight(ns))
    enddo

    return
  end subroutine CKYTX

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! mass fractions, $\wb$, from mole fractions, $\xb$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKXTY(Xk, iwrk, rwrk, Yk)
    
    real(kind=8) :: Yk(1:nspecies), Xk(1:nspecies), rwrk
    integer :: iwrk, ns
    real*8 :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Xk(ns)*molecular_weight(ns)
    enddo
    do ns = 1, nspecies
       Yk(ns) = Xk(ns)*molecular_weight(ns)/molmix
    enddo

    return
  end subroutine CKXTY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes the mixture-averaged molecular mass, $\bar{m}$, given $\wb$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKMMWY(Yk, iwrk, rwrk, molmix)

    real(kind=8) :: Yk(1:nspecies), rwrk, molmix
    integer :: iwrk, ns

    molmix = 0.0d0
    do ns = 1, nspecies
       ! CHANGING UNITS TO BE COMPATIBLE WITH compute_S
!       molmix = molmix + Yk(ns)/molecular_weight(ns)
       molmix = molmix + Yk(ns)/molmass(ns)
    enddo
    molmix = 1.0d0/molmix

    return
  end subroutine CKMMWY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! computes $\rho(P,\wb,T)$
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CKRHOY(pt,temp,Yk,IWRK,RWRK,rho)

    real(kind=8) :: pt, temp, Yk(1:nspecies), rwrk, rho
    integer :: iwrk, ns     
    real(kind=8) :: molmix

    molmix = 0.0d0
    do ns = 1, nspecies
       molmix = molmix + Yk(ns)/molecular_weight(ns)
    enddo
    molmix = 1.0d0/molmix

    rho = pt/(Runiv/molmix)/temp

    return
  end subroutine CKRHOY

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! takes $\rho,T,P,\wb$, and $\xb$ as inputs and computes the following:
  ! eta     is the dynamic viscosity, $\eta$
  ! kappa   is the thermal conductivity, $\lambda$
  ! zeta    is the bulk viscosity, $\kappa$
  ! diff_ij is the diffusion matrix, $\chi$
  ! chitil are the thermodiffusion coefficients, $\zeta$.
  !                (I modified the routine to make this so)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ideal_mixture_transport(density,temperature,pressure,Yk,Xk,eta,kappa, &
                                     zeta,diff_ij,chitil,nspec)
    
    real*8, intent(in) :: density,temperature,pressure
    integer, intent(in) :: nspec
    real*8, intent(in), dimension(1:nspec) :: Yk,Xk
    real*8, intent(out) :: eta,kappa,zeta

    real*8, intent(out), dimension(1:nspec,1:nspec) :: diff_ij
    real*8, dimension(1:nspec,1:nspec) :: old_diff_ij
    real*8, dimension(1:nspec) :: gam_matrix
    real*8, intent(out), dimension(1:nspec) :: chitil
    real(kind=8) :: ptc
    real(kind=8), dimension(1:nspec) :: tdv
    real(kind=8) :: Dbin(1:nspec,1:nspec)
    real(kind=8) :: omega11(1:nspec,1:nspec)
    real(kind=8) :: sigma11(1:nspec,1:nspec)
    real(kind=8) :: a_ij1(1:nspec,1:nspec)
    real(kind=8) :: a_ij2(1:nspec,1:nspec)


    !     Local variables                                                   
    real*8 muxp,muyp 
    real*8 lamxp,lamyp 
    real*8 kxp,kyp 
    real*8 phiflx 
    real*8 et0,et1,etx,zeta0,zeta1,zetax
    real*8 etax,etay,etaz,xzeta,yzeta,zzeta
    real*8 kap0,kap1,kapx,xkap,ykap,zkap 
    real*8 cu0, cu1,cuy,cuz, conc,temp 
    real*8 diamx,mx,omc,ctemp,Ix,Kx,K0,K1,Kb
    real*8 s0, s1
    real*8 prom,ratm,difm,summ,omcsq,comc,csq 


    real*8 :: Mwmix, wbar, sqrtT
    integer :: ii, jj, ns ,i ,j
    real(kind=8), dimension(1:nspec) :: xx
    real(kind=8), dimension(1:nspec,1:nspec) :: dd
    real(kind=8), dimension(1:nspec) :: Cv

    real(kind=8), dimension(1:nspec) :: xxtr, yytr
    real(kind=8), dimension(1:nspec) :: cpk  

    real(kind=8) :: rwrk
    integer :: iwrk
    integer :: old


    !====================================================

    !        do ii = 1, nspec
    !         do jj = 1, nspec
    !           dd(ii,jj) = 0.5d0*(dia(ii) + dia(jj))
    !         enddo
    !        enddo


    if ( use_fake_diff .eq. 1)then

       eta = 0.d0
       zeta = 0.d0
       kappa = 0.d0
       chitil = 0.d0

       diff_ij = 0.d0

       do ii = 1,nspec
          diff_ij(ii,ii) = diff_ij(ii,ii)+fake_diff_coeff
       enddo

       !         do i = 1, nspec
       !          do j = 1, nspec
       !            diff_ij(i,j) = diff_ij(i,j)*Yk(i)
       !          enddo
       !         enddo   


    else

       ! mole fractions correction - EGLIB
       do ii = 1, nspec
          xxtr(ii) = Xk(ii) + (1.0d-16)*(sum(Xk(:))/dble(nspec)-Xk(ii))
       enddo

       ! molecular weight of mixture - EGLIB
       Mwmix = 0.0d0
       do ii = 1, nspec
          MWmix = MWmix + xxtr(ii)*molecular_weight(ii)
       enddo
       wbar = MWmix

       ! mass fractions correction - EGLIB
       do ii = 1, nspec
          yytr(ii) = molecular_weight(ii)/MWmix*xxtr(ii)
       enddo

       !      old = 0

       !      if(old .eq.1)then

       !      call viscosityM(xxtr,temperature,molmass,k_B,dd,eta,nspec) 

       !   not bulk viscosity
       !      zeta = 0.d0

       !      call lambdaM(xxtr,temperature,molmass,k_B,dd,Cv,kappa,nspec)



       !            if(eta.lt.0.0d0) then
       !             print*, 'bug in eta: ', eta
       !             stop
       !            endif
       !            if(kappa.lt.0.0d0) then
       !             print*, 'bug in kappa: ', kappa, eta, Yk(:)
       !             print*, 'dd: ', dd(1,1), dd(2,2), dd(3,3)
       !             print*, 'Cv: ', Cv, temperature
       !             stop
       !            endif


       !      call diffusivityM(xxtr,yytr,temperature,molmass,k_B,dd,pressure,diff_ij,nspec)
       !      call thermalDiffM(xxtr,temperature,molmass,k_B,dd,chitil,nspec)

       !      else

       ! find binary diffusion coefficients  
       ! HCB 8.2-9   
       sqrtT = dsqrt(temperature)
       do i = 1, nspec
          do j = 1, nspec

             !            Dbin(i,j) = 3.0d0/16.0d0*sqrt(2.0d0*M_PI*k_B**3.d00*temperature**3.0d0    &
             !              *(molmass(i)+molmass(j))/molmass(i)/molmass(j))  &
             !              /(pressure*M_PI*dd(i,j)**2.0d0)

             Dbin(i,j) = Dbinbar(i,j)*temperature*sqrtT/pressure
             omega11(i,j) = omega11bar(i,j)*sqrtT
             sigma11(i,j) = sigma11bar(i,j)*sqrtT
             a_ij1(i,j) = amat1bar(i,j)/sqrtT
             a_ij2(i,j) = amat2bar(i,j)/sqrtT

          enddo
       enddo

       call visc_lin(nspec,k_B,omega11,yytr,temperature,density,molmass,eta)
       zeta = 0.d0

       call lambda_lin(nspec,k_B,Dbin,omega11,yytr,temperature,density,molmass,kappa)

       !  jmax = 3
       call D_GIOVANGIGLI(nspec,Dbin,yytr,xxtr,diff_ij)

       !  jmax = 0
       !        call PGAMMAP(nspec,yytr,xxtr,Dbin,diff_ij,molecular_weight,wbar)

       call thermalDiff(nspec,sigma11,a_ij1,a_ij2,alphabar,xxtr,sqrtT,molmass,chitil)
       chitil = chitil*fake_soret_factor

       ! multiply chitil by mole fractions so they become thermodiffusion
       ! coefficients, zeta
       do i=1,nspec
          chitil(i) = chitil(i)*Xk(i)
       end do

       !      endif

       !    print *," in ideal after",pressure,temperature,density
       !    print *," xk,yk ", Xk, Yk

       !      print *," eta, kappa ",eta,kappa
       !      print *, "Diff ", diff_ij
       !      print *, " chitil ",chitil
       !     print *,"Dbin ", Dbin
       !      stop
       !      stop

    endif

    return

  end subroutine ideal_mixture_transport

end module eos_model_module
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are local routines needed above
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------

subroutine D_GIOVANGIGLI(nspec,Dbin,Ykp,Xkp,D_tilde)

  implicit none 

  integer :: nspec, i, j, k, jj
  real(kind=8) :: Ykp(1:nspec), Xkp(1:nspec)
  real(kind=8) :: Dbin(1:nspec,1:nspec), term1, term2
  real(kind=8) :: D_tilde(1:nspec,1:nspec), Di(1:nspec), Diff_ij(1:nspec,1:nspec)
  real(kind=8) :: Deltamat(1:nspec,1:nspec), Zmat(1:nspec,1:nspec)
  real(kind=8), dimension(1:nspec,1:nspec) :: Pmat, Jmat
  real(kind=8), dimension(1:nspec) :: Minv, Mmat
  real(kind=8), dimension(1:nspec,1:nspec) :: PJ, matrix1, matrix2
  real(kind=8) :: scr

  integer :: jmax


  jmax = 3 

  ! Find Di matrix 
  do i = 1, nspec
     term2 = 0.0d0
     do j = 1, nspec
        if(j.ne.i) then
           term2 = term2 + Xkp(j)/Dbin(i,j)
        endif
     enddo
     Di(i) = (1.d0-Ykp(i))/term2 
  enddo


  ! Compute Mmat and Minv
  do i = 1, nspec
     Mmat(i) = Xkp(i)/Di(i)
     Minv(i) = Di(i)/Xkp(i)
  enddo


  ! Compute P matrix
  Pmat = 0.0d0
  do i = 1, nspec
     do j = 1, nspec
        Pmat(i,j) = - Ykp(j) 
        if(i.eq.j) then
           Pmat(i,j) =  Pmat(i,j) + 1.0d0  
        endif
     enddo
  enddo


  ! Compute Deltamat
  Deltamat = 0.0d0 
  do i = 1, nspec
     do j = 1, nspec
        if(i.eq.j) then
           term1 = 0.0d0
           do k = 1, nspec
              if(k.ne.i) then
                 term1 = term1 + Xkp(i)*Xkp(k)/Dbin(i,k)
              endif
           enddo
           Deltamat(i,i) = term1
        else
           Deltamat(i,j) = -Xkp(i)*Xkp(j)/Dbin(i,j) 
        endif
        Zmat(i,j) = -Deltamat(i,j)
     enddo
  enddo


  ! Compute Zmat
  do i = 1, nspec
     Zmat(i,i) = Zmat(i,i) + Mmat(i)
  enddo


  ! Compute Jmat
  do i = 1, nspec
     do j = 1, nspec
        Jmat(i,j) = Minv(i)*Zmat(i,j)
     enddo
  enddo

  ! Compute PJ
  PJ = 0.0d0
  do i = 1, nspec
     do j = 1, nspec
        do k = 1, nspec
           PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
        enddo
     enddo
  enddo



  ! Compute P M^-1 Pt; store it in matrix2
  do i = 1, nspec
     do j = 1, nspec
        scr = 0.d0
        do k = 1, nspec
           scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
           ! notice the change in indices for Pmat to represent Pmat^t
        enddo
        matrix2(i,j) = scr
        Diff_ij(i,j) = scr
     enddo
  enddo

  if(jmax.gt.0)then

     do jj = 1,jmax


        !         matrix1=0
        do i = 1, nspec
           do j = 1, nspec
              scr = 0.d0
              do k = 1, nspec
                 scr = scr + PJ(i,k)*Diff_ij(k,j)
              enddo
              matrix1(i,j) = scr+matrix2(i,j)
           enddo
        enddo

        Diff_ij=matrix1

     enddo

  endif



  ! Compute D_tilde
  do i = 1, nspec
     do j = 1, nspec
        D_tilde(i,j) = Diff_ij(i,j)*Ykp(i)
     enddo
  enddo


  return
end subroutine D_GIOVANGIGLI


!-------------------------------------------------

subroutine visc_lin(nspec,k_B,omega11,Ykp,T,rho,mk,etaMix)

  implicit none 

  integer :: nspec
  real(kind=8), dimension(1:nspec) ::  mk, Ykp,  nk
  real(kind=8), dimension(1:nspec,1:nspec) :: d

  real(kind=8) :: T, rho, etaMix, sum1
  integer :: ii, jj, kk
  real*8 :: k_B
  real(kind=8), dimension(1:nspec,1:nspec) :: omega11, QoR
  real(kind=8) :: rhs(1:nspec), bSonine(1:nspec), diag(1:nspec)
  integer :: ip(1:nspec)


  ! find number density
  do ii = 1, nspec
     nk(ii) = rho*Ykp(ii)/mk(ii)
  enddo


  do  ii = 1, nspec
     diag(ii) = 0.d0
     do kk = 1, nspec
        diag(ii) = diag(ii) + nk(kk)*mk(kk)/(mk(ii)+mk(kk))**2.0d0 *      &
             (5.0d0*mk(ii)*omega11(ii,kk) +  &
             3.0d0*mk(kk)*omega11(ii,kk))

     enddo
  enddo

  do  ii = 1, nspec
     do jj = 1, nspec
        sum1 = -nk(jj)*mk(jj)/(mk(ii)+mk(jj))**2.0d0 *      &
             2.0d0*mk(jj)*omega11(ii,jj) 

        if(ii.eq.jj)then
           sum1 = sum1 + diag(ii)
        endif

        QoR(ii,jj) = -(16.0d0/15.0d0)*(mk(ii)/mk(jj))*sum1  
        ! HCB 7.4-62,63
     enddo
     bsonine(ii) = -1.0d0
  enddo



  call decomp(nspec,nspec,QoR,ip)
  call solve(nspec,nspec,QoR,bsonine,ip)


  sum1 = 0.0d0
  do ii = 1, nspec
     sum1 = sum1 + nk(ii)*bSonine(ii)
     ! HCB 7.4-56
  enddo


  etaMix = 0.5d0*k_B*T*sum1




  return
end subroutine visc_lin

!-------------------------------------------------------------------------

subroutine lambda_lin(nspec,k_B,Dbin,omega11,Ykp,T,rho,mk,lammix)

  implicit none 

  integer :: nspec   
  real(kind=8) :: k_B
  real(kind=8) :: mw, T, rho, mk(1:nspec),  lammix
  real(kind=8), dimension(1:nspec) :: Ykp,  nk ,beta
  integer :: i, j, k
  real(kind=8) :: sum00, sum01, sum11, mu
  real(kind=8), dimension(1:nspec,1:nspec) :: omega11
  real(kind=8), dimension(1:2*nspec,1:2*nspec) :: QQ
  real(kind=8) ::  aSonine(1:2*nspec), D_T(1:nspec)
  real(kind=8) :: Dbin(1:nspec,1:nspec), sum1, lamdaprime, ntotal
  real(kind=8) :: diag1(1:nspec),diag2(1:nspec),diag3(1:nspec),diag4(1:nspec)
  real(kind=8) :: ratm,sqratm,scr
  integer :: ip(2*nspec)




  ! find number density
  do i = 1, nspec
     nk(i) = rho*Ykp(i)/mk(i)
  enddo


  diag1=0.d0
  diag2=0.d0
  diag3=0.d0
  diag4=0.d0

  do i=1,nspec
     do k=1,nspec

        diag1(i) =diag1(i)+ nk(k)*mk(k)/(mk(i)+mk(k)) *   &
             nk(i)*mk(i)*omega11(i,k)
        if(k.ne.i)then
           diag2(i) = diag2(i)+nk(k)*mk(k)/(mk(i)+mk(k))*omega11(i,k)
        endif
        diag3(i) = diag3(i) + nk(i)*nk(k)*mk(k)**2.0d0/(mk(i)+mk(k))**2.0d0  &
             *  .5d0*omega11(i,k)

        diag4(i) = diag4(i)+ nk(i)*nk(k)*mk(k)/(mk(i)+mk(k))**3.0d0 *    &
             ((5.0d0/4.0d0*(6.0d0*mk(i)**2.0d0+5.0d0*mk(k)**2.0d0)*omega11(i,k) &
             - 3.0d0*mk(k)**2.0d0*omega11(i,k) ) & 
             +  4.0d0*mk(i)*mk(k)*omega11(i,k) )

     enddo
  enddo

  do i = 1, nspec
     do j = 1, nspec

        sum00 = -nk(j)*mk(j)*diag2(i) - nk(j)*mk(j)/(mk(i)+mk(j))*nk(i)*mk(i) &
             *omega11(i,j)
        sum01 = - nk(i)*nk(j)*mk(j)**2.0d0/(mk(i)+mk(j))**2.0d0 *  &
             .5d0*omega11(i,j)

        sum11 =  nk(i)*nk(j)*mk(j)/(mk(i)+mk(j))**3.0d0 *         &
             ( -(5.0d0/4.0d0*(6.0d0*mk(j)**2.0d0+5.0d0*mk(j)**2.0d0)*omega11(i,j) &
             - 3.0d0*mk(j)**2.0d0*omega11(i,j) ) & 
             +  4.0d0*mk(j)*mk(j)*omega11(i,j) )

        if(i.eq.j)then
           sum00 = sum00+diag1(i)
           sum01 = sum01+diag3(i)
           sum11 = sum11+diag4(i)
        endif

        ratm = mk(i)/mk(j)
        sqratm = sqrt(ratm)
        scr = ratm*sqratm
        QQ(i,j) =  8.0d0*sqratm/mk(i) * sum00         ! HCB 7.4-50 
        QQ(i,j+nspec) = -8.0d0*scr * sum01    ! HCB 7.4-51
        QQ(i+nspec,j) = -8.d0*sqratm*sum01
        QQ(i+nspec,j+nspec) = 8.0d0*ratm*sqratm* sum11  ! HCB 7.4-53

     enddo
  enddo

  ! Build vector rhs ; see HCB 7.4-54
  do i = 1, nspec
     aSonine(i) = 0.0d0
     beta(i)= sqrt(2.0d0*k_B*T/mk(i))
     aSonine(i+nspec) = -(15.0d0/4.0d0)*nk(i)*beta(i)
  enddo

  ! NOTE: the minus sign below; see HCB p 488
  call decomp(2*nspec,2*nspec,QQ,ip)
  call solve(2*nspec,2*nspec,QQ,aSonine,ip)


  do i = 1, nspec           
     ! HCB 7.4-9
     D_T(i) = 0.5d0*nk(i)*beta(i)*mk(i)*aSonine(i)
  enddo



  sum1=0              
  ! HCB 7.4-33
  do i = 1, nspec
     sum1 = sum1 + nk(i)*beta(i)* aSonine(nspec+i)
  enddo
  lamdaprime = -5.0d0/4.0d0 * k_B * sum1


  sum1 = 0
  ! HCB 7.4-65
  do i = 1, nspec
     do j = 1, nspec 
        sum1 = sum1 + nk(i)*nk(j)/Dbin(i,j) *     & 
             (D_T(i)/(nk(i)*mk(i)) - D_T(j)/(nk(j)*mk(j)))**2.0d0
     enddo
  enddo

  ntotal = sum(nk(:))


  lammix = lamdaprime - 0.5d0*k_B/nTotal * sum1  ! HCB 7.4-65



  return
end subroutine lambda_lin
!-------------------------------------------------------------------------

subroutine thermalDiff(nspec,sigma11,a_ij1,a_ij2,alphabar_in,Xkp,sqrtT,mk,kT)

  implicit none

  integer :: nspec
  real(kind=8), dimension(1:nspec) :: mk,  Xkp
  real(kind=8), dimension(1:nspec,1:nspec) :: alphabar_in
  real(kind=8) :: sqrtT,  kT(1:nspec)
  integer :: i, j, k
  real*8 :: k_B
  real(kind=8), dimension(1:nspec,1:nspec) :: a_ij1, a_ij2,sigma11
  real(kind=8), dimension(1:nspec,1:nspec) :: Aij
  real(kind=8), dimension(1:nspec) :: AA
  real(kind=8) :: fact1, sumTemp, sum1
  real(kind=8), dimension(1:nspec,1:nspec) :: alphaij 
  integer :: ip(1:nspec)



  ! Based on Valk 1963 (Waldmann)

  do i = 1, nspec
     do j = 1, nspec
        Aij(i,j) = Xkp(j)*a_ij2(i,j)
     enddo
     AA(i) = 15.d0/4.d0
  enddo

  do i = 1, nspec
     sumTemp = 0.0d0
     do k = 1, nspec
        sumTemp = sumTemp + Xkp(k)*a_ij1(i,k)
     enddo
     Aij(i,i) = Aij(i,i) + sumTemp
  enddo

  !          call decomp(nspec,nspec,Aij,ip)
  !          call solve(nspec,nspec,Aij,AA,ip)
  call decompnp(nspec,nspec,Aij)
  call solvenp(nspec,nspec,Aij,AA)


  do i = 1, nspec
     do j = 1, nspec

        alphaij(i,j) = alphabar_in(i,j)*(AA(i)/mk(i) - AA(j)/mk(j))/sqrtT

     enddo
  enddo

  do k = 1, nspec
     sumTemp = 0.0d0
     do j = 1, nspec
        sumTemp = sumTemp + Xkp(j)*alphaij(k,j)
     enddo
     kT(k) = sumTemp
  enddo


  return
end subroutine thermalDiff
!-------------------------------------------------
