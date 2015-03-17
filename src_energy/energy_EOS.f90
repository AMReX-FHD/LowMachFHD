module energy_EOS_module

  use bl_types
  use bl_IO_module
  use bl_constants_module
  use probin_common_module, only: k_B, Runiv, molmass
  use probin_multispecies_module, only: nspecies

  implicit none

  integer, parameter  :: MAXSPECIES = 20

  REAL*8, SAVE :: AVOGADRO

  REAL*8, SAVE  :: dia_in(MAXSPECIES)          ! diameters of molecules
  REAL*8, SAVE  :: int_deg_free_in(MAXSPECIES) ! internal degrees of freedom
  integer, SAVE :: use_fake_diff = 0           ! force diffusion coefficients to constant
  REAL*8, SAVE  :: fake_diff_coeff = -1.d0     ! Adjust transport coefficients artificially
  REAL*8, SAVE  :: fake_soret_factor = 1.0d0   ! Adjust transport coefficients artificially

  NAMELIST /probin_energy_EOS/ dia_in
  NAMELIST /probin_energy_EOS/ int_deg_free_in
  NAMELIST /probin_energy_EOS/ use_fake_diff
  NAMELIST /probin_energy_EOS/ fake_diff_coeff
  NAMELIST /probin_energy_EOS/ fake_soret_factor

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

  subroutine energy_EOS_init()

!    use problem_setup

    integer :: iwrk, nfit, i, ic, ii, j
    integer :: dochem, dostrang
    double precision :: rwrk
    integer ns
    real*8 :: mu,Fij,Fijstar,fact1
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs

    narg = command_argument_count()

    ! default values to be replace from inputs file or command line
    dia_in = 0.d0
    int_deg_free_in = 0
    use_fake_diff = 0
    fake_diff_coeff = -1.d0
    fake_soret_factor = 1.d0

    ! read from input file 
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_energy_EOS)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--use_fake_diff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_fake_diff

       case ('--fake_diff_coeff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fake_diff_coeff

       case ('--fake_soret_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fake_soret_factor

       case ('--')
          farg = farg + 1
          exit

       case default

       end select
       farg = farg + 1
    end do

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
!      R_g(ns) = Runiv / molecular_weight(ns)
       cvgas(ns) = 0.5d0*(3+int_deg_free(ns))*Runiv/molecular_weight(ns)
       cpgas(ns) = 0.5d0*(5+int_deg_free(ns))*Runiv/molecular_weight(ns)

    enddo
 
    write(6,*)"e0fref",e0ref
    write(6,*)"dia",dia
    write(6,*)"molecular_weight",molecular_weight
    write(6,*)"molecular_mass",molmass
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

  end subroutine energy_EOS_init

  subroutine CKCPBS(temp,Yk,IWRK,RWRK,cpmix)

      real(kind=8) :: temp,Yk(1:nspecies),RWRK,cpmix
      integer :: IWRK, ns

      cpmix = 0.0d0
      do ns = 1, nspecies
         cpmix = cpmix + Yk(ns)*cpgas(ns)
      enddo
 
  return
  end subroutine


  subroutine CKCVBS(temp,Yk,IWRK,RWRK,cvmix)

      real(kind=8) :: temp,Yk(1:nspecies),RWRK,cvmix
      integer :: IWRK, ns

      cvmix = 0.0d0
      do ns = 1, nspecies
         cvmix = cvmix + Yk(ns)*cvgas(ns)
      enddo
 
  return
  end subroutine


                 subroutine CKUBMS(temp,Yk,IWRK,RWRK,eintmix)

               real(kind=8) :: temp,Yk(1:nspecies),RWRK,eintmix
               integer :: IWRK, ns
               real(kind=8) :: cvmix, e0                                  

               cvmix = 0.0d0; e0 = 0.0d0
               do ns = 1, nspecies
                cvmix = cvmix + Yk(ns)*cvgas(ns)
                e0 = e0 + Yk(ns)*e0ref(ns)
               enddo
 
                 eintmix = e0 + cvmix*temp 
 
                 if(temp.ge.6000.0d0) then
                   print*, 'BUG IN CKUBMS ', temp,Yk  
                 endif
 

                 return
                 end subroutine

!------------------------------------------------------------------------------------------

                 subroutine get_t_given_ey(eintmix, Yk, IWRK, RWRK, temp,ierr)     

               real(kind=8) :: temp,Yk(1:nspecies),RWRK,eintmix
               integer :: IWRK, ns , ierr
               real(kind=8) :: cvmix, e0                                  

               cvmix = 0.0d0; e0 = 0.0d0
               do ns = 1, nspecies
                cvmix = cvmix + Yk(ns)*cvgas(ns)
                e0 = e0 + Yk(ns)*e0ref(ns)
               enddo
 
                 temp = (eintmix-e0)/cvmix 

                 if(temp.ge.6000.0d0) then
                   print*, 'BUG IN feeytt ', eintmix, Yk, cvmix  
                   stop
                 endif
 

                 return
                 end subroutine

!------------------------------------------------------------------------------------------

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
                 end subroutine

!------------------------------------------------------------------------------------------

                 subroutine ckhms(T, iwrk, rwrk, hk) 

             real(kind=8) :: T, rwrk, hk(1:nspecies)
             integer :: iwrk, ns   


             do ns = 1, nspecies 
              hk(ns) = e0ref(ns) + cpgas(ns)*T
             enddo


                 return
                 end subroutine
!------------------------------------------------------------------------------------------

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
                 end subroutine

!------------------------------------------------------------------------------------------

                 subroutine CKYTCR(rho,T,Yk, iwrk, rwrk, Ck)

             real(kind=8) :: Yk(1:nspecies), Ck(1:nspecies), rwrk
             real(kind=8) :: rho,T
             integer :: iwrk, ns
             real*8 :: molmix

             do ns = 1, nspecies
               Ck(ns) = rho*Yk(ns)/molecular_weight(ns)
             enddo 


                 return
                 end subroutine

!------------------------------------------------------------------------------------------

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
                 end subroutine
!------------------------------------------------------------------------------------------

                 subroutine CKMMWY(Yk, iwrk, rwrk, molmix)

             real(kind=8) :: Yk(1:nspecies), rwrk, molmix
             integer :: iwrk, ns

             molmix = 0.0d0
             do ns = 1, nspecies
               molmix = molmix + Yk(ns)/molecular_weight(ns)
             enddo 
             molmix = 1.0d0/molmix


                 return
                 end subroutine

!------------------------------------------------------------------------------------------

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
                 end subroutine

!------------------------------------------------------------------------------------------

                 subroutine CKCPMS(T, iwrk, rwrk, cp) 

            real(kind=8) :: T, rwrk, cp(1:nspecies)
            integer :: iwrk, ns
                            
               do ns = 1, nspecies
                cp(ns) = cpgas(ns)
               enddo 
 

                 return
                 end subroutine

!------------------------------------------------------------------------------------------

   subroutine ideal_mixture_transport(density,temperature,pressure,Yk,Xk,eta,kappa,zeta,diff_ij,chitil,nspec)

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


   end subroutine


end module energy_EOS_module
   
!----------------------------------------------------------------------

          subroutine PGAMMAP(nspec,Ykp,Xkp,Dbin,D_tilde,wk,wbar)

           implicit none 
 
          integer :: nspec, i, j, k
          real(kind=8) :: wk(1:nspec),wbar
          real(kind=8) :: Ykp(1:nspec), Xkp(1:nspec)
          real(kind=8) :: Dbin(1:nspec,1:nspec), term1, term2
          real(kind=8) :: D_tilde(1:nspec,1:nspec), Di(1:nspec), GAM(1:nspec,1:nspec)
          real(kind=8), dimension(1:nspec,1:nspec) :: Pmat, Pmat1


          ! Find Di matrix 
          do i = 1, nspec
           term1 = 0.0d0 
           term2 = 0.0d0
           do j = 1, nspec
            if(j.ne.i) then
              term1 = term1 + Ykp(j)
              term2 = term2 + Xkp(j)/Dbin(i,j)
            endif
           enddo   
           Di(i) = term1/term2 
          enddo   


          ! Compute GAMMA
          GAM = 0.0d0
          do i = 1, nspec
           GAM(i,i) = wk(i)/wbar*Di(i)           
          enddo


          ! Compute P matrix
          Pmat = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
             Pmat(i,j) = - Ykp(i) 
             if(i.eq.j) then
              Pmat(i,j) =  Pmat(i,j) + 1.0d0  
             endif
           enddo
          enddo


          ! Compute Pmat*GAM and store it in Pmat1
          Pmat1 = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             Pmat1(i,j) = Pmat1(i,j) + Pmat(i,k)*GAM(k,j)
            enddo
           enddo
          enddo

          ! Compute Pmat1*P and store it in D_tilde
          D_tilde = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             D_tilde(i,j) = D_tilde(i,j) + Pmat1(i,k)*Pmat(k,j)
            enddo
           enddo
          enddo




        return
        end subroutine

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
        end subroutine



!----------------------------------------------------------------------

          subroutine D_GIOVANGIGLI_old(nspec,Dbin,Ykp,Xkp,D_tilde)

           implicit none 
 
          integer :: nspec, i, j, k, jj
          real(kind=8) :: mk(1:nspec), dk(1:nspec), wk(1:nspec)
          real(kind=8) :: Yk(1:nspec), Xk(1:nspec), Ykp(1:nspec), Xkp(1:nspec)
          real(kind=8) :: dij(1:nspec,1:nspec), T, rho, pres, mw, eps1
          real(kind=8) :: Dbin(1:nspec,1:nspec), term1, term2
          real(kind=8) :: D_tilde(1:nspec,1:nspec), Di(1:nspec), Diff_ij(1:nspec,1:nspec)
          real(kind=8) :: Deltamat(1:nspec,1:nspec), Zmat(1:nspec,1:nspec)
          real(kind=8), dimension(1:nspec,1:nspec) :: Pmat, Jmat, Minv, Mmat
          real(kind=8), dimension(1:nspec,1:nspec) :: PJ, matrix1, matrix2, matrix3, PJJ
          integer :: jmax


          jmax = 3 

         
          eps1 = 1.0d-16

          ! Find Di matrix 
          do i = 1, nspec
           term1 = 0.0d0 
           term2 = 0.0d0
           do j = 1, nspec
            if(j.ne.i) then
              term1 = term1 + Ykp(j)
              term2 = term2 + Xkp(j)/Dbin(i,j)
            endif
           enddo   
           Di(i) = term1/term2 
          enddo   


          ! Compute Mmat and Minv
          Mmat = 0.0d0
          Minv = 0.0d0
          do i = 1, nspec
           Mmat(i,i) = Xkp(i)/Di(i)
           Minv(i,i) = Di(i)/Xkp(i)
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
           enddo
          enddo  


          ! Compute Zmat
          do i = 1, nspec
           do j = 1, nspec
            Zmat(i,j) = Mmat(i,j) - Deltamat(i,j) 
           enddo 
          enddo  


          ! Compute Jmat
          Jmat = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             Jmat(i,j) = Jmat(i,j) + Minv(i,k)*Zmat(k,j)
            enddo
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
          matrix1 = 0.0d0
          matrix2 = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             matrix1(i,j) = matrix1(i,j) + Pmat(i,k)*Minv(k,j)
            enddo
           enddo
          enddo
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             matrix2(i,j) = matrix2(i,j) + matrix1(i,k)*Pmat(j,k) 
                ! notice the change in indices for Pmat to represent Pmat^t
            enddo
           enddo
          enddo


          ! Initialize PJJ with the identity matrix
          PJJ = 0.0d0
          do i = 1, nspec
           PJJ(i,i) = 1.0d0
          enddo 

          jj = 0
          
          ! matrix3 is the summation of the matrix multiplications
          matrix3 = 0.0d0

701       continue 
          matrix1 = 0.0d0 

          ! Compute PJJ*matrix2; sore it in matrix1
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             matrix1(i,j) = matrix1(i,j) + PJJ(i,k)*matrix2(k,j)
            enddo
           enddo
          enddo 
 

          matrix3 = matrix3 + matrix1
           
          jj = jj + 1 

          if(jj.le.jmax) then

          ! Re-evaluate PJJ; matrix1 is a scratch array 
          matrix1 = 0.0d0 
          do i = 1, nspec
           do j = 1, nspec
            do k = 1, nspec
             matrix1(i,j) = matrix1(i,j) + PJJ(i,k)*PJ(k,j)
            enddo
           enddo
          enddo          
          PJJ = matrix1

            goto 701
          endif




          ! Compute Diff_ij
          Diff_ij = matrix3


          ! Compute D_tilde
          D_tilde = 0.0d0
          do i = 1, nspec
           do j = 1, nspec
             D_tilde(i,j) = Diff_ij(i,j)*Ykp(i)
           enddo
          enddo   


        return
        end subroutine


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
          end subroutine

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
         end subroutine
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
          end subroutine 


!-------------------------------------------------

          subroutine viscosityM(xoriginal,T,m,k_B,d,etaMix,nspeci)

            use bl_constants_module

         integer :: nspeci
         real(kind=8), dimension(1:nspeci) :: x, xoriginal, m
         real(kind=8), dimension(1:nspeci,1:nspeci) :: d

         real(kind=8), dimension(1:nspeci) :: eta1
         real(kind=8), dimension(1:nspeci,1:nspeci) :: eta2,H
         real(kind=8), dimension(1:nspeci+1,1:nspeci+1) :: N
         real(kind=8) :: T,  sumTemp, etaMix
         integer :: ii, jj, kk
         real(kind=8) :: numerator, denominator
         real(kind=8) :: eps, sumx
         real*8 :: k_B, avo

         real*8 :: DMGT

         avo = 6.0246d23
         ! clear variables
         H(1:nspeci,1:nspeci) = 0.0d0
         N(1:nspeci+1,1:nspeci+1) = 0.0d0
         eta1(1:nspeci) = 0.0d0
         eta2(1:nspeci,1:nspeci) = 0.0d0
            
 
             x = xoriginal
             ! DANGER
           do ii = 1, nspeci
            if(x(ii).lt.1.0d-6) then
              x(ii) = x(ii) + 1.0d-6 
            endif
           enddo           
           sumx = sum(x(:))  
           x(:) = x(:)/sumx

      do ii=1,nspeci
        ! page 528 HCB
       eta1(ii) = 5.0d0/(16.0d0*d(ii,ii)**2.0d0)*sqrt(m(ii)*k_B*T/M_PI)

      enddo


      do ii=1,nspeci
       do jj=1,ii
        ! page 529 HCB
         eta2(ii,jj) = 5.0d0/(16.0d0*d(ii,jj)**2.0d0)*   &
             sqrt(2.0d0*m(ii)*m(jj)/(m(ii)+m(jj))*k_B*T/M_PI)
        eta2(jj,ii) = eta2(ii,jj)
       enddo
      enddo

    
      do ii=1,nspeci
       do jj=1,nspeci
        if( ii .eq. jj ) then
            sumTemp = 0.0d0
            do kk=1,nspeci
             if( kk .ne. ii ) then
               sumTemp = sumTemp + 2.0d0*x(ii)*x(kk)/eta2(ii,kk)    &
                * (m(ii)*m(kk)/(m(ii)+m(kk))**2.0d0)                &
                * (5.0d0/3.0d0 + m(kk)/m(ii))
             endif
            enddo
            H(ii,ii) = (x(ii)**2.0d0)/eta1(ii) + sumTemp
        else
            H(ii,jj) = -2.0d0*x(ii)*x(jj)/eta2(ii,jj)                  &
               * (m(ii)*m(jj)/(m(ii)+m(jj))**2.0d0) * (2.0d0/3.0d0)
        endif
       enddo
      enddo


       N(1:nspeci,1:nspeci) = H(1:nspeci,1:nspeci)
       N(1:nspeci,nspeci+1) = x(1:nspeci)
       N(nspeci+1,1:nspeci) = x(1:nspeci)
       N(nspeci+1,nspeci+1) = 0.0d0

      eps = 1.0d-200 ! PRECISION 
      numerator = DMGT(eps,nspeci+1,N)
      denominator = DMGT(eps,nspeci,H)

      etaMix = -numerator/denominator


      return
          end subroutine

!-------------------------------------------------

          subroutine lambdaM(xoriginal,T,m,k_B,d,Cv,lamMix,nspeci)

           use bl_constants_module

         integer :: nspeci
         real(kind=8), dimension(1:nspeci) :: x, xoriginal, m, Cv
         real(kind=8), dimension(1:nspeci,1:nspeci) :: d

         real(kind=8), dimension(1:nspeci) :: lam1
         real(kind=8), dimension(1:nspeci,1:nspeci) :: lam2
         real(kind=8), dimension(1:nspeci,1:nspeci) :: L00,L01,L10,L11
         real(kind=8), dimension(1:2*nspeci,1:2*nspeci) :: L
         real(kind=8), dimension(1:2*nspeci+1,1:2*nspeci+1) :: N
         real(kind=8) :: T, sumTemp, lamMix
         integer :: ii, jj, kk
         real(kind=8) :: numerator, denominator
         real(kind=8) :: eps, sumx
         real*8 :: k_B, avo, matscale

         real*8 :: DMGT

         ! clear variables
         L00(1:nspeci,1:nspeci) = 0.0d0
         L01(1:nspeci,1:nspeci) = 0.0d0
         L10(1:nspeci,1:nspeci) = 0.0d0
         L11(1:nspeci,1:nspeci) = 0.0d0
         L(1:2*nspeci,1:2*nspeci) = 0.0d0
         N(1:2*nspeci+1,1:2*nspeci+1) = 0.0d0

         lam1(1:nspeci) = 0.0d0
         lam2(1:nspeci,1:nspeci) = 0.0d0

         avo = 6.0246d23

            x = xoriginal 
             ! DANGER
           do ii = 1, nspeci
            if(x(ii).lt.1.0d-6) then
              x(ii) = x(ii) + 1.0d-6 
            endif
           enddo           
           sumx = sum(x(:))  
           x(:) = x(:)/sumx




      matscale = 0.d0
      do ii=1,nspeci
         ! p 534, HC!
        lam1(ii) = 75.0d0/(64.0d0*d(ii,ii)**2.0d0)* k_B*  &
             sqrt(k_B*T/m(ii)/M_PI)
        matscale = max(matscale,lam1(ii))
      enddo


           ! p 535, HCB
      do ii=1,nspeci
       do jj=1,nspeci
         lam2(ii,jj) = 75.0d0/(64.0d0*d(ii,jj)**2.0d0)*k_B*   &
             sqrt((m(ii)+m(jj))/(2.0d0*m(ii)*m(jj))*k_B*T/M_PI)
       enddo
      enddo

           ! p 538, HCB
      do ii=1,nspeci
       do jj=1,nspeci
        if( ii .eq. jj ) then
          L00(ii,ii) = 0.0d0
        else
            sumTemp = 0.0d0
            do kk=1,nspeci
             if( kk .ne. ii ) then
               sumTemp = sumTemp + 2.0d0*x(ii)*x(kk)/lam2(ii,kk)   &
                * m(jj)/m(ii)
             endif
            enddo
            L00(ii,jj) = 2.0d0*x(ii)*x(jj)/lam2(ii,jj) + sumTemp
        endif
       enddo
      enddo


      do ii=1,nspeci
       do jj=1,nspeci
        if( ii .eq. jj ) then
            sumTemp = 0.0d0
            do kk=1,nspeci
             if( kk .ne. ii ) then
               sumTemp = sumTemp + x(ii)*x(kk)*m(kk)*(1.0d0/5.0d0)    &
        /(m(ii)+m(kk))/lam2(ii,kk)
             endif
            enddo
          L01(ii,ii) = 5.0d0*sumTemp
        else
          L01(ii,jj) = -5.0d0*x(ii)*x(jj)*m(ii)*(1.0d0/5.0d0)   &
                     /(m(ii)+m(jj))/lam2(ii,jj)
        endif
       enddo
      enddo


      do ii=1,nspeci
       do jj=1,nspeci
         L10(ii,jj) = m(jj)/m(ii)*L01(ii,jj)
       enddo
      enddo

      
      do ii=1,nspeci
       do jj=1,nspeci
        if( ii .eq. jj ) then
            sumTemp = 0.0d0
            do kk=1,nspeci
             if( kk .ne. ii ) then
               sumTemp = sumTemp + 2.0d0*x(ii)*x(kk)*(15.0d0/2.0d0*m(ii)**2.0d0   &
                    + 13.0d0/4.0d0*m(kk)**2.0d0 + 4.0d0*m(ii)*m(kk))             &
                    /((m(ii)+m(kk))**2.0d0)/lam2(ii,kk)
             endif
            enddo
            L11(ii,jj) = -4.0d0*x(ii)*x(ii)/lam1(ii) - sumTemp
        else
            L11(ii,jj) = 2.0d0*x(ii)*x(jj)*m(ii)*m(jj)     &
                  /((m(ii)+m(jj))**2.0d0)/lam2(ii,jj)*(27.0d0/4.0d0)
        endif
       enddo
      enddo


       L(1:nspeci,1:nspeci) = L00(1:nspeci,1:nspeci)
       L(1:nspeci,nspeci+1:2*nspeci) = L01(1:nspeci,1:nspeci)
       L(nspeci+1:2*nspeci,1:nspeci) = L10(1:nspeci,1:nspeci)
       L(nspeci+1:2*nspeci,nspeci+1:2*nspeci) = L11(1:nspeci,1:nspeci)


       N(1:2*nspeci,1:2*nspeci) = L(1:2*nspeci,1:2*nspeci)
       N(1:nspeci,2*nspeci+1) = 0.0d0
       N(nspeci+1:2*nspeci,2*nspeci+1) = x(1:nspeci)
       N(2*nspeci+1,1:nspeci) = 0.0d0
       N(2*nspeci+1,nspeci+1:2*nspeci) = x(1:nspeci)
       N(2*nspeci+1,2*nspeci+1) = 0.0d0

       L = L*matscale
       N = N*matscale

      eps = 1.0d-200 ! PRECISION 
      numerator = DMGT(eps,2*nspeci+1,N)
      denominator = DMGT(eps,2*nspeci,L)

      lamMix = 4.0d0*numerator/(denominator*matscale)


      return
          end subroutine


!-------------------------------------------------


        subroutine diffusivityM(x,y,T,m,k_B,d,pres,Diff,nspeci)

           use bl_constants_module

        integer :: nspeci
         real(kind=8), dimension(1:nspeci) :: x, y, m
         real(kind=8), dimension(1:nspeci,1:nspeci) :: d
         real(kind=8), dimension(1:nspeci,1:nspeci) :: Diff    &
                        , Pmat, Pmatt, Mmat, D2, D3

         real(kind=8), dimension(1:nspeci,1:nspeci) :: D1

         real(kind=8) :: T, pres,  sumTemp, MWmix, sumy
         integer :: ii, jj, kk, ll, nn
         real(kind=8) :: deter
         real(kind=8) :: eps
         real(kind=8), dimension(1:nspeci) :: Dimix, FDV
         real*8 :: k_B


         ! clear variables
         D1(1:nspeci,1:nspeci) = 0.0d0
         Pmat(1:nspeci,1:nspeci) = 0.0d0
         Pmatt(1:nspeci,1:nspeci) = 0.0d0

           ! Binary diffusivity
         do ii = 1, nspeci
          do jj = 1, nspeci
           ! p 539, HCB
             D1(ii,jj) = (3.0d0/16.0d0)*sqrt(2*M_PI*(k_B**3.0d0)*   &
                      (T**3.0d0)*(m(ii)+m(jj))/m(ii)/m(jj))         &
                     /(pres*M_PI*(d(ii,jj)**2.0d0))
          enddo
         enddo

           
           MWmix = 0.0d0
             do ii = 1, nspeci
              MWmix = MWmix + y(ii)/m(ii)
             enddo
             MWmix = 1.0d0/MWmix


             Dimix(1:nspeci) = 0.0d0
             do ii = 1, nspeci
              sumTemp = 0.0d0
              sumy = 0.0d0
              do jj = 1, nspeci
                if(jj.ne.ii) then
                  sumy = sumy + y(jj)
                  sumTemp = sumTemp + x(jj)/D1(ii,jj)
                endif
              enddo
              Dimix(ii) = sumy/sumTemp
             enddo

             FDV(1:nspeci) = 0.0d0
             do ii = 1, nspeci
              FDV(ii) = m(ii)/MWmix*Dimix(ii)
             enddo

             Pmat(1:nspeci,1:nspeci) = 0.0d0
             do ii = 1, nspeci
              do jj = 1, nspeci
               Pmat(ii,jj) = -y(ii)/sum(y(:))
               if(jj.eq.ii) then
                Pmat(ii,jj) = Pmat(ii,jj) + 1.0d0
               endif
              enddo
             enddo


             Pmatt(1:nspeci,1:nspeci) = 0.0d0
             do ii = 1, nspeci
              do jj = 1, nspeci
                Pmatt(ii,jj) = Pmat(jj,ii)
              enddo
             enddo


             Mmat(1:nspeci,1:nspeci) = 0.0d0
             do ii = 1, nspeci
               Mmat(ii,ii) = FDV(ii)
             enddo


             D2(1:nspeci,1:nspeci) = 0.0d0
            ! multiply Pmat and Mmat and store it in D2
           do ii = 1, nspeci
            do jj = 1, nspeci
             do ll = 1, nspeci
              D2(ii,jj) = D2(ii,jj) + Pmat(ii,ll)*Mmat(ll,jj)
             enddo
            enddo
           enddo

            D3(1:nspeci,1:nspeci) = 0.0d0
            ! multiply current D2 and Pmatt and store it in D3
           do ii = 1, nspeci
            do jj = 1, nspeci
             do ll = 1, nspeci
              !D3(ii,jj) = D3(ii,jj) + D2(ii,ll)*Pmatt(ll,jj)
              D3(ii,jj) = D3(ii,jj) + D2(ii,ll)*Pmat(ll,jj)
             enddo
            enddo
           enddo
           ! D3 is Dtilde in EGLIB


           Diff(1:nspeci,1:nspeci) = D3(1:nspeci,1:nspeci)
           ! Diff is Dtilde in EGLIB 


          return
        end subroutine
      

!-------------------------------------------------

        subroutine thermalDiffM(x,T,m,k_B,d,kT,nspeci)

           use bl_constants_module

         integer :: nspeci
         real(kind=8), dimension(1:nspeci) :: x, m
         real(kind=8), dimension(1:nspeci) :: kT
         real(kind=8), dimension(1:nspeci,1:nspeci) :: d
         real*8 :: k_B

         real(kind=8) :: T,  sumTemp
         integer :: ii, jj, kk
         real(kind=8) :: eps, Runiv, sigma11, sigma12, Fij, Fijstar

         real(kind=8), dimension(1:nspeci,1:nspeci) :: a_ij1, a_ij2
         real(kind=8), dimension(1:nspeci,1:nspeci) :: Aij, Aij_U
         real(kind=8), dimension(1:nspeci) :: Bi
         real(kind=8), dimension(1:nspeci) :: AA
         real(kind=8) :: fact1, sum1
         real(kind=8), dimension(1:nspeci,1:nspeci) :: alphaij



         ! Based on Valk 1963 (Waldmann)
         
         a_ij1(1:nspeci,1:nspeci) = 0.0d0
         a_ij2(1:nspeci,1:nspeci) = 0.0d0
         do ii = 1, nspeci
          do jj = 1, nspeci
           sigma11 = sqrt( k_B*T/(2.0d0*M_PI*m(ii)*m(jj)/  &
                      (m(ii)+m(jj))) )*M_PI*(d(ii,jj)**2.0d0)
           Fij = (6.0d0*m(ii)*m(ii) + 13.0d0/5.0d0*m(jj)*m(jj) +   &
                  16.0d0/5.0d0*m(ii)*m(jj))/((m(ii)+m(jj))**2.0d0)
           Fijstar = -27.0d0/5.0d0

           a_ij1(ii,jj) = 5.0d0/(k_B*T)*m(ii)*m(jj)/    &
                      (m(ii)+m(jj))*sigma11*Fij
           a_ij2(ii,jj) = 5.0d0/(k_B*T)*m(ii)*m(jj)*m(ii)*m(jj)/    &
                      ((m(ii)+m(jj))**3.0d0)*sigma11*Fijstar
          enddo
         enddo


         Aij(1:nspeci,1:nspeci) = 0.0d0
         do ii = 1, nspeci
          do jj = 1, nspeci
           Aij(ii,jj) = x(jj)*a_ij2(ii,jj)
          enddo
          Bi(ii) = 15.0d0/4.0d0
         enddo

         do ii = 1, nspeci
           sumTemp = 0.0d0
           do kk = 1, nspeci
            sumTemp = sumTemp + x(kk)*a_ij1(ii,kk)
           enddo
           Aij(ii,ii) = Aij(ii,ii) + sumTemp
         enddo



          Aij_U = Aij
          do jj = 1, nspeci-1
           do ii = jj+1, nspeci
             fact1 = Aij_U(ii,jj)/Aij_U(jj,jj)
              do kk = 1, nspeci
                Aij_U(ii,kk) = Aij_U(ii,kk) - fact1*Aij_U(jj,kk)
              enddo
              Bi(ii) = Bi(ii) - fact1*Bi(jj)
           enddo
          enddo

         AA(1:nspeci) = 0.0d0
         do ii = nspeci, 1, -1
          sum1 = 0
           do jj = nspeci, ii+1, -1
             sum1 = sum1 + Aij_U(ii,jj)*AA(jj)
           enddo
           AA(ii) = (Bi(ii)-sum1)/Aij_U(ii,ii)
         enddo


          do ii = 1, nspeci
           do jj = 1, nspeci
            fact1 = m(ii)*m(jj)/(m(ii)+m(jj))
           sigma11 = sqrt( k_B*T/(2.0d0*M_PI*m(ii)*m(jj)/   &
                      (m(ii)+m(jj))) )*M_PI*(d(ii,jj)**2.0d0)
            sigma12 = 3.0d0*sigma11

            alphaij(ii,jj) = 8.0d0/(3.0d0*k_B*T)*fact1*fact1*   &
        (5.0d0/2.0d0*sigma11 - sigma12)*(AA(ii)/m(ii) - AA(jj)/m(jj))
           enddo
          enddo

          do kk = 1, nspeci
           sumTemp = 0.0d0
           do jj = 1, nspeci
            sumTemp = sumTemp + x(jj)*alphaij(kk,jj)
           enddo
!          kT(kk) = x(kk)*sumTemp
! jbb want to return chi_tilde
           kT(kk) = sumTemp
          enddo

            
      return
          end subroutine

!-------------------------------------------------
        Subroutine TSRGT(eps, nn, AA, it, CCC, Kp1, Lp)
       real(kind=8) :: eps
       integer :: nn,it
       real(kind=8) :: AA(nn,nn), CCC(nn,nn)
       integer :: Kp1(nn),Lp(nn)
       real(kind=8)  :: po,t0
       integer :: ii,jj, kk

      CCC=AA; it=1; kk=1
      do while (it==1.and.kk<nn)
       po=CCC(kk,kk); lo=kk; ko=kk
       do ii=kk, nn
        do jj=kk, nn
        if (dabs(CCC(ii,jj))>dabs(po)) then
          po=CCC(ii,jj); lo=ii; ko=jj
        end if
      end do
      end do
      Lp(kk)=lo; Kp1(kk)=ko
      if (dabs(po)<eps) then
      it=0
      else
      if (lo.ne.kk) then
        do jj=kk, nn
          t0=CCC(kk,jj); CCC(kk,jj)=CCC(lo,jj); CCC(lo,jj)=t0
        end do
      end if
      if (ko.ne.kk) then
        do ii=1, nn
          t0=CCC(ii,kk); CCC(ii,kk)=CCC(ii,ko); CCC(ii,ko)=t0
        end do
      end if
      do ii=kk+1, nn
        CCC(ii,kk)=CCC(ii,kk)/po
        do jj=kk+1, nn
          CCC(ii,jj)=CCC(ii,jj)-CCC(ii,kk)*CCC(kk,jj)
        end do
      end do
      kk=kk+1
      end if
      end do
      if (it==1.and.dabs(CCC(nn,nn))<eps)  it=0
      return
      End !TSRGT


      double precision Function DMGT(eps, nn, AA)
       real(kind=8) :: eps, AA(nn,nn)
       real(kind=8) :: d0

       real(kind=8), pointer :: CCC(:,:)
       integer,pointer :: Kp1(:), Lp(:)

!allocate local matrix C and vectors Kp, Lp
         allocate(CCC(nn,nn),STAT=ialloc)
         allocate(Kp1(nn),STAT=ialloc)
         allocate(Lp(nn),STAT=ialloc)

       call TSRGT(eps,nn,AA,it,CCC,Kp1,Lp)  !call triangularization subroutine
       if (it==0) then
        d0=0.d0  !matrix singular, det=0
       else       !matrix regular, det<>0
       d0=1.d0
         do kk=1, nn
          d0=d0*CCC(kk,kk)
        end do
       ll=0
        do kk=1, nn-1
      if (Lp(kk).ne.kk)  ll=ll+1
      if (Kp1(kk).ne.kk)  ll=ll+1
       end do
      if (MOD(ll,2).ne.0) d0=-d0  !l is odd
      end if
      DMGT=d0   !return determinant

       deallocate(CCC,Kp1,Lp) 

      return
      End

