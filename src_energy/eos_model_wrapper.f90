module eos_model_wrapper_module

  use ml_layout_module
  use eos_model_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: prob_lo, prob_hi
  use probin_energy_module, only: heating_type
  implicit none

  private

  public :: convert_conc_to_molefrac, ideal_mixture_transport_wrapper, &
            add_external_heating, compute_S_alpha, scale_deltaP, &
            compute_h, compute_hk, compute_p, compute_cp

contains

  subroutine convert_conc_to_molefrac(mla,conc,molefrac,conc_to_molefrac)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: conc(:)
    type(multifab) , intent(inout) :: molefrac(:)
    logical        , intent(in   ) :: conc_to_molefrac

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = conc(1)%ng
    ng_2 = molefrac(1)%ng

    do n=1,nlevs
       do i=1,nfabs(conc(n))
          dp1 => dataptr(conc(n), i)
          dp2 => dataptr(molefrac(n), i)
          lo = lwb(get_box(conc(n), i))
          hi = upb(get_box(conc(n), i))
          select case (dm)
          case (2)
             call convert_conc_to_molefrac_2d(dp1(:,:,1,:),ng_1,dp2(:,:,1,:),ng_2, &
                                              lo,hi,conc_to_molefrac)
          case (3)
             call convert_conc_to_molefrac_3d(dp1(:,:,:,:),ng_1,dp2(:,:,:,:),ng_2, &
                                              lo,hi,conc_to_molefrac)
          end select
       end do
    end do

  end subroutine convert_conc_to_molefrac
  
  subroutine convert_conc_to_molefrac_2d(conc,ng_1,molefrac,ng_2,lo,hi,conc_to_molefrac)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::     conc(lo(1)-ng_1:,lo(2)-ng_1:,:)
    real(kind=dp_t), intent(inout) :: molefrac(lo(1)-ng_2:,lo(2)-ng_2:,:)
    logical        , intent(in   ) :: conc_to_molefrac

    ! local
    integer :: i,j
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    if (conc_to_molefrac) then

       do j=lo(2)-ng_2,hi(2)+ng_2
          do i=lo(1)-ng_2,hi(1)+ng_2
             call CKYTX(conc(i,j,:),iwrk,rwrk,molefrac(i,j,:))
          end do
       end do

    else

       do j=lo(2)-ng_1,hi(2)+ng_1
          do i=lo(1)-ng_1,hi(1)+ng_1
             call CKXTY(molefrac(i,j,:),iwrk,rwrk,conc(i,j,:))
          end do
       end do

    end if

  end subroutine convert_conc_to_molefrac_2d
  
  subroutine convert_conc_to_molefrac_3d(conc,ng_1,molefrac,ng_2,lo,hi,conc_to_molefrac)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::     conc(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
    real(kind=dp_t), intent(inout) :: molefrac(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:)
    logical        , intent(in   ) :: conc_to_molefrac

    ! local
    integer :: i,j,k
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    if (conc_to_molefrac) then

       do k=lo(3)-ng_2,hi(3)+ng_2
          do j=lo(2)-ng_2,hi(2)+ng_2
             do i=lo(1)-ng_2,hi(1)+ng_2
                call CKYTX(conc(i,j,k,:),iwrk,rwrk,molefrac(i,j,k,:))
             end do
          end do
       end do

    else

       do k=lo(3)-ng_1,hi(3)+ng_1
          do j=lo(2)-ng_1,hi(2)+ng_1
             do i=lo(1)-ng_1,hi(1)+ng_1
                call CKXTY(molefrac(i,j,k,:),iwrk,rwrk,conc(i,j,k,:))
             end do
          end do
       end do

    end if

  end subroutine convert_conc_to_molefrac_3d

  ! takes $\rho,T,P,\wb$, and $\xb$ as inputs and computes the following:
  ! eta     is the dynamic viscosity, $\eta$
  ! kappa   is the thermal conductivity, $\lambda$
  ! zeta    is the bulk viscosity, $\kappa$
  ! diff_ij is the diffusion matrix, $\chi$
  ! chitil are the thermodiffusion coefficients, $\zeta$.
  !                (I modified the routine to make this so)
  subroutine ideal_mixture_transport_wrapper(mla,rhotot,Temp,p0,Yk,Xk,eta,kappa, &
                                             zeta,diff_ij,chitil)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: p0
    type(multifab) , intent(in   ) :: Yk(:)
    type(multifab) , intent(in   ) :: Xk(:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: zeta(:)
    type(multifab) , intent(inout) :: diff_ij(:)
    type(multifab) , intent(inout) :: chitil(:)

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,ng_8,ng_9
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp4(:,:,:,:)
    real(kind=dp_t), pointer :: dp5(:,:,:,:)
    real(kind=dp_t), pointer :: dp6(:,:,:,:)
    real(kind=dp_t), pointer :: dp7(:,:,:,:)
    real(kind=dp_t), pointer :: dp8(:,:,:,:)
    real(kind=dp_t), pointer :: dp9(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = rhotot(1)%ng
    ng_2 = Temp(1)%ng
    ng_3 = Yk(1)%ng
    ng_4 = Xk(1)%ng
    ng_5 = eta(1)%ng
    ng_6 = kappa(1)%ng
    ng_7 = zeta(1)%ng
    ng_8 = diff_ij(1)%ng
    ng_9 = chitil(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rhotot(n))
          dp1 => dataptr(rhotot(n),i)
          dp2 => dataptr(Temp(n),i)
          dp3 => dataptr(Yk(n),i)
          dp4 => dataptr(Xk(n),i)
          dp5 => dataptr(eta(n),i)
          dp6 => dataptr(kappa(n),i)
          dp7 => dataptr(zeta(n),i)
          dp8 => dataptr(diff_ij(n),i)
          dp9 => dataptr(chitil(n),i)
          lo = lwb(get_box(rhotot(n), i))
          hi = upb(get_box(rhotot(n), i))
          select case (dm)
          case (2)
             call ideal_mixture_transport_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,1),ng_2, &
                                             dp3(:,:,1,:),ng_3,dp4(:,:,1,:),ng_4, &
                                             dp5(:,:,1,1),ng_5,dp6(:,:,1,1),ng_6, &
                                             dp7(:,:,1,1),ng_7,dp8(:,:,1,:),ng_8, &
                                             dp9(:,:,1,:),ng_9,p0,lo,hi)
          case (3)
             call ideal_mixture_transport_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,1),ng_2, &
                                             dp3(:,:,:,:),ng_3,dp4(:,:,:,:),ng_4, &
                                             dp5(:,:,:,1),ng_5,dp6(:,:,:,1),ng_6, &
                                             dp7(:,:,:,1),ng_7,dp8(:,:,:,:),ng_8, &
                                             dp9(:,:,:,:),ng_9,p0,lo,hi)
          end select
       end do
    end do

  end subroutine ideal_mixture_transport_wrapper
  
  subroutine ideal_mixture_transport_2d(rhotot,ng_1,Temp,ng_2,Yk,ng_3,Xk,ng_4,eta,ng_5, &
                                        kappa,ng_6,zeta,ng_7,diff_ij,ng_8,chitil,ng_9, &
                                        p0,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,ng_8,ng_9
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::  rhotot(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(in   ) ::    Temp(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t), intent(in   ) ::      Yk(lo(1)-ng_3:,lo(2)-ng_3:,:)
    real(kind=dp_t), intent(in   ) ::      Xk(lo(1)-ng_4:,lo(2)-ng_4:,:)
    real(kind=dp_t), intent(inout) ::     eta(lo(1)-ng_5:,lo(2)-ng_5:)
    real(kind=dp_t), intent(inout) ::   kappa(lo(1)-ng_6:,lo(2)-ng_6:)
    real(kind=dp_t), intent(inout) ::    zeta(lo(1)-ng_7:,lo(2)-ng_7:)
    real(kind=dp_t), intent(inout) :: diff_ij(lo(1)-ng_8:,lo(2)-ng_8:,:)
    real(kind=dp_t), intent(inout) ::  chitil(lo(1)-ng_9:,lo(2)-ng_9:,:)
    real(kind=dp_t), intent(in   ) :: p0

    ! local
    integer :: i,j,ng

    ng = min(ng_5,ng_6,ng_7,ng_8,ng_9)

    do j=lo(2)-ng,hi(2)+ng
       do i=lo(1)-ng,hi(1)+ng

          call ideal_mixture_transport(rhotot(i,j),Temp(i,j),p0,Yk(i,j,:),Xk(i,j,:), &
                                       eta(i,j),kappa(i,j),zeta(i,j),diff_ij(i,j,:), &
                                       chitil(i,j,:),nspecies)

       end do
    end do

  end subroutine ideal_mixture_transport_2d
  
  subroutine ideal_mixture_transport_3d(rhotot,ng_1,Temp,ng_2,Yk,ng_3,Xk,ng_4,eta,ng_5, &
                                        kappa,ng_6,zeta,ng_7,diff_ij,ng_8,chitil,ng_9, &
                                        p0,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,ng_8,ng_9
    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::  rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(in   ) ::    Temp(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
    real(kind=dp_t), intent(in   ) ::      Yk(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
    real(kind=dp_t), intent(in   ) ::      Xk(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:,:)
    real(kind=dp_t), intent(inout) ::     eta(lo(1)-ng_5:,lo(2)-ng_5:,lo(3)-ng_5:)
    real(kind=dp_t), intent(inout) ::   kappa(lo(1)-ng_6:,lo(2)-ng_6:,lo(3)-ng_6:)
    real(kind=dp_t), intent(inout) ::    zeta(lo(1)-ng_7:,lo(2)-ng_7:,lo(3)-ng_7:)
    real(kind=dp_t), intent(inout) :: diff_ij(lo(1)-ng_8:,lo(2)-ng_8:,lo(3)-ng_8:,:)
    real(kind=dp_t), intent(inout) ::  chitil(lo(1)-ng_9:,lo(2)-ng_9:,lo(3)-ng_9:,:)
    real(kind=dp_t), intent(in   ) :: p0

    ! local
    integer :: i,j,k,ng

    ng = min(ng_5,ng_6,ng_7,ng_8,ng_9)

    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             call ideal_mixture_transport(rhotot(i,j,k),Temp(i,j,k),p0,Yk(i,j,k,:), &
                                          Xk(i,j,k,:),eta(i,j,k),kappa(i,j,k), &
                                          zeta(i,j,k),diff_ij(i,j,k,:), &
                                          chitil(i,j,k,:),nspecies)

          end do
       end do
    end do

  end subroutine ideal_mixture_transport_3d

  subroutine add_external_heating(mla,rhotot,rhoHext,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: rhoHext(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = rhotot(1)%ng
    ng_2 = rhoHext(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rhotot(n))
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(rhoHext(n), i)
          lo = lwb(get_box(rhotot(n), i))
          hi = upb(get_box(rhotot(n), i))
          select case (dm)
          case (2)
             call add_external_heating_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,1),ng_2,dx(n,:),lo,hi)
          case (3)
             call add_external_heating_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,1),ng_2,dx(n,:),lo,hi)
          end select
       end do
    end do

  end subroutine add_external_heating
  
  subroutine add_external_heating_2d(rhotot,ng_1,rhoHext,ng_2,dx,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::  rhotot(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(inout) :: rhoHext(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y,xcen,ycen,r

    select case (heating_type)

    case (1)

       xcen = 0.5d0*(prob_lo(1)+prob_hi(1))
       ycen = 0.5d0*(prob_lo(2)+prob_hi(2))

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(dble(i)+0.5d0)

             r = sqrt((x-xcen)**2 + (y-ycen)**2)
             
             rhoHext(i,j) = rhoHext(i,j) + 1.d10*rhotot(i,j)*exp(-100*r**2)

          end do
       end do

    case default
       ! no external heating
    end select       

  end subroutine add_external_heating_2d
  
  subroutine add_external_heating_3d(rhotot,ng_1,rhoHext,ng_2,dx,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) ::  rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(inout) :: rhoHext(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z,xcen,ycen,zcen,r

    select case (heating_type)

    case (1)

       xcen = 0.5d0*(prob_lo(1)+prob_hi(1))
       ycen = 0.5d0*(prob_lo(2)+prob_hi(2))
       zcen = 0.5d0*(prob_lo(3)+prob_hi(3))

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3)*(dble(k)+0.5d0)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1)*(dble(i)+0.5d0)

                r = sqrt((x-xcen)**2 + (y-ycen)**2 + (z-zcen)**2)
             
                rhoHext(i,j,k) = rhoHext(i,j,k) + 1.d10*rhotot(i,j,k)*exp(-100*r**2)

             end do
          end do
       end do

    case default
       ! no external heating
    end select    

       

  end subroutine add_external_heating_3d

  subroutine compute_S_alpha(mla,S,alpha,mass_fluxdiv,rhoh_fluxdiv,conc,Temp,rhotot,p0)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: S(:)
    type(multifab) , intent(inout) :: alpha(:)
    type(multifab) , intent(in   ) :: mass_fluxdiv(:)
    type(multifab) , intent(in   ) :: rhoh_fluxdiv(:)
    type(multifab) , intent(in   ) :: conc(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    real(kind=dp_t), intent(in   ) :: p0

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp4(:,:,:,:)
    real(kind=dp_t), pointer :: dp5(:,:,:,:)
    real(kind=dp_t), pointer :: dp6(:,:,:,:)
    real(kind=dp_t), pointer :: dp7(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = S(1)%ng
    ng_2 = alpha(1)%ng
    ng_3 = mass_fluxdiv(1)%ng
    ng_4 = rhoh_fluxdiv(1)%ng
    ng_5 = conc(1)%ng
    ng_6 = Temp(1)%ng
    ng_7 = rhotot(1)%ng

    do n=1,nlevs
       do i=1,nfabs(S(n))
          dp1 => dataptr(S(n), i)
          dp2 => dataptr(alpha(n), i)
          dp3 => dataptr(mass_fluxdiv(n), i)
          dp4 => dataptr(rhoh_fluxdiv(n), i)
          dp5 => dataptr(conc(n), i)
          dp6 => dataptr(Temp(n), i)
          dp7 => dataptr(rhotot(n), i)
          lo = lwb(get_box(S(n), i))
          hi = upb(get_box(S(n), i))
          select case (dm)
          case (2)
             call compute_S_alpha_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,1),ng_2, &
                                     dp3(:,:,1,:),ng_3,dp4(:,:,1,1),ng_4, &
                                     dp5(:,:,1,:),ng_5,dp6(:,:,1,1),ng_6, &
                                     dp7(:,:,1,1),ng_7,p0,lo,hi)
          case (3)
             call compute_S_alpha_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,1),ng_2, &
                                     dp3(:,:,:,:),ng_3,dp4(:,:,:,1),ng_4, &
                                     dp5(:,:,:,:),ng_5,dp6(:,:,:,1),ng_6, &
                                     dp7(:,:,:,1),ng_7,p0,lo,hi)
          end select
       end do
    end do

  end subroutine compute_S_alpha
  
  subroutine compute_S_alpha_2d(S,ng_1,alpha,ng_2,mass_fluxdiv,ng_3,rhoh_fluxdiv,ng_4, &
                                conc,ng_5,Temp,ng_6,rhotot,ng_7,p0,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::            S(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(inout) ::        alpha(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t), intent(in   ) :: mass_fluxdiv(lo(1)-ng_3:,lo(2)-ng_3:,:)
    real(kind=dp_t), intent(in   ) :: rhoh_fluxdiv(lo(1)-ng_4:,lo(2)-ng_4:)
    real(kind=dp_t), intent(in   ) ::         conc(lo(1)-ng_5:,lo(2)-ng_5:,:)
    real(kind=dp_t), intent(in   ) ::         Temp(lo(1)-ng_6:,lo(2)-ng_6:)
    real(kind=dp_t), intent(in   ) ::       rhotot(lo(1)-ng_7:,lo(2)-ng_7:)
    real(kind=dp_t), intent(in   ) :: p0

    ! local
    integer :: i,j,n
    real(kind=dp_t) :: W,beta,theta,cpmix,cvmix
    real(kind=dp_t) :: hk(nspecies),P_rho,P_T,P_w(nspecies)

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          ! mixture-averaged molecular mass, W = (sum (w_k/m_k) )^-1
          call CKMMWY(conc(i,j,:),iwrk,rwrk,W)
          ! h_k, c_p, and c_v
          call CKHMS(Temp(i,j),iwrk,rwrk,hk)
          call CKCPBS(Temp(i,j),conc(i,j,:),iwrk,rwrk,cpmix)
          call CKCVBS(Temp(i,j),conc(i,j,:),iwrk,rwrk,cvmix)
          ! P_rho, P_T, P_w
          call compute_P_rho(P_rho,rhotot(i,j),conc(i,j,:),Temp(i,j))
          call compute_P_T(P_T,rhotot(i,j),conc(i,j,:),Temp(i,j))
          call compute_P_w(P_w,rhotot(i,j),conc(i,j,:),Temp(i,j))

          S(i,j) = 0.d0
          do n=1,nspecies
             beta = P_w(n) / (rhotot(i,j)**2*P_rho) - hk(n)*P_T/(rhotot(i,j)**2*cpmix*P_rho)
             S(i,j) = S(i,j) + beta*mass_fluxdiv(i,j,n)
          end do
          theta = P_T / (rhotot(i,j)**2*P_rho*cpmix)
          S(i,j) = S(i,j) + theta*rhoh_fluxdiv(i,j)
          alpha(i,j) = cvmix/(cpmix*p0)

       end do
    end do
       
  end subroutine compute_S_alpha_2d
  
  subroutine compute_S_alpha_3d(S,ng_1,alpha,ng_2,mass_fluxdiv,ng_3,rhoh_fluxdiv,ng_4, &
                                conc,ng_5,Temp,ng_6,rhotot,ng_7,p0,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::            S(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(inout) ::        alpha(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
    real(kind=dp_t), intent(in   ) :: mass_fluxdiv(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
    real(kind=dp_t), intent(in   ) :: rhoh_fluxdiv(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:)
    real(kind=dp_t), intent(in   ) ::         conc(lo(1)-ng_5:,lo(2)-ng_5:,lo(3)-ng_5:,:)
    real(kind=dp_t), intent(in   ) ::         Temp(lo(1)-ng_6:,lo(2)-ng_6:,lo(3)-ng_6:)
    real(kind=dp_t), intent(in   ) ::       rhotot(lo(1)-ng_7:,lo(2)-ng_7:,lo(3)-ng_7:)
    real(kind=dp_t), intent(in   ) :: p0

    ! local
    integer :: i,j,k,n
    real(kind=dp_t) :: W,beta,theta,cpmix,cvmix
    real(kind=dp_t) :: hk(nspecies),P_rho,P_T,P_w(nspecies)

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! mixture-averaged molecular mass, W = (sum (w_k/m_k) )^-1
             call CKMMWY(conc(i,j,k,:),iwrk,rwrk,W)
             ! h_k, c_p, and c_v
             call CKHMS(Temp(i,j,k),iwrk,rwrk,hk)
             call CKCPBS(Temp(i,j,k),conc(i,j,k,:),iwrk,rwrk,cpmix)
             call CKCVBS(Temp(i,j,k),conc(i,j,k,:),iwrk,rwrk,cvmix)
             ! P_rho, P_T, P_w
             call compute_P_rho(P_rho,rhotot(i,j,k),conc(i,j,k,:),Temp(i,j,k))
             call compute_P_T(P_T,rhotot(i,j,k),conc(i,j,k,:),Temp(i,j,k))
             call compute_P_w(P_w,rhotot(i,j,k),conc(i,j,k,:),Temp(i,j,k))

             S(i,j,k) = 0.d0
             do n=1,nspecies
                beta = P_w(n) / (rhotot(i,j,k)**2*P_rho) - hk(n)*P_T/(rhotot(i,j,k)**2*cpmix*P_rho)
                S(i,j,k) = S(i,j,k) + beta*mass_fluxdiv(i,j,k,n)
             end do
             theta = P_T / (rhotot(i,j,k)**2*P_rho*cpmix)
             S(i,j,k) = S(i,j,k) + theta * rhoh_fluxdiv(i,j,k)
             alpha(i,j,k) = cvmix/(cpmix*p0)

          end do
       end do
    end do
       
  end subroutine compute_S_alpha_3d

  subroutine scale_deltaP(mla,deltaP,rhotot,Temp,conc,p0,dt,factor)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: deltaP(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: conc(:)
    real(kind=dp_t), intent(in   ) :: p0,dt,factor

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2,ng_3,ng_4
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp4(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = deltaP(1)%ng
    ng_2 = rhotot(1)%ng
    ng_3 = Temp(1)%ng
    ng_4 = conc(1)%ng

    do n=1,nlevs
       do i=1,nfabs(deltaP(n))
          dp1 => dataptr(deltaP(n), i)
          dp2 => dataptr(rhotot(n), i)
          dp3 => dataptr(Temp(n), i)
          dp4 => dataptr(conc(n), i)
          lo = lwb(get_box(deltaP(n), i))
          hi = upb(get_box(deltaP(n), i))
          select case (dm)
          case (2)
             call scale_deltaP_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,1),ng_2, &
                                  dp3(:,:,1,1),ng_3,dp4(:,:,1,:),ng_4, &
                                  p0,dt,lo,hi,factor)
          case (3)
             call bl_error("scale_deltaP_3d not written yet")
          end select
       end do
    end do

  end subroutine scale_deltaP
  
  subroutine scale_deltaP_2d(deltaP,ng_1,rhotot,ng_2,Temp,ng_3,conc,ng_4, &
                             p0,dt,lo,hi,factor)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: deltaP(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(inout) :: rhotot(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t), intent(in   ) ::   Temp(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t), intent(in   ) ::   conc(lo(1)-ng_4:,lo(2)-ng_4:,:)
    real(kind=dp_t), intent(in   ) :: p0,dt,factor

    ! local
    integer :: i,j
    real(kind=dp_t) :: P_rho

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          call compute_P_rho(P_rho,rhotot(i,j),conc(i,j,:),Temp(i,j))

          deltaP(i,j) = factor * deltaP(i,j) / (rhotot(i,j)*P_rho*dt)

       end do
    end do
       
  end subroutine scale_deltaP_2d

  subroutine compute_h(mla,Temp,h,conc)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(inout) :: h(:)
    type(multifab) , intent(in   ) :: conc(:)

    ! local
    integer :: n,nlevs,comp

    type(multifab) :: hk(mla%nlevel)

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(hk(n),mla%la(n),nspecies,h(n)%ng)
    end do

    ! compute w_k h_k
    call compute_hk(mla,Temp,hk)
    do n=1,nlevs
       call multifab_mult_mult_c(hk(n),1,conc(n),1,nspecies,hk(n)%ng)
    end do

    do n=1,nlevs
       call multifab_setval(h(n),0.d0,all=.true.)
       do comp=1,nspecies
          ! sum (w_k h_k)
          call multifab_plus_plus_c(h(n),1,hk(n),comp,1,h(n)%ng)
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(hk(n))
    end do

  end subroutine compute_h

  subroutine compute_hk(mla,Temp,hk)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(inout) :: hk(:)

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = Temp(1)%ng
    ng_2 = hk(1)%ng

    do n=1,nlevs
       do i=1,nfabs(Temp(n))
          dp1 => dataptr(Temp(n), i)
          dp2 => dataptr(hk(n), i)
          lo = lwb(get_box(Temp(n), i))
          hi = upb(get_box(Temp(n), i))
          select case (dm)
          case (2)
             call compute_hk_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,:),ng_2,lo,hi)
          case (3)
             call compute_hk_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,:),ng_2,lo,hi)
          end select
       end do
    end do

  end subroutine compute_hk
  
  subroutine compute_hk_2d(Temp,ng_1,hk,ng_2,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) :: Temp(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(inout) ::   hk(lo(1)-ng_2:,lo(2)-ng_2:,:)

    ! local
    integer :: i,j
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do j=lo(2)-ng_2,hi(2)+ng_2
       do i=lo(1)-ng_2,hi(1)+ng_2
          call CKHMS(Temp(i,j),iwrk,rwrk,hk(i,j,:))
       end do
    end do

  end subroutine compute_hk_2d
  
  subroutine compute_hk_3d(Temp,ng_1,hk,ng_2,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) :: Temp(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(inout) ::   hk(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:)

    ! local
    integer :: i,j,k
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do k=lo(3)-ng_2,hi(3)+ng_2
       do j=lo(2)-ng_2,hi(2)+ng_2
          do i=lo(1)-ng_2,hi(1)+ng_2
             call CKHMS(Temp(i,j,k),iwrk,rwrk,hk(i,j,k,:))
          end do
       end do
    end do

  end subroutine compute_hk_3d

  subroutine compute_p(mla,rhotot,Temp,conc,pres)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: conc(:)
    type(multifab) , intent(inout) :: pres(:)

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2,ng_3,ng_4
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp4(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = rhotot(1)%ng
    ng_2 = Temp(1)%ng
    ng_3 = conc(1)%ng
    ng_4 = pres(1)%ng

    do n=1,nlevs
       do i=1,nfabs(pres(n))
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(Temp(n), i)
          dp3 => dataptr(conc(n), i)
          dp4 => dataptr(pres(n), i)
          lo = lwb(get_box(pres(n), i))
          hi = upb(get_box(pres(n), i))
          select case (dm)
          case (2)
             call compute_p_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,1),ng_2, &
                               dp3(:,:,1,:),ng_3,dp4(:,:,1,1),ng_4,lo,hi)
          case (3)
             call compute_p_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,1),ng_2, &
                               dp3(:,:,:,:),ng_3,dp4(:,:,:,1),ng_4,lo,hi)
          end select
       end do
    end do

  end subroutine compute_p
  
  subroutine compute_p_2d(rhotot,ng_1,Temp,ng_2,conc,ng_3,pres,ng_4,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) :: rhotot(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(in   ) ::   Temp(lo(1)-ng_2:,lo(2)-ng_2:)
    real(kind=dp_t), intent(in   ) ::   conc(lo(1)-ng_3:,lo(2)-ng_3:,:)
    real(kind=dp_t), intent(inout) ::   pres(lo(1)-ng_4:,lo(2)-ng_4:)

    ! local
    integer :: i,j
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call CKPY(rhotot(i,j),Temp(i,j),conc(i,j,:),iwrk,rwrk,pres(i,j))
       end do
    end do

  end subroutine compute_p_2d
  
  subroutine compute_p_3d(rhotot,ng_1,Temp,ng_2,conc,ng_3,pres,ng_4,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,ng_4,lo(:),hi(:)
    real(kind=dp_t), intent(in   ) :: rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(in   ) ::   Temp(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
    real(kind=dp_t), intent(in   ) ::   conc(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:,:)
    real(kind=dp_t), intent(inout) ::   pres(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:)

    ! local
    integer :: i,j,k
    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             call CKPY(rhotot(i,j,k),Temp(i,j,k),conc(i,j,k,:),iwrk,rwrk,pres(i,j,k))
          end do
       end do
    end do

  end subroutine compute_p_3d

  subroutine compute_cp(mla,cp,conc,Temp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: cp(:)
    type(multifab) , intent(in   ) :: conc(:)
    type(multifab) , intent(in   ) :: Temp(:)

    ! local
    integer :: n,nlevs,i,dm
    integer :: ng_1,ng_2,ng_3
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_1 = cp(1)%ng
    ng_2 = conc(1)%ng
    ng_3 = Temp(1)%ng

    do n=1,nlevs
       do i=1,nfabs(cp(n))
          dp1 => dataptr(cp(n), i)
          dp2 => dataptr(conc(n), i)
          dp3 => dataptr(Temp(n), i)
          lo = lwb(get_box(cp(n), i))
          hi = upb(get_box(cp(n), i))
          select case (dm)
          case (2)
             call compute_cp_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,:),ng_2, &
                                dp3(:,:,1,1),ng_3,lo,hi)
          case (3)
             call compute_cp_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,:),ng_2, &
                                dp3(:,:,:,1),ng_3,lo,hi)
          end select
       end do
    end do

  end subroutine compute_cp
  
  subroutine compute_cp_2d(cp,ng_1,conc,ng_2,Temp,ng_3,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_1:,lo(2)-ng_1:)
    real(kind=dp_t), intent(in   ) :: conc(lo(1)-ng_2:,lo(2)-ng_2:,:)
    real(kind=dp_t), intent(in   ) :: Temp(lo(1)-ng_3:,lo(2)-ng_3:)

    ! local
    integer :: i,j

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do j=lo(2)-ng_1,hi(2)+ng_1
       do i=lo(1)-ng_1,hi(1)+ng_1
          call CKCPBS(Temp(i,j),conc(i,j,:),iwrk,rwrk,cp(i,j))
       end do
    end do       

  end subroutine compute_cp_2d
  
  subroutine compute_cp_3d(cp,ng_1,conc,ng_2,Temp,ng_3,lo,hi)

    integer        , intent(in   ) :: ng_1,ng_2,ng_3,lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::   cp(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
    real(kind=dp_t), intent(in   ) :: conc(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_1:,:)
    real(kind=dp_t), intent(in   ) :: Temp(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_1:)

    ! local
    integer :: i,j,k

    integer :: iwrk
    real(kind=dp_t) :: rwrk

    do k=lo(3)-ng_1,hi(3)+ng_1
       do j=lo(2)-ng_1,hi(2)+ng_1
          do i=lo(1)-ng_1,hi(1)+ng_1
             call CKCPBS(Temp(i,j,k),conc(i,j,k,:),iwrk,rwrk,cp(i,j,k))
          end do
       end do
    end do       

  end subroutine compute_cp_3d
  
end module eos_model_wrapper_module

