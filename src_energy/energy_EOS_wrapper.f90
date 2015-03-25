module energy_eos_wrapper_module

  use ml_layout_module
  use energy_EOS_module
  use probin_multispecies_module, only: nspecies
  implicit none

  private

  public :: convert_conc_to_molefrac, ideal_mixture_transport_wrapper

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
                                             dp3(:,:,1,:),ng_3,dp4(:,:,1,:),ng_2, &
                                             dp5(:,:,1,1),ng_5,dp6(:,:,1,1),ng_2, &
                                             dp7(:,:,1,1),ng_7,dp8(:,:,1,:),ng_2, &
                                             dp9(:,:,1,:),ng_9,p0,lo,hi)
          case (3)
             call ideal_mixture_transport_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,1),ng_2, &
                                             dp3(:,:,:,:),ng_3,dp4(:,:,:,:),ng_2, &
                                             dp5(:,:,:,1),ng_5,dp6(:,:,:,1),ng_2, &
                                             dp7(:,:,:,1),ng_7,dp8(:,:,:,:),ng_2, &
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
  
end module energy_eos_wrapper_module

