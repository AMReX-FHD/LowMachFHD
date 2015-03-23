module energy_eos_wrapper_module

  use ml_layout_module
  use energy_EOS_module
  implicit none

  private

  public :: convert_conc_to_molefrac

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
  
end module energy_eos_wrapper_module
