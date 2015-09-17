module average_to_faces_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use probin_reactdiff_module, only: avg_type

  implicit none

  private

  public :: average_to_faces

contains



  subroutine average_to_faces(mla,n_cc,n_fc,incomp,outcomp,numcomp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_cc(:)
    type(multifab) , intent(inout) :: n_fc(:,:)
    integer        , intent(in   ) :: incomp,outcomp,numcomp

    integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
    integer :: ng_c,ng_f

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: fx(:,:,:,:)
    real(kind=dp_t), pointer :: fy(:,:,:,:)
    real(kind=dp_t), pointer :: fz(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_c = n_cc(1)%ng
    ng_f = n_fc(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(n_cc(n))
          cp => dataptr(n_cc(n),i)
          fx => dataptr(n_fc(n,1),i)
          fy => dataptr(n_fc(n,2),i)
          lo = lwb(get_box(n_cc(n),i))
          hi = upb(get_box(n_cc(n),i))
          select case (dm)
          case (2)
             call average_to_faces_2d(cp(:,:,1,:),ng_c, &
                                      fx(:,:,1,:),fy(:,:,1,:),ng_f, lo,hi, &
                                      incomp,outcomp,numcomp)
          case (3)
             fz => dataptr(n_fc(n,3),i)
             call average_to_faces_3d(cp(:,:,:,:),ng_c, &
                                      fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),ng_f, lo,hi, &
                                      incomp,outcomp,numcomp)
          end select
       end do
    end do

    ! sync the n_fces at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(n_fc(n,i))  
       end do
    end do

  end subroutine average_to_faces

  subroutine average_to_faces_2d(n_cc,ng_c,n_fcx,n_fcy,ng_f,lo,hi,incomp,outcomp,numcomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_f,incomp,outcomp,numcomp
    real(kind=dp_t), intent(in   ) ::   n_cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t), intent(inout) ::  n_fcx(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  n_fcy(lo(1)-ng_f:,lo(2)-ng_f:,:)

    integer :: i,j,comp

    do comp=0,numcomp-1

       if (avg_type .eq. 1) then
          ! arithmetic averaging

          ! x-faces
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                n_fcx(i,j,outcomp+comp) = 0.5d0*(n_cc(i-1,j,incomp+comp) + n_cc(i,j,incomp+comp))
             end do
          end do

          ! y-faces
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                n_fcy(i,j,outcomp+comp) = 0.5d0*(n_cc(i,j-1,incomp+comp) + n_cc(i,j,incomp+comp))
             end do
          end do

       else if (avg_type .eq. 2) then
          ! geometric averaging

          ! x-faces
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                n_fcx(i,j,outcomp+comp) = sqrt(max(n_cc(i-1,j,incomp+comp),0.d0)*max(n_cc(i,j,incomp+comp),0.d0))
             end do
          end do

          ! y-faces
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                n_fcy(i,j,outcomp+comp) = sqrt(max(n_cc(i,j-1,incomp+comp),0.d0)*max(n_cc(i,j,incomp+comp),0.d0))
             end do
          end do

       else if (avg_type .eq. 3) then
          ! harmonic averaging

          ! x-faces
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                n_fcx(i,j,outcomp+comp) = 2.d0 / (1.d0/n_cc(i-1,j,incomp+comp) + 1.d0/n_cc(i,j,incomp+comp))
             end do
          end do

          ! y-faces
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                n_fcy(i,j,outcomp+comp) = 2.d0 / (1.d0/n_cc(i,j-1,incomp+comp) + 1.d0/n_cc(i,j,incomp+comp))
             end do
          end do

       else
          call bl_error("average_to_faces_2d: invalid avg_type")
       end if

    end do

  end subroutine average_to_faces_2d

  subroutine average_to_faces_3d(n_cc,ng_c,n_fcx,n_fcy,n_fcz,ng_f,lo,hi,incomp,outcomp,numcomp)

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_f,incomp,outcomp,numcomp
    real(kind=dp_t), intent(in   ) ::   n_cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(inout) ::  n_fcx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  n_fcy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  n_fcz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)

    integer :: i,j,k,comp

    do comp=0,numcomp-1

       if (avg_type .eq. 1) then
          ! arithmetic averaging

          ! x-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   n_fcx(i,j,k,outcomp+comp) = 0.5d0*(n_cc(i-1,j,k,incomp+comp) + n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

          ! y-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   n_fcy(i,j,k,outcomp+comp) = 0.5d0*(n_cc(i,j-1,k,incomp+comp) + n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

          ! z-faces
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   n_fcz(i,j,k,outcomp+comp) = 0.5d0*(n_cc(i,j,k-1,incomp+comp) + n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

       else if (avg_type .eq. 2) then
          ! geometric averaging

          ! x-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   n_fcx(i,j,k,outcomp+comp) = sqrt(max(n_cc(i-1,j,k,incomp+comp),0.d0)*max(n_cc(i,j,k,incomp+comp),0.d0))
                end do
             end do
          end do

          ! y-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   n_fcy(i,j,k,outcomp+comp) = sqrt(max(n_cc(i,j-1,k,incomp+comp),0.d0)*max(n_cc(i,j,k,incomp+comp),0.d0))
                end do
             end do
          end do

          ! z-faces
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   n_fcz(i,j,k,outcomp+comp) = sqrt(max(n_cc(i,j,k-1,incomp+comp),0.d0)*max(n_cc(i,j,k,incomp+comp),0.d0))
                end do
             end do
          end do

       else if (avg_type .eq. 3) then
          ! harmonic averaging

          ! x-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1
                   n_fcx(i,j,k,outcomp+comp) = 2.d0 / (1.d0/n_cc(i-1,j,k,incomp+comp) + 1.d0/n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

          ! y-faces
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)
                   n_fcy(i,j,k,outcomp+comp) = 2.d0 / (1.d0/n_cc(i,j-1,k,incomp+comp) + 1.d0/n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

          ! z-faces
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   n_fcz(i,j,k,outcomp+comp) = 2.d0 / (1.d0/n_cc(i,j,k-1,incomp+comp) + 1.d0/n_cc(i,j,k,incomp+comp))
                end do
             end do
          end do

       else
          call bl_error("average_to_faces_3d: invalid avg_type")
       end if

    end do

  end subroutine average_to_faces_3d

end module average_to_faces_module
