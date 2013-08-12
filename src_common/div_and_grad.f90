module div_and_grad_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: compute_grad, compute_div

contains

  subroutine compute_grad(mla,phi,gradp,dx, &
                          start_incomp,start_bccomp,start_outcomp,num_comp,the_bc_level)

    ! compute the face-centered gradient of a cell-centered field

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(inout) :: gradp(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: start_incomp,start_bccomp,start_outcomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,i,dm,nlevs,ng_p,ng_g,comp,bccomp,outcomp
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: gpx(:,:,:,:)
    real(kind=dp_t), pointer :: gpy(:,:,:,:)
    real(kind=dp_t), pointer :: gpz(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_g = gradp(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(phi(n))
          pp  => dataptr(phi(n), i)
          gpx => dataptr(gradp(n,1), i)
          gpy => dataptr(gradp(n,2), i)
          lo = lwb(get_box(phi(n), i))
          hi = upb(get_box(phi(n), i))
          do comp=start_incomp,start_incomp+num_comp-1
             bccomp = start_bccomp + (comp-start_incomp)
             outcomp = start_outcomp + (comp-start_incomp)
             select case (dm)
             case (2)
                call compute_grad_2d(pp(:,:,1,comp), ng_p, &
                                     gpx(:,:,1,outcomp), gpy(:,:,1,outcomp), ng_g, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             case (3)
                gpz => dataptr(gradp(n,3), i)
                call compute_grad_3d(pp(:,:,:,comp), ng_p, &
                                     gpx(:,:,:,outcomp), gpy(:,:,:,outcomp), gpz(:,:,:,outcomp), ng_g, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end select
          end do
       end do
    end do

  contains
    
    subroutine compute_grad_2d(phi,ng_p,gpx,gpy,ng_g,lo,hi,dx,bc)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) :: gpx(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpy(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      ! local
      integer :: i,j

      if(.false.) then
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            print*, phi(i,j)
         end do
      end do
      endif

      ! x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) / dx(1)
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
         i=lo(1)
         do j=lo(2),hi(2)
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) / (0.5d0*dx(1))
         end do
      end if
      if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
         i=hi(1)+1
         do j=lo(2),hi(2)
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) / (0.5d0*dx(1))
         end do
      end if

      ! y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) / dx(2)
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
         j=lo(2)
         do i=lo(1),hi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) / (0.5d0*dx(2))
         end do
      end if

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
         j=hi(2)+1
         do i=lo(1),hi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) / (0.5d0*dx(2))
         end do
      end if

    end subroutine compute_grad_2d

    subroutine compute_grad_3d(phi,ng_p,gpx,gpy,gpz,ng_g,lo,hi,dx,bc)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) :: gpx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpy(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      ! local
      integer :: i,j,k
      
      ! x-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) / dx(1)
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
         i=lo(1)
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) / (0.5d0*dx(1))
            end do
         end do
      end if
      if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
         i=hi(1)+1
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) / (0.5d0*dx(1))
            end do
         end do
      end if

      ! y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) / dx(2)
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
         j=lo(2)
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) / (0.5d0*dx(2))
            end do
         end do
      end if
      if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
         j=hi(2)+1
         do k=lo(3),hi(3)
            do i=lo(1),hi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) / (0.5d0*dx(2))
            end do
         end do
      end if

      ! z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) / dx(3)
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. EXT_DIR) then
         k=lo(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) / (0.5d0*dx(3))
            end do
         end do
      end if

      ! alter stencil at boundary since ghost value represents value at boundary
      if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. EXT_DIR) then
         k=hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) / (0.5d0*dx(3))
            end do
         end do
      end if

    end subroutine compute_grad_3d

  end subroutine compute_grad

  subroutine compute_div(mla,phi_fc,div,dx,start_incomp,start_outcomp,num_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_fc(:,:)
    type(multifab) , intent(inout) :: div(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: start_incomp, start_outcomp, num_comp

    real(kind=dp_t), pointer :: pxp(:,:,:,:) 
    real(kind=dp_t), pointer :: pyp(:,:,:,:) 
    real(kind=dp_t), pointer :: pzp(:,:,:,:) 
    real(kind=dp_t), pointer :: dp(:,:,:,:) 
    integer :: i,n,nlevs,dm,ng_p,ng_d,comp,outcomp
    integer :: lo(mla%dim),hi(mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_p = phi_fc(1,1)%ng
    ng_d = div(1)%ng

    do n = 1,nlevs
       do i = 1, nfabs(div(n))
          pxp => dataptr(phi_fc(n,1), i)
          pyp => dataptr(phi_fc(n,2), i)
          dp => dataptr(div(n), i)
          lo =  lwb(get_box(div(n), i))
          hi =  upb(get_box(div(n), i))
          do comp=start_incomp,start_incomp+num_comp-1
             outcomp = start_outcomp + (comp-start_incomp)
             select case (dm)
             case (2)
                call compute_div_2d(pxp(:,:,1,comp), pyp(:,:,1,comp), ng_p, &
                                    dp(:,:,1,outcomp), ng_d, dx(n,:),lo,hi)
             case (3) 
                pzp => dataptr(phi_fc(n,3), i)
                call compute_div_3d(pxp(:,:,:,comp), pyp(:,:,:,comp), pzp(:,:,:,comp), ng_p, &
                                    dp(:,:,:,outcomp), ng_d, dx(n,:),lo,hi)
             end select
          end do
       end do
    end do

  contains

    subroutine compute_div_2d(phix,phiy,ng_p,div,ng_d,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_d
      real(kind=dp_t), intent(in   ) :: phix(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiy(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) ::  div(lo(1)-ng_d:,lo(2)-ng_d:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            div(i,j) = &
                 (phix(i+1,j) - phix(i,j)) / dx(1) + &
                 (phiy(i,j+1) - phiy(i,j)) / dx(2)
         end do
      end do

    end subroutine compute_div_2d

    subroutine compute_div_3d(phix,phiy,phiz,ng_p,div,ng_d,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_d
      real(kind=dp_t), intent(in   ) :: phix(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiz(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) ::  div(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               div(i,j,k) = &
                    (phix(i+1,j,k) - phix(i,j,k)) / dx(1) + &
                    (phiy(i,j+1,k) - phiy(i,j,k)) / dx(2) + &
                    (phiz(i,j,k+1) - phiz(i,j,k)) / dx(3)
            end do
         end do
      end do

    end subroutine compute_div_3d

  end subroutine compute_div

end module div_and_grad_module
