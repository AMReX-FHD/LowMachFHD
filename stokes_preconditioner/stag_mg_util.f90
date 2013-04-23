module stag_mg_util_module

  use ml_layout_module
  use multifab_module

  implicit none

  private

  public :: cc_restriction, stag_restriction, nodal_restriction, edge_restriction, &
       stag_prolongation, stag_mg_update

contains

  ! coarsen a cell-centered quantity
  subroutine cc_restriction(la,phi_c,phi_f,the_bc_level)

    use bl_error_module
    use define_bc_module

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(inout) :: phi_c
    type(multifab), intent(in   ) :: phi_f
    type(bc_level), intent(in   ) :: the_bc_level

    ! local
    integer :: i,dm,ng_c
    integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
    integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

    real(kind=dp_t), pointer :: acp(:,:,:,:)
    real(kind=dp_t), pointer :: afp(:,:,:,:)

    dm = get_dim(la)

    ng_c = phi_c%ng

    if (ng_c .ne. 1) then
       call bl_error("cc_restriction assumes only 1 ghost cell for non-periodic boundary conditions")
    end if

    do i=1,nfabs(phi_c)
       acp => dataptr(phi_c, i)
       afp => dataptr(phi_f, i)
       lo_c = lwb(get_box(phi_c, i))
       hi_c = upb(get_box(phi_c, i))
       lo_f = lwb(get_box(phi_f, i))
       hi_f = upb(get_box(phi_f, i))
       select case (dm)
       case (2)
          call cc_restriction_2d(acp(:,:,1,1),afp(:,:,1,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, &
                                   the_bc_level%adv_bc_level_array(i,:,:,dm+2))
       case (3)
          call cc_restriction_3d(acp(:,:,:,1),afp(:,:,:,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, &
                                   the_bc_level%adv_bc_level_array(i,:,:,dm+2))
       end select
    end do

    call multifab_fill_boundary(phi_c)

  end subroutine cc_restriction

  subroutine cc_restriction_2d(phi_c,phi_f,ng_c,lo_c,hi_c,lo_f,hi_f,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_c
    real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_c:,lo_f(2)-ng_c:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j

    do j=lo_c(2),hi_c(2)
       do i=lo_c(1),hi_c(1)

          phi_c(i,j) = 0.25d0*(  phi_f(2*i,2*j  ) + phi_f(2*i+1,2*j  ) &
                               + phi_f(2*i,2*j+1) + phi_f(2*i+1,2*j+1) )
       end do
    end do

    ! overwrite x-ghost cells
    ! value in ghost cells represents boundary value
    ! average 2 overlying fine values
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do j=lo_c(2),hi_c(2)
          phi_c(lo_c(1)-1,j) = ( phi_f(lo_f(1)-1,2*j) + phi_f(lo_f(1)-1,2*j+1) ) / 2.d0
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do j=lo_c(2),hi_c(2)
          phi_c(hi_c(1)+1,j) = ( phi_f(hi_f(1)+1,2*j) + phi_f(hi_f(1)+1,2*j+1) ) / 2.d0
       end do
    end if

    ! overwrite y-ghost cells
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do i=lo_c(1),hi_c(1)
          phi_c(i,lo_c(2)-1) = ( phi_f(2*i,lo_f(2)-1) + phi_f(2*i+1,lo_f(2)-1) ) / 2.d0
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do i=lo_c(1),hi_c(1)
          phi_c(i,hi_c(2)+1) = ( phi_f(2*i,hi_f(2)+1) + phi_f(2*i+1,hi_f(2)+1) ) / 2.d0
       end do

    end if

  end subroutine cc_restriction_2d

  subroutine cc_restriction_3d(phi_c,phi_f,ng_c,lo_c,hi_c,lo_f,hi_f,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_c
    real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_c:,lo_f(2)-ng_c:,lo_f(3)-ng_c:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j,k

    do k=lo_c(3),hi_c(3)
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)

             phi_c(i,j,k) = 0.125d0*(  phi_f(2*i,2*j  ,2*k  ) + phi_f(2*i+1,2*j  ,2*k  ) &
                                     + phi_f(2*i,2*j+1,2*k  ) + phi_f(2*i+1,2*j+1,2*k  ) &
                                     + phi_f(2*i,2*j  ,2*k+1) + phi_f(2*i+1,2*j  ,2*k+1) &
                                     + phi_f(2*i,2*j+1,2*k+1) + phi_f(2*i+1,2*j+1,2*k+1) )
          end do
       end do
    end do

    ! overwrite x-ghost cells
    ! value in ghost cells represents boundary value
    ! average 4 overlying fine values
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)
             phi_c(lo_c(1)-1,j,k) = ( phi_f(lo_f(1)-1,2*j,2*k)   + phi_f(lo_f(1)-1,2*j+1,2*k  ) &
                                     +phi_f(lo_f(1)-1,2*j,2*k+1) + phi_f(lo_f(1)-1,2*j+1,2*k+1) ) / 4.d0
          end do
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)
             phi_c(hi_c(1)+1,j,k) = ( phi_f(hi_f(1)+1,2*j,2*k)   + phi_f(hi_f(1)+1,2*j+1,2*k  ) &
                                     +phi_f(hi_f(1)+1,2*j,2*k+1) + phi_f(hi_f(1)+1,2*j+1,2*k+1) ) / 4.d0
          end do
       end do
    end if

    ! overwrite y-ghost cells
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do k=lo_c(3),hi_c(3)
          do i=lo_c(1),hi_c(1)
             phi_c(i,lo_c(2)-1,k) = ( phi_f(2*i,lo_f(2)-1,2*k)   + phi_f(2*i+1,lo_f(2)-1,2*k  ) &
                                     +phi_f(2*i,lo_f(2)-1,2*k+1) + phi_f(2*i+1,lo_f(2)-1,2*k+1) ) / 4.d0
          end do
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do k=lo_c(3),hi_c(3)
          do i=lo_c(1),hi_c(1)
             phi_c(i,hi_c(2)+1,k) = ( phi_f(2*i,hi_f(2)+1,2*k)   + phi_f(2*i+1,hi_f(2)+1,2*k  ) &
                                     +phi_f(2*i,hi_f(2)+1,2*k+1) + phi_f(2*i+1,hi_f(2)+1,2*k+1) ) / 4.d0
          end do
       end do
    end if

    ! overwrite z-ghost cells
    if (adv_bc(3,1) .eq. FOEXTRAP .or. adv_bc(3,1) .eq. EXT_DIR) then
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)
             phi_c(i,j,lo_c(3)-1) = ( phi_f(2*i,2*j  ,lo_f(3)-1) + phi_f(2*i+1,2*j  ,lo_f(3)-1) &
                                     +phi_f(2*i,2*j+1,lo_f(3)-1) + phi_f(2*i+1,2*j+1,lo_f(3)-1) ) / 4.d0
          end do
       end do
    end if
    if (adv_bc(3,2) .eq. FOEXTRAP .or. adv_bc(3,2) .eq. EXT_DIR) then
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)
             phi_c(i,j,hi_c(3)+1) = ( phi_f(2*i,2*j  ,hi_f(3)+1) + phi_f(2*i+1,2*j  ,hi_f(3)+1) &
                                     +phi_f(2*i,2*j+1,hi_f(3)+1) + phi_f(2*i+1,2*j+1,hi_f(3)+1) ) / 4.d0
          end do
       end do
    end if

  end subroutine cc_restriction_3d

  ! coarsen a staggered quantity
  subroutine stag_restriction(la,phi_f,phi_c,simple_stencil_in)

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(in   ) :: phi_f(:) ! face-centered
    type(multifab), intent(inout) :: phi_c(:) ! face-centered
    logical, intent(in), optional :: simple_stencil_in

    ! local
    integer :: i,dm,ng_f,ng_c
    integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
    integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

    logical :: simple_stencil

    real(kind=dp_t), pointer :: fpx(:,:,:,:)
    real(kind=dp_t), pointer :: fpy(:,:,:,:)
    real(kind=dp_t), pointer :: fpz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)

    simple_stencil = .false.
    if (present(simple_stencil_in)) then
       simple_stencil = simple_stencil_in
    end if

    dm = get_dim(la)

    ng_f = phi_f(1)%ng
    ng_c = phi_c(1)%ng

    do i=1,nfabs(phi_c(1))
       fpx => dataptr(phi_f(1), i)
       fpy => dataptr(phi_f(2), i)
       cpx => dataptr(phi_c(1), i)
       cpy => dataptr(phi_c(2), i)
       lo_c = lwb(get_box(phi_c(1), i))
       hi_c = upb(get_box(phi_c(1), i))
       lo_f = lwb(get_box(phi_f(1), i))
       hi_f = upb(get_box(phi_f(1), i))
       select case (dm)
       case (2)
          call stag_restriction_2d(fpx(:,:,1,1),fpy(:,:,1,1),ng_f, &
                                   cpx(:,:,1,1),cpy(:,:,1,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, simple_stencil)
       case (3)
          fpz => dataptr(phi_f(3), i)
          cpz => dataptr(phi_c(3), i)
          call stag_restriction_3d(fpx(:,:,:,1),fpy(:,:,:,1),fpz(:,:,:,1),ng_f, &
                                   cpx(:,:,:,1),cpy(:,:,:,1),cpz(:,:,:,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f, simple_stencil)
       end select
    end do

  end subroutine stag_restriction

  subroutine stag_restriction_2d(phix_f,phiy_f,ng_f,phix_c,phiy_c,ng_c, &
                                 lo_c,hi_c,lo_f,hi_f,simple_stencil)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(in   ) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
    real(kind=dp_t), intent(inout) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
    real(kind=dp_t), intent(inout) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
    logical        , intent(in   ) :: simple_stencil

    ! local
    integer :: i,j
    
    if (simple_stencil) then

       ! 2 point stencils
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)+1
             phix_c(i,j) = 0.5d0*(phix_f(2*i,2*j) + phix_f(2*i,2*j+1))
          end do
       end do

       do j=lo_c(2),hi_c(2)+1
          do i=lo_c(1),hi_c(1)
             phiy_c(i,j) = 0.5d0*(phiy_f(2*i,2*j) + phiy_f(2*i+1,2*j))
          end do
       end do

    else

       ! 6 point stencils
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)+1
             phix_c(i,j) = 0.25d0*(phix_f(2*i,2*j) + phix_f(2*i,2*j+1)) &
                        + 0.125d0*( phix_f(2*i+1,2*j) + phix_f(2*i+1,2*j+1) &
                                   +phix_f(2*i-1,2*j) + phix_f(2*i-1,2*j+1))
          end do
       end do

       do j=lo_c(2),hi_c(2)+1
          do i=lo_c(1),hi_c(1)
             phiy_c(i,j) = 0.25d0*(phiy_f(2*i,2*j) + phiy_f(2*i+1,2*j)) &
                        + 0.125d0*( phiy_f(2*i,2*j+1) + phiy_f(2*i+1,2*j+1) &
                                   +phiy_f(2*i,2*j-1) + phiy_f(2*i+1,2*j-1))
          end do
       end do

    end if


  end subroutine stag_restriction_2d

  subroutine stag_restriction_3d(phix_f,phiy_f,phiz_f,ng_f, &
                                 phix_c,phiy_c,phiz_c,ng_c, &
                                 lo_c,hi_c,lo_f,hi_f,simple_stencil)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(in   ) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phiz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: phiz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    logical        , intent(in   ) :: simple_stencil

    ! local
    integer :: i,j,k

    if (simple_stencil) then

       ! 4 point stencils
       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)
             do i=lo_c(1),hi_c(1)+1
                phix_c(i,j,k) = 0.25d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                                         +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) )
             end do
          end do
       end do

       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)+1
             do i=lo_c(1),hi_c(1)
                phiy_c(i,j,k) = 0.25d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                                         +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) )
             end do
          end do
       end do

       do k=lo_c(3),hi_c(3)+1
          do j=lo_c(2),hi_c(2)
             do i=lo_c(1),hi_c(1)
                phiz_c(i,j,k) = 0.25d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                                         +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) )
             end do
          end do
       end do

    else

       ! 12 point stencils
       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)
             do i=lo_c(1),hi_c(1)+1
                phix_c(i,j,k) = 0.125d0* ( phix_f(2*i,2*j,2*k  ) + phix_f(2*i,2*j+1,2*k  ) &
                                          +phix_f(2*i,2*j,2*k+1) + phix_f(2*i,2*j+1,2*k+1) ) &
                               + 0.0625* ( phix_f(2*i+1,2*j,2*k  ) + phix_f(2*i+1,2*j+1,2*k  ) &
                                          +phix_f(2*i+1,2*j,2*k+1) + phix_f(2*i+1,2*j+1,2*k+1) ) &
                               + 0.0625* ( phix_f(2*i-1,2*j,2*k  ) + phix_f(2*i-1,2*j+1,2*k  ) &
                                          +phix_f(2*i-1,2*j,2*k+1) + phix_f(2*i-1,2*j+1,2*k+1) )
             end do
          end do
       end do

       do k=lo_c(3),hi_c(3)
          do j=lo_c(2),hi_c(2)+1
             do i=lo_c(1),hi_c(1)
                phiy_c(i,j,k) = 0.125d0* ( phiy_f(2*i,2*j,2*k  ) + phiy_f(2*i+1,2*j,2*k  ) &
                                          +phiy_f(2*i,2*j,2*k+1) + phiy_f(2*i+1,2*j,2*k+1) ) &
                               + 0.0625* ( phiy_f(2*i,2*j+1,2*k  ) + phiy_f(2*i+1,2*j+1,2*k  ) &
                                          +phiy_f(2*i,2*j+1,2*k+1) + phiy_f(2*i+1,2*j+1,2*k+1) ) &
                               + 0.0625* ( phiy_f(2*i,2*j-1,2*k  ) + phiy_f(2*i+1,2*j-1,2*k  ) &
                                          +phiy_f(2*i,2*j-1,2*k+1) + phiy_f(2*i+1,2*j-1,2*k+1) )
             end do
          end do
       end do

       do k=lo_c(3),hi_c(3)+1
          do j=lo_c(2),hi_c(2)
             do i=lo_c(1),hi_c(1)
                phiz_c(i,j,k) = 0.125d0* ( phiz_f(2*i,2*j  ,2*k) + phiz_f(2*i+1,2*j  ,2*k) &
                                          +phiz_f(2*i,2*j+1,2*k) + phiz_f(2*i+1,2*j+1,2*k) ) &
                             + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k+1) + phiz_f(2*i+1,2*j  ,2*k+1) &
                                          +phiz_f(2*i,2*j+1,2*k+1) + phiz_f(2*i+1,2*j+1,2*k+1) ) &
                             + 0.0625d0* ( phiz_f(2*i,2*j  ,2*k-1) + phiz_f(2*i+1,2*j  ,2*k-1) &
                                          +phiz_f(2*i,2*j+1,2*k-1) + phiz_f(2*i+1,2*j+1,2*k-1) )
             end do
          end do
       end do

    end if

  end subroutine stag_restriction_3d

  ! coarsen a nodal quantity
  subroutine nodal_restriction(la,phi_f,phi_c)

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(in   ) :: phi_f
    type(multifab), intent(inout) :: phi_c

    ! local
    integer :: i,dm,ng_f,ng_c
    integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
    integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    dm = get_dim(la)

    ng_f = phi_f%ng
    ng_c = phi_c%ng

    if (ng_f .ne. 0 .or. ng_c .ne. 0) then
       call bl_error("nodal_restriction assumes 0 ghost cells")
    end if

    do i=1,nfabs(phi_c)
       fp => dataptr(phi_f, i)
       cp => dataptr(phi_c, i)
       lo_c = lwb(get_box(phi_c, i))
       hi_c = upb(get_box(phi_c, i))
       lo_f = lwb(get_box(phi_f, i))
       hi_f = upb(get_box(phi_f, i))
       select case (dm)
       case (2)
          call nodal_restriction_2d(fp(:,:,1,1),ng_f,cp(:,:,1,1),ng_c, &
                                    lo_c, hi_c, lo_f, hi_f)
       case (3)
          call bl_error("as of 2/25/13, 3D does not require nodal_restriction")
          call nodal_restriction_3d(fp(:,:,:,1),ng_f,cp(:,:,:,1),ng_c, &
                                    lo_c, hi_c, lo_f, hi_f)
       end select
    end do

    call multifab_internal_sync(phi_c)

  end subroutine nodal_restriction

  subroutine nodal_restriction_2d(phi_f,ng_f,phi_c,ng_c,lo_c,hi_c,lo_f,hi_f)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
    real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)

    ! local
    integer :: i,j

    do j=lo_c(2),hi_c(2)+1
       do i=lo_c(1),hi_c(1)+1
          phi_c(i,j) = phi_f(2*i,2*j)
       end do
    end do

  end subroutine nodal_restriction_2d

  subroutine nodal_restriction_3d(phi_f,ng_f,phi_c,ng_c,lo_c,hi_c,lo_f,hi_f)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(in   ) :: phi_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: phi_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

    ! local
    integer :: i,j,k

    do k=lo_c(3),hi_c(3)+1
       do j=lo_c(2),hi_c(2)+1
          do i=lo_c(1),hi_c(1)+1
             phi_c(i,j,k) = phi_f(2*i,2*j,2*k)
          end do
       end do
    end do

  end subroutine nodal_restriction_3d

  ! coarsen an edge-based quantity (3d only, no such thing as edges in 2d)
  subroutine edge_restriction(la,phi_f,phi_c)

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(in   ) :: phi_f(:) ! 3 components (xy, xz, yz edges)
    type(multifab), intent(inout) :: phi_c(:) ! 3 components (xy, xz, yz edges)

    ! local
    integer :: i,dm,ng_f,ng_c
    integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
    integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

    real(kind=dp_t), pointer :: fp1(:,:,:,:)
    real(kind=dp_t), pointer :: fp2(:,:,:,:)
    real(kind=dp_t), pointer :: fp3(:,:,:,:)
    real(kind=dp_t), pointer :: cp1(:,:,:,:)
    real(kind=dp_t), pointer :: cp2(:,:,:,:)
    real(kind=dp_t), pointer :: cp3(:,:,:,:)

    dm = get_dim(la)

    ng_f = phi_f(1)%ng
    ng_c = phi_c(1)%ng

    if (ng_f .ne. 0 .or. ng_c .ne. 0) then
       call bl_error("edge_restriction assumes 0 ghost cells")
    end if

    do i=1,nfabs(phi_c(1))
       fp1 => dataptr(phi_f(1), i)
       cp1 => dataptr(phi_c(1), i)
       fp2 => dataptr(phi_f(2), i)
       cp2 => dataptr(phi_c(2), i)
       fp3 => dataptr(phi_f(3), i)
       cp3 => dataptr(phi_c(3), i)
       lo_c = lwb(get_box(phi_c(1), i))
       hi_c = upb(get_box(phi_c(1), i))
       lo_f = lwb(get_box(phi_f(1), i))
       hi_f = upb(get_box(phi_f(1), i))
       select case (dm)
       case (2)
          call bl_error("no such thing as edge_restriction in 2d")
       case (3)
          call edge_restriction_3d(fp1(:,:,:,1),fp2(:,:,:,1),fp3(:,:,:,1),ng_f, &
                                   cp1(:,:,:,1),cp2(:,:,:,1),cp3(:,:,:,1),ng_c, &
                                   lo_c, hi_c, lo_f, hi_f)
       end select
    end do

    do i=1,3
       call multifab_internal_sync(phi_c(i))
    end do

  end subroutine edge_restriction

  subroutine edge_restriction_3d(phixy_f,phixz_f,phiyz_f,ng_f,phixy_c,phixz_c,phiyz_c,ng_c, &
                                 lo_c,hi_c,lo_f,hi_f)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(in   ) :: phixy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phixz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phiyz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: phixy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: phixz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: phiyz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

    ! local
    integer :: i,j,k

    ! xy edges
    do k=lo_c(3),hi_c(3)
       do j=lo_c(2),hi_c(2)+1
          do i=lo_c(1),hi_c(1)+1
             phixy_c(i,j,k) = 0.5d0*(phixy_f(2*i,2*j,2*k)+phixy_f(2*i,2*j,2*k+1))
          end do
       end do
    end do

    ! xz edges
    do k=lo_c(3),hi_c(3)+1
       do j=lo_c(2),hi_c(2)
          do i=lo_c(1),hi_c(1)+1
             phixz_c(i,j,k) =  0.5d0*(phixz_f(2*i,2*j,2*k)+phixz_f(2*i,2*j+1,2*k))
          end do
       end do
    end do

    ! yz edges
    do k=lo_c(3),hi_c(3)+1
       do j=lo_c(2),hi_c(2)+1
          do i=lo_c(1),hi_c(1)
             phiyz_c(i,j,k) =  0.5d0*(phiyz_f(2*i,2*j,2*k)+phiyz_f(2*i+1,2*j,2*k))
          end do
       end do
    end do

  end subroutine edge_restriction_3d

  ! staggered prolongation from coarser grid to finer
  subroutine stag_prolongation(la,phi_f,phi_c)

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(inout) :: phi_f(:) ! face-centered
    type(multifab), intent(in   ) :: phi_c(:) ! face-centered

    ! local
    integer :: i,dm,ng_f,ng_c
    integer :: lo_c(get_dim(la)), hi_c(get_dim(la))
    integer :: lo_f(get_dim(la)), hi_f(get_dim(la))

    real(kind=dp_t), pointer :: fpx(:,:,:,:)
    real(kind=dp_t), pointer :: fpy(:,:,:,:)
    real(kind=dp_t), pointer :: fpz(:,:,:,:)
    real(kind=dp_t), pointer :: cpx(:,:,:,:)
    real(kind=dp_t), pointer :: cpy(:,:,:,:)
    real(kind=dp_t), pointer :: cpz(:,:,:,:)

    dm = get_dim(la)

    ng_f = phi_f(1)%ng
    ng_c = phi_c(1)%ng

    do i=1,nfabs(phi_f(1))
       fpx => dataptr(phi_f(1), i)
       fpy => dataptr(phi_f(2), i)
       cpx => dataptr(phi_c(1), i)
       cpy => dataptr(phi_c(2), i)
       lo_c = lwb(get_box(phi_c(1), i))
       hi_c = upb(get_box(phi_c(1), i))
       lo_f = lwb(get_box(phi_f(1), i))
       hi_f = upb(get_box(phi_f(1), i))
       select case (dm)
       case (2)
          call stag_prolongation_2d(fpx(:,:,1,1),fpy(:,:,1,1),ng_f, &
                                    cpx(:,:,1,1),cpy(:,:,1,1),ng_c, &
                                    lo_c, hi_c, lo_f, hi_f)
       case (3)
          fpz => dataptr(phi_f(3), i)
          cpz => dataptr(phi_c(3), i)
          call stag_prolongation_3d(fpx(:,:,:,1),fpy(:,:,:,1),fpz(:,:,:,1),ng_f, &
                                    cpx(:,:,:,1),cpy(:,:,:,1),cpz(:,:,:,1),ng_c, &
                                    lo_c, hi_c, lo_f, hi_f)
       end select
    end do

  end subroutine stag_prolongation

  subroutine stag_prolongation_2d(phix_f,phiy_f,ng_f,phix_c,phiy_c,ng_c, &
                                  lo_c,hi_c,lo_f,hi_f)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(inout) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
    real(kind=dp_t), intent(inout) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:)

    ! local
    integer :: i,j,ioff,joff

    do j=lo_f(2),hi_f(2)
       do i=lo_f(1),hi_f(1)+1
          
          if (mod(j,2) .eq. 0) then
             joff = -1
          else
             joff = 1
          end if
          
          if (mod(i,2) .eq. 0) then
             
             ! linear interpolation
             phix_f(i,j) = phix_f(i,j) + 0.75d0*phix_c(i/2,j/2) + 0.25d0*phix_c(i/2,j/2+joff)

          else

             ! bilinear interpolation
             phix_f(i,j) = phix_f(i,j) + 0.375d0*phix_c(i/2  ,j/2) &
                                       + 0.125d0*phix_c(i/2  ,j/2+joff) &
                                       + 0.375d0*phix_c(i/2+1,j/2) &
                                       + 0.125d0*phix_c(i/2+1,j/2+joff)

          end if

       end do
    end do
       
    do j=lo_f(2),hi_f(2)+1
       do i=lo_f(1),hi_f(1)

          if (mod(i,2) .eq. 0) then
             ioff = -1
          else
             ioff = 1
          end if

          if (mod(j,2) .eq. 0) then

             ! linear interpolation
             phiy_f(i,j) = phiy_f(i,j) + 0.75d0*phiy_c(i/2,j/2) + 0.25d0*phiy_c(i/2+ioff,j/2)

          else

             ! bilinear interpolation
             phiy_f(i,j) = phiy_f(i,j) + 0.375d0*phiy_c(i/2,j/2  ) &
                                       + 0.125d0*phiy_c(i/2+ioff,j/2  ) &
                                       + 0.375d0*phiy_c(i/2,j/2+1) &
                                       + 0.125d0*phiy_c(i/2+ioff,j/2+1)

             end if

          end do
       end do

  end subroutine stag_prolongation_2d

  subroutine stag_prolongation_3d(phix_f,phiy_f,phiz_f,ng_f,phix_c,phiy_c,phiz_c,ng_c, &
                                  lo_c,hi_c,lo_f,hi_f)

    integer        , intent(in   ) :: lo_c(:),hi_c(:),lo_f(:),hi_f(:),ng_f,ng_c
    real(kind=dp_t), intent(inout) :: phix_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: phiy_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: phiz_f(lo_f(1)-ng_f:,lo_f(2)-ng_f:,lo_f(3)-ng_f:)
    real(kind=dp_t), intent(in   ) :: phix_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: phiy_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: phiz_c(lo_c(1)-ng_c:,lo_c(2)-ng_c:,lo_c(3)-ng_c:)

    ! local
    integer :: i,j,k,ioff,joff,koff

    real(kind=dp_t), parameter :: nine16 = 9.d0/16.d0
    real(kind=dp_t), parameter :: three16 = 3.d0/16.d0
    real(kind=dp_t), parameter :: one16 = 1.d0/16.d0
    real(kind=dp_t), parameter :: nine32 = 9.d0/32.d0
    real(kind=dp_t), parameter :: three32 = 3.d0/32.d0
    real(kind=dp_t), parameter :: one32 = 1.d0/32.d0
    
    do k=lo_f(3),hi_f(3)
       do j=lo_f(2),hi_f(2)
          do i=lo_f(1),hi_f(1)+1

             if (mod(j,2) .eq. 0) then
                joff = -1
             else
                joff = 1
             end if

             if (mod(k,2) .eq. 0) then
                koff = -1
             else
                koff = 1
             end if

             if (mod(i,2) .eq. 0) then
                ! bilinear in the yz plane
                phix_f(i,j,k) = phix_f(i,j,k) &
                     + nine16*phix_c(i/2,j/2,k/2) &
                     + three16*phix_c(i/2,j/2+joff,k/2) &
                     + three16*phix_c(i/2,j/2,k/2+koff) &
                     + one16*phix_c(i/2,j/2+joff,k/2+koff)
             else
                ! bilinear in the yz plane, linear in x
                phix_f(i,j,k) = phix_f(i,j,k) &
                     + nine32*phix_c(i/2,j/2,k/2) &
                     + three32*phix_c(i/2,j/2+joff,k/2) &
                     + three32*phix_c(i/2,j/2,k/2+koff) &
                     + one32*phix_c(i/2,j/2+joff,k/2+koff)&
                     + nine32*phix_c(i/2+1,j/2,k/2) &
                     + three32*phix_c(i/2+1,j/2+joff,k/2) &
                     + three32*phix_c(i/2+1,j/2,k/2+koff) &
                     + one32*phix_c(i/2+1,j/2+joff,k/2+koff)
             end if

          end do
       end do
    end do

    do k=lo_f(3),hi_f(3)
       do j=lo_f(2),hi_f(2)+1
          do i=lo_f(1),hi_f(1)

             if (mod(i,2) .eq. 0) then
                ioff = -1
             else
                ioff = 1
             end if

             if (mod(k,2) .eq. 0) then
                koff = -1
             else
                koff = 1
             end if

             if (mod(j,2) .eq. 0) then
                ! bilinear in the xz plane
                phiy_f(i,j,k) = phiy_f(i,j,k) &
                     + nine16*phiy_c(i/2,j/2,k/2) &
                     + three16*phiy_c(i/2+ioff,j/2,k/2) &
                     + three16*phiy_c(i/2,j/2,k/2+koff) &
                     + one16*phiy_c(i/2+ioff,j/2,k/2+koff)
             else
                ! bilinear in the yz plane, linear in y
                phiy_f(i,j,k) = phiy_f(i,j,k) &
                     + nine32*phiy_c(i/2,j/2,k/2) &
                     + three32*phiy_c(i/2+ioff,j/2,k/2) &
                     + three32*phiy_c(i/2,j/2,k/2+koff) &
                     + one32*phiy_c(i/2+ioff,j/2,k/2+koff)&
                     + nine32*phiy_c(i/2,j/2+1,k/2) &
                     + three32*phiy_c(i/2+ioff,j/2+1,k/2) &
                     + three32*phiy_c(i/2,j/2+1,k/2+koff) &
                     + one32*phiy_c(i/2+ioff,j/2+1,k/2+koff)
             end if

          end do
       end do
    end do

    do k=lo_f(3),hi_f(3)+1
       do j=lo_f(2),hi_f(2)
          do i=lo_f(1),hi_f(1)

             if (mod(i,2) .eq. 0) then
                ioff = -1
             else
                ioff = 1
             end if

             if (mod(j,2) .eq. 0) then
                joff = -1
             else
                joff = 1
             end if

             if (mod(k,2) .eq. 0) then
                ! bilinear in the xy plane
                phiz_f(i,j,k) = phiz_f(i,j,k) &
                     + nine16*phiz_c(i/2,j/2,k/2) &
                     + three16*phiz_c(i/2+ioff,j/2,k/2) &
                     + three16*phiz_c(i/2,j/2+joff,k/2) &
                     + one16*phiz_c(i/2+ioff,j/2+joff,k/2)
             else
                ! bilinear in the xy plane, linear in z
                phiz_f(i,j,k) = phiz_f(i,j,k) &
                     + nine32*phiz_c(i/2,j/2,k/2) &
                     + three32*phiz_c(i/2+ioff,j/2,k/2) &
                     + three32*phiz_c(i/2,j/2+joff,k/2) &
                     + one32*phiz_c(i/2+ioff,j/2+joff,k/2)&
                     + nine32*phiz_c(i/2,j/2,k/2+1) &
                     + three32*phiz_c(i/2+ioff,j/2,k/2+1) &
                     + three32*phiz_c(i/2,j/2+joff,k/2+1) &
                     + one32*phiz_c(i/2+ioff,j/2+joff,k/2+1)
             end if

          end do
       end do
    end do

  end subroutine stag_prolongation_3d

  ! finish the Jacobi iteration by multiplying the residual by the inverse
  ! of the diagonal-element-only matrix
  subroutine stag_mg_update(la,phi_fc,rhs_fc,Lphi_fc, &
                            alpha_fc,beta_cc,beta_nd,beta_ed,gamma_cc,dx,color_in)
    
    use probin_module, only: visc_type

    type(layout)  , intent(in   ) :: la
    type(multifab), intent(inout) :: phi_fc(:)   ! face-centered
    type(multifab), intent(in   ) :: rhs_fc(:)   ! face-centered
    type(multifab), intent(in   ) :: Lphi_fc(:)  ! face-centered
    type(multifab), intent(in   ) :: alpha_fc(:) ! face-centered
    type(multifab), intent(in   ) :: beta_cc     ! cell-centered
    type(multifab), intent(in   ) :: beta_nd     ! nodal
    type(multifab), intent(in   ) :: beta_ed(:)  ! edge-based
    type(multifab), intent(in   ) :: gamma_cc    ! cell-centered
    real(kind=dp_t),intent(in   ) :: dx(:)
    integer        , intent(in   ), optional :: color_in
    
    ! local
    integer :: i,dm,ng_p,ng_r,ng_l,ng_a,ng_b,ng_n,ng_g,ng_e
    integer :: lo(get_dim(la)), hi(get_dim(la))
    integer :: color
 
    real(kind=dp_t), pointer :: ppx(:,:,:,:)
    real(kind=dp_t), pointer :: ppy(:,:,:,:)
    real(kind=dp_t), pointer :: ppz(:,:,:,:)
    real(kind=dp_t), pointer :: rpx(:,:,:,:)
    real(kind=dp_t), pointer :: rpy(:,:,:,:)
    real(kind=dp_t), pointer :: rpz(:,:,:,:)

    real(kind=dp_t), pointer :: lpx(:,:,:,:)
    real(kind=dp_t), pointer :: lpy(:,:,:,:)
    real(kind=dp_t), pointer :: lpz(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:)
    real(kind=dp_t), pointer :: apy(:,:,:,:)
    real(kind=dp_t), pointer :: apz(:,:,:,:)
    real(kind=dp_t), pointer ::  bp(:,:,:,:)
    real(kind=dp_t), pointer :: bp1(:,:,:,:)
    real(kind=dp_t), pointer :: bp2(:,:,:,:)
    real(kind=dp_t), pointer :: bp3(:,:,:,:)
    real(kind=dp_t), pointer :: bnp(:,:,:,:)
    real(kind=dp_t), pointer ::  kp(:,:,:,:)

    dm = get_dim(la)

    if (present(color_in)) then
       color = color_in
    else
       color = 0
    end if

    ng_p = phi_fc(1)%ng
    ng_r = rhs_fc(1)%ng
    ng_l = Lphi_fc(1)%ng
    ng_a = alpha_fc(1)%ng
    ng_b = beta_cc%ng
    ng_g = gamma_cc%ng

    do i=1,nfabs(Lphi_fc(1))
       ppx => dataptr(phi_fc(1), i)
       ppy => dataptr(phi_fc(2), i)
       rpx => dataptr(rhs_fc(1), i)
       rpy => dataptr(rhs_fc(2), i)
       lpx => dataptr(Lphi_fc(1), i)
       lpy => dataptr(Lphi_fc(2), i)
       apx => dataptr(alpha_fc(1), i)
       apy => dataptr(alpha_fc(2), i)
       bp  => dataptr(beta_cc, i)
       kp  => dataptr(gamma_cc, i)
       lo = lwb(get_box(Lphi_fc(1), i))
       hi = upb(get_box(Lphi_fc(1), i))
       select case(dm)
       case (2)
          ng_n = beta_nd%ng
          bnp => dataptr(beta_nd, i)
          call stag_mg_update_2d(ppx(:,:,1,1),ppy(:,:,1,1),ng_p, &
                                 rpx(:,:,1,1),rpy(:,:,1,1),ng_r, &
                                 lpx(:,:,1,1),lpy(:,:,1,1),ng_l, &
                                 apx(:,:,1,1),apy(:,:,1,1),ng_a, &
                                 bp(:,:,1,1),ng_b,bnp(:,:,1,1),ng_n, &
                                 kp(:,:,1,1), ng_g,lo,hi,dx,color)
       case (3)
          ng_e = beta_ed(1)%ng
          ppz => dataptr(phi_fc(3), i)
          rpz => dataptr(rhs_fc(3), i)
          lpz => dataptr(Lphi_fc(3), i)
          apz => dataptr(alpha_fc(3), i)
          bp1 => dataptr(beta_ed(1), i)
          bp2 => dataptr(beta_ed(2), i)
          bp3 => dataptr(beta_ed(3), i)
          call stag_mg_update_3d(ppx(:,:,:,1),ppy(:,:,:,1),ppz(:,:,:,1),ng_p, &
                                 rpx(:,:,:,1),rpy(:,:,:,1),rpz(:,:,:,1),ng_r, &
                                 lpx(:,:,:,1),lpy(:,:,:,1),lpz(:,:,:,1),ng_l, &
                                 apx(:,:,:,1),apy(:,:,:,1),apz(:,:,:,1),ng_a, &
                                 bp(:,:,:,1),ng_b, &
                                 bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                                 kp(:,:,:,1),ng_g,lo,hi,dx,color)
          
       end select
    end do

  end subroutine stag_mg_update

  subroutine stag_mg_update_2d(phix,phiy,ng_p,rhsx,rhsy,ng_r, &
                               Lpx,Lpy,ng_l,alphax,alphay,ng_a,beta,ng_b, &
                               beta_nd,ng_n,gamma,ng_g,lo,hi,dx,color)

    use probin_module, only: visc_type, stag_mg_omega

    integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_r,ng_l,ng_a,ng_b,ng_n,ng_g
    real(kind=dp_t), intent(inout) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(inout) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    rhsx(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(in   ) ::    rhsy(lo(1)-ng_r:,lo(2)-ng_r:)
    real(kind=dp_t), intent(in   ) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(in   ) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_nd(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: color

    ! local
    integer :: i,j

    real(kind=dp_t) :: fac, dxsq, fourthirds
    real(kind=dp_t) :: b,c

    ! coloring parameters
    logical :: do_x, do_y
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       offset = 2
    end if

    dxsq = dx(1)**2
    fourthirds = 4.d0/3.d0
    
    if (visc_type .eq. -1) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j) + &
                     (beta(i,j)+beta(i-1,j)+beta_nd(i,j)+beta_nd(i,j+1))/dxsq

                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j) + &
                     (beta(i,j)+beta(i,j-1)+beta_nd(i,j)+beta_nd(i+1,j))/dxsq

                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2))

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j) + 4.d0*b/dxsq
                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j) + 4.d0*b/dxsq
                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_nd(i,j)+beta_nd(i,j+1))/dxsq

                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_nd(i,j)+beta_nd(i+1,j))/dxsq

                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2))

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j)+6.d0*b/dxsq
                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j)+6.d0*b/dxsq
                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j) &
                     +fourthirds*beta(i-1,j)+gamma(i-1,j) &
                     +beta_nd(i,j)+beta_nd(i,j+1))/dxsq

                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j) &
                     +fourthirds*beta(i,j-1)+gamma(i,j-1) &
                     +beta_nd(i,j)+beta_nd(i+1,j))/dxsq

                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2))
       c = gamma(lo(1),lo(2))

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                fac = alphax(i,j) + (14.d0/3.d0*b+2.d0*c)/dxsq
                phix(i,j) = phix(i,j) + stag_mg_omega*(rhsx(i,j)-Lpx(i,j)) / fac

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                fac = alphay(i,j) + (14.d0/3.d0*b+2.d0*c)/dxsq
                phiy(i,j) = phiy(i,j) + stag_mg_omega*(rhsy(i,j)-Lpy(i,j)) / fac

             end do
          end do

       end if

    end if

  end subroutine stag_mg_update_2d

  subroutine stag_mg_update_3d(phix,phiy,phiz,ng_p,rhsx,rhsy,rhsz,ng_r, &
                               Lpx,Lpy,Lpz,ng_l,alphax,alphay,alphaz,ng_a,beta,ng_b, &
                               beta_xy,beta_xz,beta_yz,ng_e,gamma,ng_g, &
                               lo,hi,dx,color)

    use probin_module, only: visc_type, stag_mg_omega

    integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_r,ng_l,ng_a,ng_b,ng_e,ng_g
    real(kind=dp_t), intent(inout) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(inout) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(inout) ::    phiz(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    rhsx(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) ::    rhsy(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) ::    rhsz(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    real(kind=dp_t), intent(in   ) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(in   ) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(in   ) ::     Lpz(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphaz(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: color

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: fac, dxsq, fourthirds
    real(kind=dp_t) :: b,c

    ! coloring parameters
    logical :: do_x, do_y, do_z
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    do_z = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 5 .or. color .eq. 6) then
       do_x = .false.
       do_y = .false.
       offset = 2
    end if
    
    dxsq = dx(1)**2
    fourthirds = 4.d0/3.d0

    if (visc_type .eq. -1) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( beta(i,j,k)+beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + 6.d0*b/dxsq
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + 6.d0*b/dxsq
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + 6.d0*b/dxsq
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + 8.d0*b/dxsq
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + 8.d0*b/dxsq
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + 8.d0*b/dxsq
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2),lo(3))
       c = gamma(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   fac = alphax(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphay(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

                   fac = alphaz(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq
                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac

                end do
             end do
          end do

       end if

    end if

  end subroutine stag_mg_update_3d

end module stag_mg_util_module
