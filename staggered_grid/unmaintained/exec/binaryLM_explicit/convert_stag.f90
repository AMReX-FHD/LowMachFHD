module convert_module

  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: convert_m_to_umac, convert_cons_to_prim, shift_face_to_cc, average_face_to_cc, &
       average_cc_to_face, average_cc_to_node, average_cc_to_edge
  
contains

  subroutine convert_m_to_umac(mla,s_face,m,umac,m_to_umac)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) ::    s_face(:,:)
    type(multifab) , intent(inout) ::    m(:,:)  
    type(multifab) , intent(inout) :: umac(:,:)
    logical        , intent(in   ) :: m_to_umac

    ! local
    integer :: n,i,dm,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    if (m_to_umac) then

       ! compute umac = m / rho in valid region
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac(n,i), 1, m(n,i), 1, 1, 0)
             call multifab_div_div_c(umac(n,i), 1, s_face(n,i), 1, 1, 0)
          end do
       end do

    else

       ! compute m = rho * umac in valid plus ghost region
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(m(n,i), 1, umac(n,i), 1, 1, 1)
             call multifab_mult_mult_c(m(n,i), 1, s_face(n,i), 1, 1, 1)
          end do
       end do

    end if

  end subroutine convert_m_to_umac

  subroutine convert_cons_to_prim(mla,s,prim,cons_to_prim)
    
    use probin_module, only: nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) ::    s(:)
    type(multifab) , intent(inout) :: prim(:)
    logical        , intent(in   ) :: cons_to_prim

    ! local
    integer :: nlevs,n,i

    nlevs = mla%nlevel

    if (cons_to_prim) then

       ! cons to prim - excluding ghost cells
       do n=1,nlevs
          call multifab_copy_c(prim(n),1,s(n),1,nscal,0)
          do i=2,nscal
             call multifab_div_div_c(prim(n),i,s(n),1,1,0)
          end do
       end do

    else

       ! prim to cons - including ghost cells
       do n=1,nlevs
          call multifab_copy_c(s(n),1,prim(n),1,nscal,s(n)%ng)
          do i=2,nscal
             call multifab_mult_mult_c(s(n),i,s(n),1,1,s(n)%ng)
          end do
       end do

    end if
    
  end subroutine convert_cons_to_prim

  subroutine shift_face_to_cc(mla,face,face_comp,cc,cc_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: face(:)
    type(multifab) , intent(inout) ::   cc(:)
    integer        , intent(in   ) :: face_comp,cc_comp

    ! local
    integer :: n,i,dm,nlevs,ng_f,ng_c
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_f = face(1)%ng
    ng_c = cc(1)%ng

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp => dataptr(cc(n), i)
          ep => dataptr(face(n), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          select case (dm)
          case (2)
             call shift_face_to_cc_2d(ep(:,:,1,:), ng_f, cp(:,:,1,:), ng_c, &
                                      face_comp, cc_comp, lo, hi)
          case (3)
             call shift_face_to_cc_3d(ep(:,:,:,:), ng_f, cp(:,:,:,:), ng_c, &
                                      face_comp, cc_comp, lo, hi)
          end select
       end do
    end do

  contains
    
    subroutine shift_face_to_cc_2d(face,ng_f,cc,ng_c,face_comp,cc_comp,lo,hi)

      integer        , intent(in   ) :: ng_f,ng_c,face_comp,cc_comp,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,:)

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            cc(i,j,cc_comp) = face(i,j,face_comp)
         end do
      end do

    end subroutine shift_face_to_cc_2d

    subroutine shift_face_to_cc_3d(face,ng_f,cc,ng_c,face_comp,cc_comp,lo,hi)

      integer        , intent(in   ) :: ng_f,ng_c,face_comp,cc_comp,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k,cc_comp) = face(i,j,k,face_comp)
            end do
         end do
      end do

    end subroutine shift_face_to_cc_3d

  end subroutine shift_face_to_cc

  ! A. Donev
  subroutine average_face_to_cc(mla,face,face_comp,cc,cc_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: face(:)
    type(multifab) , intent(inout) ::   cc(:)
    integer        , intent(in   ) :: face_comp,cc_comp

    ! local
    integer :: n,i,dm,nlevs,ng_f,ng_c
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: av_dim ! Along which dimension to do the average

    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_f = face(1)%ng
    ng_c = cc(1)%ng

    do n=1,nlevs
       if(count(face(n)%nodal)/=1) then
          stop "average_face_to_cc passed wrong type of face multifab"
       else if(face(n)%nodal(1)) then
          av_dim=1
       else if(face(n)%nodal(2)) then
          av_dim=2
       else if(face(n)%nodal(3)) then
          av_dim=3            
       end if

       do i=1,nfabs(cc(n))
          cp => dataptr(cc(n), i)
          ep => dataptr(face(n), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          select case (dm)
          case (2)
             call average_face_to_cc_2d(ep(:,:,1,:), ng_f, cp(:,:,1,:), ng_c, &
                                       face_comp, cc_comp, lo, hi, av_dim)
          case (3)
             call average_face_to_cc_3d(ep(:,:,:,:), ng_f, cp(:,:,:,:), ng_c, &
                                       face_comp, cc_comp, lo, hi, av_dim)
          end select
          
         
       end do
    end do

  contains

    subroutine average_face_to_cc_2d(face,ng_f,cc,ng_c,face_comp,cc_comp,lo,hi,av_dim)

      integer        , intent(in   ) :: ng_f,ng_c,face_comp,cc_comp,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
      integer, intent(in) :: av_dim

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            select case(av_dim)
            case(1)
               cc(i,j,cc_comp) = 0.5d0*(face(i,j,face_comp)+face(i+1,j,face_comp))
            case(2)   
               cc(i,j,cc_comp) = 0.5d0*(face(i,j,face_comp)+face(i,j+1,face_comp))
            case default
               stop "av_dim>2 in 2d"
            end select
         end do
      end do

    end subroutine average_face_to_cc_2d

    subroutine average_face_to_cc_3d(face,ng_f,cc,ng_c,face_comp,cc_comp,lo,hi,av_dim)

      integer        , intent(in   ) :: ng_f,ng_c,face_comp,cc_comp,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
      integer, intent(in) :: av_dim

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               select case(av_dim)
               case(1)
                  cc(i,j,k,cc_comp) = 0.5d0*(face(i,j,k,face_comp)+face(i+1,j,k,face_comp))
               case(2)   
                  cc(i,j,k,cc_comp) = 0.5d0*(face(i,j,k,face_comp)+face(i,j+1,k,face_comp))
               case(3)   
                  cc(i,j,k,cc_comp) = 0.5d0*(face(i,j,k,face_comp)+face(i,j,k+1,face_comp))
               case default
                  stop "av_dim>3 in 3d"
               end select
            end do
         end do
      end do

    end subroutine average_face_to_cc_3d

  end subroutine average_face_to_cc

  subroutine average_cc_to_face(nlevs,cc,face,start_scomp,start_bccomp,num_comp,the_bc_level)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: face(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,i,dm,ng_c,ng_f,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: epx(:,:,:,:)
    real(kind=dp_t), pointer :: epy(:,:,:,:)
    real(kind=dp_t), pointer :: epz(:,:,:,:)

    dm = cc(1)%dim

    ng_c = cc(1)%ng
    ng_f = face(1,1)%ng

    if (ng_f .ge. ng_c) then
       call bl_error("average_cc_to_face requires ng_f < ng_c")
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp  => dataptr(cc(n), i)
          epx => dataptr(face(n,1), i)
          epy => dataptr(face(n,2), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          select case (dm)
          case (2)
             do scomp=start_scomp,start_scomp+num_comp-1
                bccomp = start_bccomp + scomp - start_scomp
                call average_cc_to_face_2d(cp(:,:,1,scomp),ng_c, &
                                           epx(:,:,1,scomp),epy(:,:,1,scomp),ng_f, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end do
          case (3)
             do scomp=start_scomp,start_scomp+num_comp-1
                bccomp = start_bccomp + scomp - start_scomp
                epz => dataptr(face(n,3), i)
                call average_cc_to_face_3d(cp(:,:,:,scomp),ng_c, &
                                           epx(:,:,:,scomp),epy(:,:,:,scomp),epz(:,:,:,scomp),ng_f, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end do
          end select
       end do
    end do

  end subroutine average_cc_to_face
    
  subroutine average_cc_to_face_2d(cc,ng_c,facex,facey,ng_f,lo,hi,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_f
    real(kind=dp_t), intent(in   ) ::    cc(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(inout) :: facex(lo(1)-ng_f:,lo(2)-ng_f:)
    real(kind=dp_t), intent(inout) :: facey(lo(1)-ng_f:,lo(2)-ng_f:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j

    ! x-faces
    do j=lo(2)-ng_f,hi(2)+ng_f
       do i=lo(1)-ng_f,hi(1)+ng_f+1
          facex(i,j) = 0.5d0*(cc(i,j)+cc(i-1,j))
       end do
    end do

    ! overwrite x-boundary faces
    ! value in ghost cells represents boundary value
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do i=lo(1)-ng_f,lo(1)
          facex(i,lo(2)-ng_f:hi(2)+ng_f) = cc(lo(1)-1,lo(2)-ng_f:hi(2)+ng_f)
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do i=hi(1)+1,hi(1)+ng_f+1
          facex(i,lo(2)-ng_f:hi(2)+ng_f) = cc(hi(1)+1,lo(2)-ng_f:hi(2)+ng_f)
       end do
    end if

    ! y-faces
    do j=lo(2)-ng_f,hi(2)+ng_f+1
       do i=lo(1)-ng_f,hi(1)+ng_f
          facey(i,j) = 0.5d0*(cc(i,j)+cc(i,j-1))
       end do
    end do

    ! overwrite y-boundary faces
    ! value in ghost cells represents boundary value
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_f,lo(2)
          facey(lo(1)-ng_f:hi(1)+ng_f,j) = cc(lo(1)-ng_f:hi(1)+ng_f,lo(2)-1)
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do j=hi(2)+1,hi(2)+ng_f+1
          facey(lo(1)-ng_f:hi(1)+ng_f,j) = cc(lo(1)-ng_f:hi(1)+ng_f,hi(2)+1)
       end do
    end if

  end subroutine average_cc_to_face_2d

  subroutine average_cc_to_face_3d(cc,ng_c,facex,facey,facez,ng_f,lo,hi,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_f
    real(kind=dp_t), intent(in   ) ::    cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: facex(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: facey(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    real(kind=dp_t), intent(inout) :: facez(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j,k

    ! x-faces
    do k=lo(3)-ng_f,hi(3)+ng_f
       do j=lo(2)-ng_f,hi(2)+ng_f
          do i=lo(1)-ng_f,hi(1)+ng_f+1
             facex(i,j,k) = 0.5d0*(cc(i,j,k)+cc(i-1,j,k))
          end do
       end do
    end do

    ! y-faces
    do k=lo(3)-ng_f,hi(3)+ng_f
       do j=lo(2)-ng_f,hi(2)+ng_f+1
          do i=lo(1)-ng_f,hi(1)+ng_f
             facey(i,j,k) = 0.5d0*(cc(i,j,k)+cc(i,j-1,k))
          end do
       end do
    end do

    ! z-faces
    do k=lo(3)-ng_f,hi(3)+ng_f+1
       do j=lo(2)-ng_f,hi(2)+ng_f
          do i=lo(1)-ng_f,hi(1)+ng_f
             facez(i,j,k) = 0.5d0*(cc(i,j,k)+cc(i,j,k-1))
          end do
       end do
    end do

    ! Boundary Conditions
    ! Note: At physical boundaries, the value in ghost cells represents the boundary value

    ! overwrite x-boundary faces
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do i=lo(1)-ng_f,lo(1)
          facex(i,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f) = &
               cc(lo(1)-1,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f)
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do i=hi(1)+1,hi(1)+ng_f+1
          facex(i,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f) = &
               cc(hi(1)+1,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f)
       end do
    end if

    ! overwrite y-boundary faces
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_f,lo(2)
          facey(lo(1)-ng_f:hi(1)+ng_f,j,lo(3)-ng_f:hi(3)+ng_f) = &
               cc(lo(1)-ng_f:hi(1)+ng_f,lo(2)-1,lo(3)-ng_f:hi(3)+ng_f)
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do j=hi(2)+1,hi(2)+ng_f+1
          facey(lo(1)-ng_f:hi(1)+ng_f,j,lo(3)-ng_f:hi(3)+ng_f) = &
               cc(lo(1)-ng_f:hi(1)+ng_f,hi(2)+1,lo(3)-ng_f:hi(3)+ng_f)
       end do
    end if

    ! overwrite z-boundary faces
    if (adv_bc(3,1) .eq. FOEXTRAP .or. adv_bc(3,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_f,lo(3)
          facez(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,k) = &
               cc(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,lo(3)-1)
       end do
    end if
    if (adv_bc(3,2) .eq. FOEXTRAP .or. adv_bc(3,2) .eq. EXT_DIR) then
       do k=hi(3)+1,hi(3)+ng_f+1
          facez(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,k) = &
               cc(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,hi(3)+1)
       end do
    end if

  end subroutine average_cc_to_face_3d

  subroutine average_cc_to_node(nlevs,cc,node,start_scomp,start_bccomp,num_comp,the_bc_level)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: node(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,i,dm,ng_c,ng_n,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    dm = cc(1)%dim

    ng_c = cc(1)%ng
    ng_n = node(1)%ng

    if (ng_n .ge. ng_c) then
       call bl_error("average_cc_to_node requires ng_n < ng_c")
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp => dataptr(cc(n), i)
          np => dataptr(node(n), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          select case (dm)
          case (2)
             do scomp=start_scomp,start_scomp+num_comp-1
                bccomp = start_bccomp + scomp - start_scomp
                call average_cc_to_node_2d(cp(:,:,1,scomp),ng_c, &
                                           np(:,:,1,scomp),ng_n, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
          end do
             case (3)
                call bl_error("should not have to call average_cc_to_node_3d")
                do scomp=start_scomp,start_scomp+num_comp-1
                   bccomp = start_bccomp + scomp - start_scomp
                   call average_cc_to_node_3d(cp(:,:,:,scomp),ng_c, &
                                              np(:,:,:,scomp),ng_n, &
                                              lo,hi, &
                                              the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
                end do
             end select
       end do
    end do

  end subroutine average_cc_to_node
 
  subroutine average_cc_to_node_2d(cc,ng_c,node,ng_n,lo,hi,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n
    real(kind=dp_t), intent(in   ) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(inout) :: node(lo(1)-ng_n:,lo(2)-ng_n:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j

    ! average to nodes
    do j=lo(2)-ng_n,hi(2)+1+ng_n
       do i=lo(1)-ng_n,hi(1)+1+ng_n
          node(i,j) = 0.25d0*(cc(i,j)+cc(i-1,j)+cc(i,j-1)+cc(i-1,j-1))
       end do
    end do

    ! Boundary Conditions
    ! Note: At physical boundaries, the value in ghost cells represents the boundary value

    ! overwrite x-boundary nodes
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_n,hi(2)+1+ng_n
          do i=lo(1)-ng_n,lo(1)
             node(i,j) = 0.5d0*(cc(lo(1)-1,j)+cc(lo(1)-1,j-1))
          end do
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do j=lo(2)-ng_n,hi(2)+1+ng_n
          do i=hi(1)+1,hi(1)+1+ng_n
             node(i,j) = 0.5d0*(cc(hi(1)+1,j)+cc(hi(1)+1,j-1))
          end do
       end do
    end if

    ! overwrite y-boundary nodes
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_n,lo(2)
          do i=lo(1)-ng_n,hi(1)+1+ng_n
             node(i,j) = 0.5d0*(cc(i,lo(2)-1)+cc(i-1,lo(2)-1))
          end do
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do j=hi(2)+1,hi(2)+1+ng_n
          do i=lo(1)-ng_n,hi(1)+1+ng_n
             node(i,j) = 0.5d0*(cc(i,hi(2)+1)+cc(i-1,hi(2)+1))
          end do
       end do
    end if

  end subroutine average_cc_to_node_2d

  subroutine average_cc_to_node_3d(cc,ng_c,node,ng_n,lo,hi,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n
    real(kind=dp_t), intent(in   ) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: node(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j,k

    ! average to nodes
    do k=lo(3)-ng_n,hi(3)+1+ng_n
       do j=lo(2)-ng_n,hi(2)+1+ng_n
          do i=lo(1)-ng_n,hi(1)+1+ng_n
             node(i,j,k) = 0.125d0*( cc(i,j,k)+cc(i-1,j,k)+cc(i,j-1,k)+cc(i,j,k-1) &
                                    +cc(i-1,j-1,k)+cc(i-1,j,k-1) &
                                    +cc(i,j-1,k-1)+cc(i-1,j-1,k-1) )
          end do
       end do
    end do

    ! Boundary Conditions
    ! Note: At physical boundaries, the value in ghost cells represents the boundary value

    ! overwrite x-boundary nodes
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_n,hi(3)+ng_n
          do j=lo(2)-ng_n,hi(2)+ng_n
             do i=lo(1)-ng_n,lo(1)
                node(i,j,k) = 0.25d0*( cc(lo(1)-1,j,k)+cc(lo(1)-1,j-1,k) &
                                      +cc(lo(1)-1,j,k-1)+cc(lo(1)-1,j-1,k-1))
             end do
          end do
       end do
    end if
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then
       do k=lo(3)-ng_n,hi(3)+ng_n
          do j=lo(2)-ng_n,hi(2)+ng_n
             do i=hi(1)+1,hi(1)+1+ng_n
                node(i,j,k) = 0.25d0*( cc(hi(1)+1,j,k)+cc(hi(1)+1,j-1,k) &
                                      +cc(hi(1)+1,j,k-1)+cc(hi(1)+1,j-1,k-1))
             end do
          end do
       end do
    end if

    ! overwrite y-boundary nodes
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_n,hi(3)+ng_n
          do j=lo(2)-ng_n,lo(2)
             do i=lo(1)-ng_n,hi(1)+ng_n
                node(i,j,k) = 0.25d0*( cc(i,lo(2)-1,k)+cc(i-1,lo(2)-1,k) &
                                      +cc(i,lo(2)-1,k-1)+cc(i-1,lo(2)-1,k-1))
             end do
          end do
       end do
    end if
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then
       do k=lo(3)-ng_n,hi(3)+ng_n
          do j=hi(2)+1,hi(2)+1+ng_n
             do i=lo(1)-ng_n,hi(1)+ng_n
                node(i,j,k) = 0.25d0*( cc(i,hi(2)+1,k)+cc(i-1,hi(2)+1,k) &
                                      +cc(i,hi(2)+1,k-1)+cc(i-1,hi(2)+1,k-1))
             end do
          end do
       end do
    end if

    ! overwrite z-boundary nodes
    if (adv_bc(3,1) .eq. FOEXTRAP .or. adv_bc(3,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_n,lo(3)
          do j=lo(2)-ng_n,hi(2)+ng_n
             do i=lo(1)-ng_n,hi(1)+ng_n
                node(i,j,k) = 0.5d0*( cc(i,j,lo(3)-1)+cc(i-1,j,lo(3)-1) &
                                     +cc(i,j-1,lo(3)-1)+cc(i-1,j-1,lo(3)-1))
             end do
          end do
       end do
    end if
    if (adv_bc(3,2) .eq. FOEXTRAP .or. adv_bc(3,2) .eq. EXT_DIR) then
       do k=hi(3)+1,hi(3)+1+ng_n
          do j=lo(2)-ng_n,hi(2)+ng_n
             do i=lo(1)-ng_n,hi(1)+ng_n
                node(i,j,k) = 0.5d0*( cc(i,j,hi(3)+1)+cc(i-1,j,hi(3)+1) &
                                     +cc(i,j-1,hi(3)+1)+cc(i-1,j-1,hi(3)+1))
             end do
          end do
       end do
    end if

  end subroutine average_cc_to_node_3d

  subroutine average_cc_to_edge(nlevs,cc,edge,start_scomp,start_bccomp, &
                                num_comp,the_bc_level)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: edge(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,i,dm,ng_c,ng_e,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:)
    real(kind=dp_t), pointer :: ep2(:,:,:,:)
    real(kind=dp_t), pointer :: ep3(:,:,:,:)

    dm = cc(1)%dim

    ng_c = cc(1)%ng
    ng_e = edge(1,1)%ng

    if (ng_e .ge. ng_c) then
       call bl_error("average_cc_to_edge requires ng_e < ng_c")
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp  => dataptr(cc(n), i)
          ep1 => dataptr(edge(n,1), i)
          ep2 => dataptr(edge(n,2), i)
          ep3 => dataptr(edge(n,3), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          select case (dm)
          case (2)
             call bl_error("no such thing as average_cc_to_edge in 2d")
          case (3)
             do scomp=start_scomp,start_scomp+num_comp-1
                bccomp = start_bccomp + scomp - start_scomp
                call average_cc_to_edge_3d(cp(:,:,:,scomp),ng_c, &
                                           ep1(:,:,:,scomp),ep2(:,:,:,scomp), &
                                           ep3(:,:,:,scomp),ng_e, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end do
          end select
       end do
    end do

  end subroutine average_cc_to_edge

  subroutine average_cc_to_edge_3d(cc,ng_c,edge_xy,edge_xz,edge_yz,ng_e,lo,hi,adv_bc)

    use bc_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_e
    real(kind=dp_t), intent(in   ) ::      cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(inout) :: edge_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edge_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edge_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    integer        , intent(in   ) :: adv_bc(:,:)

    ! local
    integer :: i,j,k

    ! xy edges
    do k=lo(3)-ng_e,hi(3)+ng_e
       do j=lo(2)-ng_e,hi(2)+1+ng_e
          do i=lo(1)-ng_e,hi(1)+1+ng_e
             edge_xy(i,j,k) = 0.25d0*(cc(i,j,k)+cc(i-1,j,k)+cc(i,j-1,k)+cc(i-1,j-1,k))
          end do
       end do
    end do

    ! xz edges
    do k=lo(3)-ng_e,hi(3)+1+ng_e
       do j=lo(2)-ng_e,hi(2)+ng_e
          do i=lo(1)-ng_e,hi(1)+1+ng_e
             edge_xz(i,j,k) = 0.25d0*(cc(i,j,k)+cc(i-1,j,k)+cc(i,j,k-1)+cc(i-1,j,k-1))
          end do
       end do
    end do

    ! yz edges
    do k=lo(3)-ng_e,hi(3)+1+ng_e
       do j=lo(2)-ng_e,hi(2)+1+ng_e
          do i=lo(1)-ng_e,hi(1)+ng_e
             edge_yz(i,j,k) = 0.25d0*(cc(i,j,k)+cc(i,j-1,k)+cc(i,j,k-1)+cc(i,j-1,k-1))
          end do
       end do
    end do

    ! Boundary Conditions
    ! Note: At physical boundaries, the value in ghost cells represents the boundary value

    ! lo-x boundary
    if (adv_bc(1,1) .eq. FOEXTRAP .or. adv_bc(1,1) .eq. EXT_DIR) then

       ! overwrite edge_xy
       do k=lo(3)-ng_e,hi(3)+ng_e
          do j=lo(2)-ng_e,hi(2)+1+ng_e
             do i=lo(1)-ng_e,lo(1)
                edge_xy(i,j,k) = 0.5d0*(cc(lo(1)-1,j,k)+cc(lo(1)-1,j-1,k))
             end do
          end do
       end do

       ! overwrite edge_xz
       do k=lo(3)-ng_e,hi(3)+1+ng_e
          do j=lo(2)-ng_e,hi(2)+ng_e
             do i=lo(1)-ng_e,lo(1)
                edge_xz(i,j,k) = 0.5d0*(cc(lo(1)-1,j,k)+cc(lo(1)-1,j,k-1))
             end do
          end do
       end do

    end if

    ! hi-x boundary
    if (adv_bc(1,2) .eq. FOEXTRAP .or. adv_bc(1,2) .eq. EXT_DIR) then

       ! overwrite edge_xy
       do k=lo(3)-ng_e,hi(3)+ng_e
          do j=lo(2)-ng_e,hi(2)+1+ng_e
             do i=hi(1)+1,hi(1)+1+ng_e
                edge_xy(i,j,k) = 0.5d0*(cc(hi(1)+1,j,k)+cc(hi(1)+1,j-1,k))
             end do
          end do
       end do

       ! overwrite edge_xz
       do k=lo(3)-ng_e,hi(3)+1+ng_e
          do j=lo(2)-ng_e,hi(2)+ng_e
             do i=hi(1)+1,hi(1)+1+ng_e
                edge_xz(i,j,k) = 0.5d0*(cc(hi(1)+1,j,k)+cc(hi(1)+1,j,k-1))
             end do
          end do
       end do

    end if

    ! lo-y boundary
    if (adv_bc(2,1) .eq. FOEXTRAP .or. adv_bc(2,1) .eq. EXT_DIR) then

       ! overwrite edge_xy
       do k=lo(3)-ng_e,hi(3)+ng_e
          do j=lo(2)-ng_e,lo(2)
             do i=lo(1)-ng_e,hi(1)+1+ng_e
                edge_xy(i,j,k) = 0.5d0*(cc(i,lo(2)-1,k)+cc(i-1,lo(2)-1,k))
             end do
          end do
       end do

       ! overwrite edge_yz
       do k=lo(3)-ng_e,hi(3)+1+ng_e
          do j=lo(2)-ng_e,lo(2)
             do i=lo(1)-ng_e,hi(1)+ng_e
                edge_yz(i,j,k) = 0.5d0*(cc(i,lo(2)-1,k)+cc(i,lo(2)-1,k-1))
             end do
          end do
       end do

    end if

    ! hi-y boundary
    if (adv_bc(2,2) .eq. FOEXTRAP .or. adv_bc(2,2) .eq. EXT_DIR) then

       ! overwrite edge_xy
       do k=lo(3)-ng_e,hi(3)+ng_e
          do j=hi(2)+1,hi(2)+1+ng_e
             do i=lo(1)-ng_e,hi(1)+1+ng_e
                edge_xy(i,j,k) = 0.5d0*(cc(i,hi(2)+1,k)+cc(i-1,hi(2)+1,k))
             end do
          end do
       end do

       ! overwrite edge_yz
       do k=lo(3)-ng_e,hi(3)+1+ng_e
          do j=hi(2)+1,hi(2)+1+ng_e
             do i=lo(1)-ng_e,hi(1)+ng_e
                edge_yz(i,j,k) = 0.5d0*(cc(i,hi(2)+1,k)+cc(i,hi(2)+1,k-1))
             end do
          end do
       end do

    end if

    ! lo-z boundary
    if (adv_bc(3,1) .eq. FOEXTRAP .or. adv_bc(3,1) .eq. EXT_DIR) then

       ! overwrite edge_xz
       do k=lo(3)-ng_e,lo(3)
          do j=lo(2)-ng_e,hi(2)+ng_e
             do i=lo(1)-ng_e,hi(1)+1+ng_e
                edge_xz(i,j,k) = 0.5d0*(cc(i,j,lo(3)-1)+cc(i-1,j,lo(3)-1))
             end do
          end do
       end do

       ! overwrite edge_yz
       do k=lo(3)-ng_e,lo(3)
          do j=lo(2)-ng_e,hi(2)+1+ng_e
             do i=lo(1)-ng_e,hi(1)+ng_e
                edge_yz(i,j,k) = 0.5d0*(cc(i,j,lo(3)-1)+cc(i,j-1,lo(3)-1))
             end do
          end do
       end do

    end if

    ! hi-z boundary
    if (adv_bc(3,2) .eq. FOEXTRAP .or. adv_bc(3,2) .eq. EXT_DIR) then

       ! overwrite edge_xz
       do k=hi(3)+1,hi(3)+1+ng_e
          do j=lo(2)-ng_e,hi(2)+ng_e
             do i=lo(1)-ng_e,hi(1)+1+ng_e
                edge_xz(i,j,k) = 0.5d0*(cc(i,j,hi(3)+1)+cc(i-1,j,hi(3)+1))
             end do
          end do
       end do

       ! overwrite edge_yz
       do k=hi(3)+1,hi(3)+1+ng_e
          do j=lo(2)-ng_e,hi(2)+1+ng_e
             do i=lo(1)-ng_e,hi(1)+ng_e
                edge_yz(i,j,k) = 0.5d0*(cc(i,j,hi(3)+1)+cc(i,j-1,hi(3)+1))
             end do
          end do
       end do

    end if

  end subroutine average_cc_to_edge_3d

end module convert_module
