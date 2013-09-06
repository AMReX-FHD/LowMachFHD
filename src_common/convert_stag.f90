module convert_stag_module

  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: average_face_to_cc, average_cc_to_face, average_cc_to_node, &
       average_cc_to_edge, shift_face_to_cc
  
contains

  subroutine average_face_to_cc(mla,face,start_face_comp,cc,start_cc_comp,num_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: face(:)
    type(multifab) , intent(inout) ::   cc(:)
    integer        , intent(in   ) :: start_face_comp,start_cc_comp,num_comp

    ! local
    integer :: n,i,dm,nlevs,ng_f,ng_c
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: av_dim ! Along which dimension to do the average
    integer :: face_comp,cc_comp

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
          do face_comp=start_face_comp,start_face_comp+num_comp-1
             cc_comp = start_cc_comp + (face_comp-start_face_comp)
             select case (dm)
             case (2)
                call average_face_to_cc_2d(ep(:,:,1,face_comp), ng_f, cp(:,:,1,cc_comp), ng_c, &
                                           lo, hi, av_dim)
             case (3)
                call average_face_to_cc_3d(ep(:,:,:,face_comp), ng_f, cp(:,:,:,cc_comp), ng_c, &
                                           lo, hi, av_dim)
             end select
          end do
       end do
    end do

  contains

    subroutine average_face_to_cc_2d(face,ng_f,cc,ng_c,lo,hi,av_dim)

      integer        , intent(in   ) :: ng_f,ng_c,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:)
      integer, intent(in) :: av_dim

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            select case(av_dim)
            case(1)
               cc(i,j) = 0.5d0*(face(i,j)+face(i+1,j))
            case(2)   
               cc(i,j) = 0.5d0*(face(i,j)+face(i,j+1))
            case default
               stop "av_dim>2 in 2d"
            end select
         end do
      end do

    end subroutine average_face_to_cc_2d

    subroutine average_face_to_cc_3d(face,ng_f,cc,ng_c,lo,hi,av_dim)

      integer        , intent(in   ) :: ng_f,ng_c,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
      integer, intent(in) :: av_dim

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               select case(av_dim)
               case(1)
                  cc(i,j,k) = 0.5d0*(face(i,j,k)+face(i+1,j,k))
               case(2)   
                  cc(i,j,k) = 0.5d0*(face(i,j,k)+face(i,j+1,k))
               case(3)   
                  cc(i,j,k) = 0.5d0*(face(i,j,k)+face(i,j,k+1))
               case default
                  stop "av_dim>3 in 3d"
               end select
            end do
         end do
      end do

    end subroutine average_face_to_cc_3d

  end subroutine average_face_to_cc

  subroutine average_cc_to_face(nlevs,cc,face,start_scomp,start_bccomp,num_comp, &
                                the_bc_level,increment_bccomp_in)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: face(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical,  intent(in), optional :: increment_bccomp_in

    ! local
    integer :: n,i,dm,ng_c,ng_f,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)
    logical :: increment_bccomp

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

    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp  => dataptr(cc(n), i)
          epx => dataptr(face(n,1), i)
          epy => dataptr(face(n,2), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          do scomp=start_scomp,start_scomp+num_comp-1
             if (increment_bccomp) then
                bccomp = start_bccomp + scomp - start_scomp
             else
                bccomp = start_bccomp
             end if
             select case (dm)
             case (2)
                call average_cc_to_face_2d(cp(:,:,1,scomp),ng_c, &
                                           epx(:,:,1,scomp),epy(:,:,1,scomp),ng_f, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             case (3)
                epz => dataptr(face(n,3), i)
                call average_cc_to_face_3d(cp(:,:,:,scomp),ng_c, &
                                           epx(:,:,:,scomp),epy(:,:,:,scomp),epz(:,:,:,scomp),ng_f, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end select
          end do
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
          if(.false.) then
            if(i.eq.9 .and. j.eq.11) then
               print*, facex(i,j)
            endif
          endif
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

  subroutine average_cc_to_node(nlevs,cc,node,start_scomp,start_bccomp,num_comp, &
                                the_bc_level,increment_bccomp_in)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: node(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical,  intent(in), optional :: increment_bccomp_in

    ! local
    integer :: n,i,dm,ng_c,ng_n,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)
    logical :: increment_bccomp

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    dm = cc(1)%dim

    ng_c = cc(1)%ng
    ng_n = node(1)%ng

    if (ng_n .ge. ng_c) then
       call bl_error("average_cc_to_node requires ng_n < ng_c")
    end if
    
    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp => dataptr(cc(n), i)
          np => dataptr(node(n), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          do scomp=start_scomp,start_scomp+num_comp-1
             if (increment_bccomp) then
                bccomp = start_bccomp + scomp - start_scomp
             else
                bccomp = start_bccomp
             end if
             select case (dm)
             case (2)
                call average_cc_to_node_2d(cp(:,:,1,scomp),ng_c, &
                                           np(:,:,1,scomp),ng_n, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             case (3)
                call average_cc_to_node_3d(cp(:,:,:,scomp),ng_c, &
                                           np(:,:,:,scomp),ng_n, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end select
          end do
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
                                num_comp,the_bc_level,increment_bccomp_in)

    integer        , intent(in   ) :: nlevs,start_scomp,start_bccomp,num_comp
    type(multifab) , intent(in   ) :: cc(:)
    type(multifab) , intent(inout) :: edge(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical,  intent(in), optional :: increment_bccomp_in

    ! local
    integer :: n,i,dm,ng_c,ng_e,scomp,bccomp
    integer :: lo(cc(1)%dim),hi(cc(1)%dim)
    logical :: increment_bccomp

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

    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if

    do n=1,nlevs
       do i=1,nfabs(cc(n))
          cp  => dataptr(cc(n), i)
          ep1 => dataptr(edge(n,1), i)
          ep2 => dataptr(edge(n,2), i)
          ep3 => dataptr(edge(n,3), i)
          lo = lwb(get_box(cc(n), i))
          hi = upb(get_box(cc(n), i))
          do scomp=start_scomp,start_scomp+num_comp-1
             if (increment_bccomp) then
                bccomp = start_bccomp + scomp - start_scomp
             else
                bccomp = start_bccomp
             end if
             select case (dm)
             case (2)
                call bl_error("no such thing as average_cc_to_edge in 2d")
             case (3)
                call average_cc_to_edge_3d(cp(:,:,:,scomp),ng_c, &
                                           ep1(:,:,:,scomp),ep2(:,:,:,scomp), &
                                           ep3(:,:,:,scomp),ng_e, &
                                           lo,hi, &
                                           the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp))
             end select
          end do
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

  subroutine shift_face_to_cc(mla,face,start_face_comp,cc,start_cc_comp,num_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: face(:)
    type(multifab) , intent(inout) ::   cc(:)
    integer        , intent(in   ) :: start_face_comp,start_cc_comp,num_comp

    ! local
    integer :: n,i,dm,nlevs,ng_f,ng_c
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: face_comp,cc_comp

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
          do face_comp=start_face_comp,start_face_comp+num_comp-1
             cc_comp = start_cc_comp + (face_comp-start_face_comp)
             select case (dm)
             case (2)
                call shift_face_to_cc_2d(ep(:,:,1,face_comp), ng_f, cp(:,:,1,cc_comp), ng_c, &
                                         lo, hi)
             case (3)
                call shift_face_to_cc_3d(ep(:,:,:,face_comp), ng_f, cp(:,:,:,cc_comp), ng_c, &
                                         lo, hi)
             end select
          end do
       end do
    end do

  contains
    
    subroutine shift_face_to_cc_2d(face,ng_f,cc,ng_c,lo,hi)

      integer        , intent(in   ) :: ng_f,ng_c,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:)

      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            cc(i,j) = face(i,j)
         end do
      end do

    end subroutine shift_face_to_cc_2d

    subroutine shift_face_to_cc_3d(face,ng_f,cc,ng_c,lo,hi)

      integer        , intent(in   ) :: ng_f,ng_c,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: face(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(  out) ::   cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)

      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               cc(i,j,k) = face(i,j,k)
            end do
         end do
      end do

    end subroutine shift_face_to_cc_3d

  end subroutine shift_face_to_cc

end module convert_stag_module
