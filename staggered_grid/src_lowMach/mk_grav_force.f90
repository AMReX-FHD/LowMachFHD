module mk_grav_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use zero_edgeval_module
  use probin_common_module, only: grav, rhobar, nspecies, rho0
  use eos_check_module, only: compute_rho_eos

  implicit none

  private

  public :: mk_grav_force, mk_grav_force_bousq

contains

  subroutine mk_grav_force(mla,m_force,increment,rhotot_fc_old,rhotot_fc_new,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_force(:,:)
    logical        , intent(in   ) :: increment
    type(multifab) , intent(in   ) :: rhotot_fc_old(:,:)
    type(multifab) , intent(in   ) :: rhotot_fc_new(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer :: sox(:,:,:,:)
    real(kind=dp_t), pointer :: soy(:,:,:,:)
    real(kind=dp_t), pointer :: soz(:,:,:,:)
    real(kind=dp_t), pointer :: snx(:,:,:,:)
    real(kind=dp_t), pointer :: sny(:,:,:,:)
    real(kind=dp_t), pointer :: snz(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"mk_grav_force")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s  = rhotot_fc_old(1,1)%ng
    ng_u  = m_force(1,1)%ng

    do n=1,nlevs
       do i=1, nfabs(m_force(n,1))
          fxp => dataptr(m_force(n,1), i)
          fyp => dataptr(m_force(n,2), i)
          sox => dataptr(rhotot_fc_old(n,1), i)
          soy => dataptr(rhotot_fc_old(n,2), i)
          snx => dataptr(rhotot_fc_new(n,1), i)
          sny => dataptr(rhotot_fc_new(n,2), i)
          lo = lwb(get_box(m_force(n,1), i))
          hi = upb(get_box(m_force(n,1), i))
          select case (dm)
          case (2)
             call mk_grav_force_2d(fxp(:,:,1,1), fyp(:,:,1,1), ng_u, &
                                   sox(:,:,1,1), soy(:,:,1,1), &
                                   snx(:,:,1,1), sny(:,:,1,1), ng_s, lo, hi, &
                                   increment)
          case (3)
             fzp => dataptr(m_force(n,3), i)
             soz => dataptr(rhotot_fc_old(n,3), i)
             snz => dataptr(rhotot_fc_new(n,3), i)
             call mk_grav_force_3d(fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_u, &
                                   sox(:,:,:,1), soy(:,:,:,1), soz(:,:,:,1), &
                                   snx(:,:,:,1), sny(:,:,:,1), snz(:,:,:,1), ng_s, lo, hi, &
                                   increment)
          end select
       end do

       ! zero wall boundary values
       call zero_edgeval_physical(m_force(n,:),1,1,the_bc_tower%bc_tower_array(n))

    enddo

    call destroy(bpt)

  end subroutine mk_grav_force

  subroutine mk_grav_force_2d(m_forcex,m_forcey,ng_u,rhotot_oldx,rhotot_oldy,rhotot_newx,rhotot_newy, &
                              ng_s,lo,hi,increment)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) ::    m_forcex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) ::    m_forcey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rhotot_oldx(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_oldy(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_newx(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_newy(lo(1)-ng_s:,lo(2)-ng_s:)
    logical        , intent(in   ) :: increment

    ! local
    integer i,j

    if (increment) then

       do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          m_forcex(i,j) = m_forcex(i,j) + 0.5d0*grav(1)*(rhotot_oldx(i,j)+rhotot_newx(i,j))
       end do
       end do

       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          m_forcey(i,j) = m_forcey(i,j) + 0.5d0*grav(2)*(rhotot_oldy(i,j)+rhotot_newy(i,j))
       end do
       end do

    else

       do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          m_forcex(i,j) = 0.5d0*grav(1)*(rhotot_oldx(i,j)+rhotot_newx(i,j))
       end do
       end do

       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          m_forcey(i,j) = 0.5d0*grav(2)*(rhotot_oldy(i,j)+rhotot_newy(i,j))
       end do
       end do

    end if

  end subroutine mk_grav_force_2d

  subroutine mk_grav_force_3d(m_forcex,m_forcey,m_forcez,ng_u, &
                              rhotot_oldx,rhotot_oldy,rhotot_oldz, &
                              rhotot_newx,rhotot_newy,rhotot_newz,ng_s,lo,hi,increment)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) ::    m_forcex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::    m_forcey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::    m_forcez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rhotot_oldx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_oldy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_oldz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_newx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_newy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: rhotot_newz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    logical        , intent(in   ) :: increment

    ! local
    integer i,j,k

    if (increment) then

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          m_forcex(i,j,k) = m_forcex(i,j,k) + 0.5d0*grav(1)*(rhotot_oldx(i,j,k)+rhotot_newx(i,j,k))
       end do
       end do
       end do

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          m_forcey(i,j,k) = m_forcey(i,j,k) + 0.5d0*grav(2)*(rhotot_oldy(i,j,k)+rhotot_newy(i,j,k))
       end do
       end do
       end do

       do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          m_forcez(i,j,k) = m_forcez(i,j,k) + 0.5d0*grav(3)*(rhotot_oldz(i,j,k)+rhotot_newz(i,j,k))
       end do
       end do
       end do

    else

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          m_forcex(i,j,k) = 0.5d0*grav(1)*(rhotot_oldx(i,j,k)+rhotot_newx(i,j,k))
       end do
       end do
       end do

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          m_forcey(i,j,k) = 0.5d0*grav(2)*(rhotot_oldy(i,j,k)+rhotot_newy(i,j,k))
       end do
       end do
       end do

       do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          m_forcez(i,j,k) = 0.5d0*grav(3)*(rhotot_oldz(i,j,k)+rhotot_newz(i,j,k))
       end do
       end do
       end do

    end if

  end subroutine mk_grav_force_3d

  subroutine mk_grav_force_bousq(mla,m_force,increment,rho_fc,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_force(:,:)
    logical        , intent(in   ) :: increment
    type(multifab) , intent(in   ) :: rho_fc(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer :: sox(:,:,:,:)
    real(kind=dp_t), pointer :: soy(:,:,:,:)
    real(kind=dp_t), pointer :: soz(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"mk_grav_force_bousq")

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_s  = rho_fc(1,1)%ng
    ng_u  = m_force(1,1)%ng

    do n=1,nlevs
       do i=1, nfabs(m_force(n,1))
          fxp => dataptr(m_force(n,1), i)
          fyp => dataptr(m_force(n,2), i)
          sox => dataptr(rho_fc(n,1), i)
          soy => dataptr(rho_fc(n,2), i)
          lo = lwb(get_box(m_force(n,1), i))
          hi = upb(get_box(m_force(n,1), i))
          select case (dm)
          case (2)
             call mk_grav_force_bousq_2d(fxp(:,:,1,1), fyp(:,:,1,1), ng_u, &
                                         sox(:,:,1,:), soy(:,:,1,:), &
                                         ng_s, lo, hi, &
                                         increment)
          case (3)
             fzp => dataptr(m_force(n,3), i)
             soz => dataptr(rho_fc(n,3), i)
             call mk_grav_force_bousq_3d(fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_u, &
                                         sox(:,:,:,:), soy(:,:,:,:), soz(:,:,:,:), &
                                         ng_s, lo, hi, &
                                         increment)
          end select
       end do

       ! zero wall boundary values
       call zero_edgeval_physical(m_force(n,:),1,1,the_bc_tower%bc_tower_array(n))

    enddo

    call destroy(bpt)

  end subroutine mk_grav_force_bousq

  subroutine mk_grav_force_bousq_2d(m_forcex,m_forcey,ng_u,rhox,rhoy, &
                                    ng_s,lo,hi,increment)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) ::     rhox(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::     rhoy(lo(1)-ng_s:,lo(2)-ng_s:,:)
    logical        , intent(in   ) :: increment

    ! local
    integer i,j,comp
    real(kind=dp_t) :: rhotot

    if(.not.increment) then
       m_forcex=0.0d0
       m_forcey=0.0d0
    end if
    
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1

       call compute_rho_eos(rhox(i,j,1:nspecies), rhotot) 
       m_forcex(i,j) = m_forcex(i,j) + grav(1)*(rhotot-rho0)

    end do
    end do

    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
       
       call compute_rho_eos(rhoy(i,j,1:nspecies), rhotot) 
       m_forcey(i,j) = m_forcey(i,j) + grav(2)*(rhotot-rho0)

    end do
    end do


  end subroutine mk_grav_force_bousq_2d

  subroutine mk_grav_force_bousq_3d(m_forcex,m_forcey,m_forcez,ng_u, &
                                    rhox,rhoy,rhoz,ng_s,lo,hi,increment)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) ::     rhox(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::     rhoy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::     rhoz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    logical        , intent(in   ) :: increment

    ! local
    integer i,j,k,comp
    real(kind=dp_t) :: rhotot

    if(.not.increment) then
       m_forcex=0.0d0
       m_forcey=0.0d0
       m_forcez=0.0d0
    end if


    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1

       call compute_rho_eos(rhox(i,j,k,1:nspecies), rhotot)
       m_forcex(i,j,k) = m_forcex(i,j,k) + grav(1)*(rhotot-rho0)

    end do
    end do
    end do

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)

       call compute_rho_eos(rhoy(i,j,k,1:nspecies), rhotot)
       m_forcey(i,j,k) = m_forcey(i,j,k) + grav(2)*(rhotot-rho0)

    end do
    end do
    end do

    do k=lo(3),hi(3)+1
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       call compute_rho_eos(rhoz(i,j,k,1:nspecies), rhotot)
       m_forcez(i,j,k) = m_forcez(i,j,k) + grav(3)*(rhotot-rho0)

    end do
    end do
    end do

  end subroutine mk_grav_force_bousq_3d

end module mk_grav_force_module
