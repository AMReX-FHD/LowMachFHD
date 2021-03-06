module stag_applyop_module

  use ml_layout_module
  use multifab_module

  implicit none

  private

  public :: stag_applyop, stag_applyop_2d, stag_applyop_3d

contains
  
  ! compute Lphi
  subroutine stag_applyop(mla,the_bc_tower,phi_fc,Lphi_fc,alpha_cc,beta_cc, &
                          gamma_cc,dx,comp_in)

    use define_bc_module
    use convert_module

    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: phi_fc(:,:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:,:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_cc(:)   ! face-centered
    type(multifab) , intent(in   ) :: beta_cc(:)    ! cell-centered
    type(multifab) , intent(in   ) :: gamma_cc(:)   ! cell-centered
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ), optional :: comp_in

    ! local
    integer :: i,n,comp,dm,nlevs

    type(multifab) :: alpha_fc(mla%nlevel,mla%dim)
    type(multifab) :: beta_nd(mla%nlevel)
    type(multifab) :: beta_ed(mla%nlevel,3)

    logical :: nodal_temp(mla%dim)

    if ( present(comp_in) ) then
       comp = comp_in
    else
       comp = -999
    end if

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(alpha_fc(n,i),mla%la(n),1,0,i)
       end do
    end do
    call average_cc_to_face(nlevs,alpha_cc,alpha_fc,1,dm+1,1,the_bc_tower%bc_tower_array)
       
    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_build_nodal(beta_nd(n),mla%la(n),1,0)
       end do
       call average_cc_to_node(nlevs,beta_cc,beta_nd,1,dm+1,1,the_bc_tower%bc_tower_array)
    else
       do n=1,nlevs
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(beta_ed(n,1),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(beta_ed(n,2),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(beta_ed(n,3),mla%la(n),1,0,nodal_temp)
       end do
       call average_cc_to_edge(nlevs,beta_cc,beta_ed,1,dm+1,1,the_bc_tower%bc_tower_array)
    end if

    if (dm .eq. 2) then
       do n=1,nlevs
          call stag_applyop_2d(mla%la(n),the_bc_tower%bc_tower_array(n), &
                               phi_fc(n,:),Lphi_fc(n,:),alpha_fc(n,:),beta_cc(n), &
                               beta_nd(n),gamma_cc(n),dx(n,:),comp)
       end do
    else if (dm .eq. 3) then
       do n=1,nlevs
          call stag_applyop_3d(mla%la(n),the_bc_tower%bc_tower_array(n), &
                               phi_fc(n,:),Lphi_fc(n,:),alpha_fc(n,:),beta_cc(n),beta_ed(n,:), &
                               gamma_cc(n),dx(n,:),comp)
       end do
    end if
    
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(alpha_fc(n,i))
       end do
       if (dm .eq. 2) then
          call multifab_destroy(beta_nd(n))
       else if (dm .eq. 3) then
          do i=1,3
             call multifab_destroy(beta_ed(n,i))
          end do
       end if
    end do

  end subroutine stag_applyop

  ! compute Lphi
  subroutine stag_applyop_2d(la,the_bc_level,phi_fc,Lphi_fc,alpha_fc,beta_cc, &
                             beta_nd,gamma_cc,dx,comp_in)
    
    use define_bc_module
    use multifab_physbc_module
    use probin_module, only: visc_type

    type(layout)   , intent(in   ) :: la
    type(bc_level) , intent(in   ) :: the_bc_level
    type(multifab) , intent(in   ) :: phi_fc(:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_fc(:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc     ! cell-centered
    type(multifab) , intent(in   ) :: beta_nd     ! nodal
    type(multifab) , intent(in   ) :: gamma_cc    ! cell-centered
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ), optional :: comp_in
    
    ! local
    integer :: comp,i,dm,ng_p,ng_l,ng_a,ng_b,ng_n,ng_g
    integer :: lo(get_dim(la)), hi(get_dim(la))
    
    real(kind=dp_t), pointer :: ppx(:,:,:,:)
    real(kind=dp_t), pointer :: ppy(:,:,:,:)
    real(kind=dp_t), pointer :: lpx(:,:,:,:)
    real(kind=dp_t), pointer :: lpy(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:)
    real(kind=dp_t), pointer :: apy(:,:,:,:)
    real(kind=dp_t), pointer ::  bp(:,:,:,:)
    real(kind=dp_t), pointer :: bnp(:,:,:,:)
    real(kind=dp_t), pointer :: kp(:,:,:,:)

    if ( present(comp_in) .and. abs(visc_type) .ne. 1) then
       call bl_error("stag_applyop_2d: comp not needed for abs(visc_type) != 1")
    end if

    if ( (.not. present(comp_in)) .and. abs(visc_type) .eq. 1) then
       call bl_error("stag_applyop_2d: comp required if abs(visc_type) = 1")
    end if

    if (dx(1) .ne. dx(2)) then
       call bl_error("stag_applyop_2d requires the same dx in all directions")
    end if

    if ( present(comp_in) ) then
       comp = comp_in
    else
       comp = -999
    end if

    dm = get_dim(la)

    if (dm .ne. 2) then
       call bl_error("stag_applyop_2d being called with dm != 2")
    end if

    ng_p = phi_fc(1)%ng
    ng_l = Lphi_fc(1)%ng
    ng_a = alpha_fc(1)%ng
    ng_b = beta_cc%ng
    ng_n = beta_nd%ng
    ng_g = gamma_cc%ng

    do i=1,nfabs(phi_fc(1))
       ppx => dataptr(phi_fc(1), i)
       ppy => dataptr(phi_fc(2), i)
       lpx => dataptr(Lphi_fc(1), i)
       lpy => dataptr(Lphi_fc(2), i)
       apx => dataptr(alpha_fc(1), i)
       apy => dataptr(alpha_fc(2), i)
       bp  => dataptr(beta_cc, i)
       bnp => dataptr(beta_nd, i)
       kp  => dataptr(gamma_cc, i)
       lo = lwb(get_box(phi_fc(1), i))
       hi = upb(get_box(phi_fc(1), i))
       call stag_applyop_2d_work(ppx(:,:,1,1),ppy(:,:,1,1),ng_p, &
                                 lpx(:,:,1,1),lpy(:,:,1,1),ng_l, &
                                 apx(:,:,1,1),apy(:,:,1,1),ng_a, &
                                 bp(:,:,1,1),ng_b, &
                                 bnp(:,:,1,1),ng_n, &
                                 kp(:,:,1,1),ng_g, &
                                 lo,hi,dx,comp)
    end do

    do i=1,dm
       ! set Lphi on physical domain boundaries to zero
       call multifab_physbc_domainvel(Lphi_fc(i),1,i,1,the_bc_level)
    end do

  end subroutine stag_applyop_2d

  subroutine stag_applyop_2d_work(phix,phiy,ng_p,Lpx,Lpy,ng_l, &
                                  alphax,alphay,ng_a,beta,ng_b,beta_nd,ng_n, &
                                  gamma,ng_g,lo,hi,dx,comp)

    use probin_module, only: visc_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_l,ng_a,ng_b,ng_n,ng_g
    real(kind=dp_t), intent(in   ) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(inout) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_nd(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: comp

    ! local
    integer :: i,j

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds
    real(kind=dp_t) :: b,c

    dxsq = dx(1)**2
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0

    if (visc_type .eq. -1) then
       
       if (comp .eq. 1) then

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     (beta(i,j)+beta(i-1,j)+beta_nd(i,j+1)+beta_nd(i,j))/dxsq) &
                     - ( phix(i+1,j)*beta(i,j) &
                        +phix(i-1,j)*beta(i-1,j) &
                        +phix(i,j+1)*beta_nd(i,j+1) &
                        +phix(i,j-1)*beta_nd(i,j) )/dxsq

             end do
          end do

       else if (comp .eq. 2) then

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     (beta(i,j)+beta(i,j-1)+beta_nd(i+1,j)+beta_nd(i,j))/dxsq) &
                     - ( phiy(i,j+1)*beta(i,j) &
                        +phiy(i,j-1)*beta(i,j-1) &
                        +phiy(i+1,j)*beta_nd(i+1,j) &
                        +phiy(i-1,j)*beta_nd(i,j) )/dxsq

             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2))

       if (comp .eq. 1) then

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + 4.d0*b/dxsq) &
                     -(phix(i+1,j)+phix(i-1,j)+phix(i,j+1)+phix(i,j-1))*b/dxsq

             end do
          end do

       else if (comp .eq. 2) then

          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 4.d0*b/dxsq) &
                     -(phiy(i,j+1)+phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j))*b/dxsq

             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                   (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_nd(i,j+1)+beta_nd(i,j))/dxsq) &

                   -( 2.d0*phix(i+1,j)*beta(i,j) &
                     +2.d0*phix(i-1,j)*beta(i-1,j) &
                     +phix(i,j+1)*beta_nd(i,j+1) &
                     +phix(i,j-1)*beta_nd(i,j) &

                     +phiy(i,j+1)*beta_nd(i,j+1) &
                     -phiy(i,j)*beta_nd(i,j) &
                     -phiy(i-1,j+1)*beta_nd(i,j+1) &
                     +phiy(i-1,j)*beta_nd(i,j) )/dxsq

          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                   (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_nd(i+1,j)+beta_nd(i,j))/dxsq) &

                   -( 2.d0*phiy(i,j+1)*beta(i,j) &
                     +2.d0*phiy(i,j-1)*beta(i,j-1) &
                     +phiy(i+1,j)*beta_nd(i+1,j) &
                     +phiy(i-1,j)*beta_nd(i,j) &

                     +phix(i+1,j)*beta_nd(i+1,j) &
                     -phix(i,j)*beta_nd(i,j) &
                     -phix(i+1,j-1)*beta_nd(i+1,j) &
                     +phix(i,j-1)*beta_nd(i,j) )/dxsq

          end do
       end do

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2))

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             Lpx(i,j) = phix(i,j)*(alphax(i,j) + 6.d0*b/dxsq) &
                  -(2.d0*phix(i+1,j)+2.d0*phix(i-1,j)+phix(i,j+1)+phix(i,j-1) &
                   +phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*b/dxsq

          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 6.d0*b/dxsq) &
                  -(2.d0*phiy(i,j+1)+2.d0*phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j) &
                   +phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*b/dxsq

          end do
       end do

    else if (visc_type .eq. -3) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                   ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i-1,j)+gamma(i-1,j) &
                    +beta_nd(i,j+1)+beta_nd(i,j) )/dxsq) &

                   -( phix(i+1,j)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phix(i-1,j)*(fourthirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phix(i,j+1)*beta_nd(i,j+1) &
                     +phix(i,j-1)*beta_nd(i,j) &

                     +phiy(i,j+1)*(beta_nd(i,j+1)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i,j)*(beta_nd(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i-1,j+1)*(beta_nd(i,j+1)-twothirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phiy(i-1,j)*(beta_nd(i,j)-twothirds*beta(i-1,j)+gamma(i-1,j)) )/dxsq

          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                   ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i,j-1)+gamma(i,j-1) &
                    +beta_nd(i+1,j)+beta_nd(i,j) )/dxsq) &

                   -( phiy(i,j+1)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phiy(i,j-1)*(fourthirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phiy(i+1,j)*beta_nd(i+1,j) &
                     +phiy(i-1,j)*beta_nd(i,j) &

                     +phix(i+1,j)*(beta_nd(i+1,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i,j)*(beta_nd(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i+1,j-1)*(beta_nd(i+1,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phix(i,j-1)*(beta_nd(i,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) )/dxsq

          end do
       end do

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2))
       c = gamma(lo(1),lo(2))

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             Lpx(i,j) = phix(i,j)*(alphax(i,j)+(14.d0*b/3.d0+2.d0*c)/dxsq) &
                  -((phix(i+1,j)+phix(i-1,j))*(fourthirds*b+c) &
                    +(phix(i,j+1)+phix(i,j-1))*b &
                    +(phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*(onethird*b+c))/dxsq

          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             Lpy(i,j) = phiy(i,j)*(alphay(i,j)+(14.d0*b/3.d0+2.d0*c)/dxsq) &
                  -((phiy(i,j+1)+phiy(i,j-1))*(fourthirds*b+c) &
                    +(phiy(i+1,j)+phiy(i-1,j))*b &
                    +(phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*(onethird*b+c))/dxsq

          end do
       end do

    end if

  end subroutine stag_applyop_2d_work
  
  ! compute Lphi
  subroutine stag_applyop_3d(la,the_bc_level,phi_fc,Lphi_fc,alpha_fc,beta_cc,beta_ed, &
                             gamma_cc,dx,comp_in)
    
    use define_bc_module
    use multifab_physbc_module
    use probin_module, only: visc_type

    type(layout)   , intent(in   ) :: la
    type(bc_level) , intent(in   ) :: the_bc_level
    type(multifab) , intent(in   ) :: phi_fc(:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_fc(:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc     ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:)  ! edge-centered
    type(multifab) , intent(in   ) :: gamma_cc    ! cell-centered
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ), optional :: comp_in
    
    ! local
    integer :: comp,i,dm,ng_p,ng_l,ng_a,ng_b,ng_e,ng_g
    integer :: lo(get_dim(la)), hi(get_dim(la))
    
    real(kind=dp_t), pointer :: ppx(:,:,:,:)
    real(kind=dp_t), pointer :: ppy(:,:,:,:)
    real(kind=dp_t), pointer :: ppz(:,:,:,:)
    real(kind=dp_t), pointer :: lpx(:,:,:,:)
    real(kind=dp_t), pointer :: lpy(:,:,:,:)
    real(kind=dp_t), pointer :: lpz(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:)
    real(kind=dp_t), pointer :: apy(:,:,:,:)
    real(kind=dp_t), pointer :: apz(:,:,:,:)
    real(kind=dp_t), pointer :: bp(:,:,:,:)
    real(kind=dp_t), pointer :: bp1(:,:,:,:)
    real(kind=dp_t), pointer :: bp2(:,:,:,:)
    real(kind=dp_t), pointer :: bp3(:,:,:,:)
    real(kind=dp_t), pointer :: kp(:,:,:,:)

    if ( present(comp_in) .and. abs(visc_type) .ne. 1) then
       call bl_error("stag_applyop_3d: comp not needed for abs(visc_type) != 1")
    end if

    if ( (.not. present(comp_in)) .and. abs(visc_type) .eq. 1) then
       call bl_error("stag_applyop_3d: comp required if abs(visc_type) = 1")
    end if

    if (dx(1) .ne. dx(2) .or. dx(1) .ne. dx(3)) then
       call bl_error("stag_applyop_3d requires the same dx in all directions")
    end if

    if ( present(comp_in) ) then
       comp = comp_in
    else
       comp = -999
    end if

    dm = get_dim(la)

    if (dm .ne. 3) then
       call bl_error("stag_applyop_3d being called with dm != 3")
    end if

    ng_p = phi_fc(1)%ng
    ng_l = Lphi_fc(1)%ng
    ng_a = alpha_fc(1)%ng
    ng_b = beta_cc%ng
    ng_e = beta_ed(1)%ng
    ng_g = gamma_cc%ng

    do i=1,nfabs(phi_fc(1))
       ppx => dataptr(phi_fc(1), i)
       ppy => dataptr(phi_fc(2), i)
       ppz => dataptr(phi_fc(3), i)
       lpx => dataptr(Lphi_fc(1), i)
       lpy => dataptr(Lphi_fc(2), i)
       lpz => dataptr(Lphi_fc(3), i)
       apx => dataptr(alpha_fc(1), i)
       apy => dataptr(alpha_fc(2), i)
       apz => dataptr(alpha_fc(3), i)
       bp  => dataptr(beta_cc, i)
       bp1 => dataptr(beta_ed(1), i)
       bp2 => dataptr(beta_ed(2), i)
       bp3 => dataptr(beta_ed(3), i)
       kp  => dataptr(gamma_cc, i)
       lo = lwb(get_box(phi_fc(1), i))
       hi = upb(get_box(phi_fc(1), i))
       call stag_applyop_3d_work(ppx(:,:,:,1),ppy(:,:,:,1),ppz(:,:,:,1),ng_p, &
                                 lpx(:,:,:,1),lpy(:,:,:,1),lpz(:,:,:,1),ng_l, &
                                 apx(:,:,:,1),apy(:,:,:,1),apz(:,:,:,1),ng_a, &
                                 bp(:,:,:,1),ng_b, &
                                 bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                                 kp(:,:,:,1),ng_g, &
                                 lo,hi,dx,comp)
    end do

    do i=1,dm
       ! set Lphi on physical domain boundaries to zero
       call multifab_physbc_domainvel(Lphi_fc(i),1,i,1,the_bc_level)
    end do

  end subroutine stag_applyop_3d

  subroutine stag_applyop_3d_work(phix,phiy,phiz,ng_p,Lpx,Lpy,Lpz,ng_l, &
                                  alphax,alphay,alphaz,ng_a,beta,ng_b, &
                                  beta_xy,beta_xz,beta_yz,ng_e, &
                                  gamma,ng_g,lo,hi,dx,comp)

    use probin_module, only: visc_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_l,ng_a,ng_b,ng_e,ng_g
    real(kind=dp_t), intent(in   ) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiz(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(inout) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpz(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphaz(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: comp

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds
    real(kind=dp_t) :: b,c

    dxsq = dx(1)**2
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0

    if (visc_type .eq. -1) then

       if (comp .eq. 1) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + &
                        ( beta(i,j,k)+beta(i-1,j,k) &
                         +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                         +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq) &
                        - ( phix(i+1,j,k)*beta(i,j,k) &
                           +phix(i-1,j,k)*beta(i-1,j,k) &
                           +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                           +phix(i,j-1,k)*beta_xy(i,j,k) &
                           +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                           +phix(i,j,k-1)*beta_xz(i,j,k) )/dxsq

                end do
             end do
          end do

       else if (comp .eq. 2) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j-1,k) &
                         +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                         +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq) &
                        - ( phiy(i,j+1,k)*beta(i,j,k) &
                           +phiy(i,j-1,k)*beta(i,j-1,k) &
                           +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                           +phiy(i-1,j,k)*beta_xy(i,j,k) &
                           +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                           +phiy(i,j,k-1)*beta_yz(i,j,k) )/dxsq

                end do
             end do
          end do

       else if (comp .eq. 3) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j,k-1) &
                         +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                         +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq) &
                        - ( phiz(i,j,k+1)*beta(i,j,k) &
                           +phiz(i,j,k-1)*beta(i,j,k-1) &
                           +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                           +phiz(i-1,j,k)*beta_xz(i,j,k) &
                           +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                           +phiz(i,j-1,k)*beta_yz(i,j,k) )/dxsq

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2),lo(3))

       if (comp .eq. 1) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)+1

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 6.d0*b/dxsq) &
                        -( phix(i+1,j,k)+phix(i-1,j,k) &
                          +phix(i,j+1,k)+phix(i,j-1,k) &
                          +phix(i,j,k+1)+phix(i,j,k-1))*b/dxsq

                end do
             end do
          end do

       else if (comp .eq. 2) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                do i=lo(1),hi(1)

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 6.d0*b/dxsq) &
                        -( phiy(i+1,j,k)+phiy(i-1,j,k) &
                          +phiy(i,j+1,k)+phiy(i,j-1,k) &
                          +phiy(i,j,k+1)+phiy(i,j,k-1))*b/dxsq

                end do
             end do
          end do

       else if (comp .eq. 3) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 6.d0*b/dxsq) &
                        -( phiz(i+1,j,k)+phiz(i-1,j,k) &
                          +phiz(i,j+1,k)+phiz(i,j-1,k) &
                          +phiz(i,j,k+1)+phiz(i,j,k-1))*b/dxsq

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j,k) = phix(i,j,k)*( alphax(i,j,k) + &
                     ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                      +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                      +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq ) &

                     -( 2.d0*phix(i+1,j,k)*beta(i,j,k) &
                       +2.d0*phix(i-1,j,k)*beta(i-1,j,k) &
                       +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                       +phix(i,j-1,k)*beta_xy(i,j,k) &
                       +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                       +phix(i,j,k-1)*beta_xz(i,j,k) &
                       
                       +phiy(i,j+1,k)*beta_xy(i,j+1,k) &
                       -phiy(i,j,k)*beta_xy(i,j,k) &
                       -phiy(i-1,j+1,k)*beta_xy(i,j+1,k) &
                       +phiy(i-1,j,k)*beta_xy(i,j,k) &
                       
                       +phiz(i,j,k+1)*beta_xz(i,j,k+1) &
                       -phiz(i,j,k)*beta_xz(i,j,k) &
                       -phiz(i-1,j,k+1)*beta_xz(i,j,k+1) &
                       +phiz(i-1,j,k)*beta_xz(i,j,k) )/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j,k) = phiy(i,j,k)*( alphay(i,j,k) + &
                     ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                      +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                      +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq ) &

                     -( 2.d0*phiy(i,j+1,k)*beta(i,j,k) &
                       +2.d0*phiy(i,j-1,k)*beta(i,j-1,k) &
                       +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                       +phiy(i-1,j,k)*beta_xy(i,j,k) &
                       +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                       +phiy(i,j,k-1)*beta_yz(i,j,k) &
                       
                       +phix(i+1,j,k)*beta_xy(i+1,j,k) &
                       -phix(i,j,k)*beta_xy(i,j,k) &
                       -phix(i+1,j-1,k)*beta_xy(i+1,j,k) &
                       +phix(i,j-1,k)*beta_xy(i,j,k) &
                       
                       +phiz(i,j,k+1)*beta_yz(i,j,k+1) &
                       -phiz(i,j,k)*beta_yz(i,j,k) &
                       -phiz(i,j-1,k+1)*beta_yz(i,j,k+1) &
                       +phiz(i,j-1,k)*beta_yz(i,j,k) )/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Lpz(i,j,k) = phiz(i,j,k)*( alphaz(i,j,k) + &
                     ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                      +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                      +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq ) &

                     -( 2.d0*phiz(i,j,k+1)*beta(i,j,k) &
                       +2.d0*phiz(i,j,k-1)*beta(i,j,k-1) &
                       +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                       +phiz(i-1,j,k)*beta_xz(i,j,k) &
                       +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                       +phiz(i,j-1,k)*beta_yz(i,j,k) &
                       
                       +phix(i+1,j,k)*beta_xz(i+1,j,k) &
                       -phix(i,j,k)*beta_xz(i,j,k) &
                       -phix(i+1,j,k-1)*beta_xz(i+1,j,k) &
                       +phix(i,j,k-1)*beta_xz(i,j,k) &
                       
                       +phiy(i,j+1,k)*beta_yz(i,j+1,k) &
                       -phiy(i,j,k)*beta_yz(i,j,k) &
                       -phiy(i,j+1,k-1)*beta_yz(i,j+1,k) &
                       +phiy(i,j,k-1)*beta_yz(i,j,k) )/dxsq

             end do
          end do
       end do

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),lo(3))

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 8.d0*b/dxsq) &
                     -( 2.d0*phix(i+1,j,k)+2.d0*phix(i-1,j,k) &
                       +phix(i,j+1,k)+phix(i,j-1,k) &
                       +phix(i,j,k+1)+phix(i,j,k-1) &
                       +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                       +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k) )*b/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 8.d0*b/dxsq) &
                     -( 2.d0*phiy(i,j+1,k)+2.d0*phiy(i,j-1,k) &
                       +phiy(i+1,j,k)+phiy(i-1,j,k) &
                       +phiy(i,j,k+1)+phiy(i,j,k-1) &
                       +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                       +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k) )*b/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 8.d0*b/dxsq) &
                     -( 2.d0*phiz(i,j,k+1)+2.d0*phiz(i,j,k-1) &
                       +phiz(i+1,j,k)+phiz(i-1,j,k) &
                       +phiz(i,j+1,k)+phiz(i,j-1,k) &
                       +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j,k-1)+phix(i,j,k-1) &
                       +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i,j+1,k-1)+phiy(i,j,k-1) )*b/dxsq

             end do
          end do
       end do

    else if (visc_type .eq. -3) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j,k) = phix(i,j,k)*( alphax(i,j,k) + &
                     ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                      +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                      +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                      +beta_xz(i,j,k)+beta_xz(i,j,k+1) )/dxsq ) &

                     -( phix(i+1,j,k)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                       +phix(i-1,j,k)*(fourthirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                       +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                       +phix(i,j-1,k)*beta_xy(i,j,k) &
                       +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                       +phix(i,j,k-1)*beta_xz(i,j,k) &
                       
                       +phiy(i,j+1,k)*(beta_xy(i,j+1,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiy(i,j,k)*(beta_xy(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiy(i-1,j+1,k)*(beta_xy(i,j+1,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                       +phiy(i-1,j,k)*(beta_xy(i,j,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                       
                       +phiz(i,j,k+1)*(beta_xz(i,j,k+1)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiz(i,j,k)*(beta_xz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiz(i-1,j,k+1)*(beta_xz(i,j,k+1)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                       +phiz(i-1,j,k)*(beta_xz(i,j,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) )/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j,k) = phiy(i,j,k)*( alphay(i,j,k) + &
                     ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                      +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                      +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                      +beta_yz(i,j,k)+beta_yz(i,j,k+1) )/dxsq ) &

                     -( phiy(i,j+1,k)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                       +phiy(i,j-1,k)*(fourthirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                       +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                       +phiy(i-1,j,k)*beta_xy(i,j,k) &
                       +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                       +phiy(i,j,k-1)*beta_yz(i,j,k) &
                       
                       +phix(i+1,j,k)*(beta_xy(i+1,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phix(i,j,k)*(beta_xy(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phix(i+1,j-1,k)*(beta_xy(i+1,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                       +phix(i,j-1,k)*(beta_xy(i,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                       
                       +phiz(i,j,k+1)*(beta_yz(i,j,k+1)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiz(i,j,k)*(beta_yz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiz(i,j-1,k+1)*(beta_yz(i,j,k+1)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                       +phiz(i,j-1,k)*(beta_yz(i,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) )/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Lpz(i,j,k) = phiz(i,j,k)*( alphaz(i,j,k) + &
                     ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                      +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1) &
                      +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                      +beta_yz(i,j,k)+beta_yz(i,j+1,k) )/dxsq ) &

                     -( phiz(i,j,k+1)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                       +phiz(i,j,k-1)*(fourthirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                       +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                       +phiz(i-1,j,k)*beta_xz(i,j,k) &
                       +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                       +phiz(i,j-1,k)*beta_yz(i,j,k) &
                       
                       +phix(i+1,j,k)*(beta_xz(i+1,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phix(i,j,k)*(beta_xz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phix(i+1,j,k-1)*(beta_xz(i+1,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                       +phix(i,j,k-1)*(beta_xz(i,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                       
                       +phiy(i,j+1,k)*(beta_yz(i,j+1,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiy(i,j,k)*(beta_yz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                       -phiy(i,j+1,k-1)*(beta_yz(i,j+1,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                       +phiy(i,j,k-1)*(beta_yz(i,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) )/dxsq

             end do
          end do
       end do

    else if (visc_type .eq. 3) then

       b = beta(lo(1),lo(2),lo(3))
       c = gamma(lo(1),lo(2),lo(3))

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1

                Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq) &
                     - ((phix(i+1,j,k)+phix(i-1,j,k))*(fourthirds*b+c) &
                        +(phix(i,j+1,k)+phix(i,j-1,k)+phix(i,j,k+1)+phix(i,j,k-1))*b &
                        +( phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                          +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k)) &
                          *(onethird*b+c))/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)

                Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq) &
                     - ((phiy(i,j+1,k)+phiy(i,j-1,k))*(fourthirds*b+c) &
                        +(phiy(i+1,j,k)+phiy(i-1,j,k)+phiy(i,j,k+1)+phiy(i,j,k-1))*b &
                        +( phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                          +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k)) &
                          *(onethird*b+c))/dxsq

             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq) &
                     - ((phiz(i,j,k+1)+phiz(i,j,k-1))*(fourthirds*b+c) &
                        +(phiz(i+1,j,k)+phiz(i-1,j,k)+phiz(i,j+1,k)+phiz(i,j-1,k))*b &
                        +( phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j,k-1)+phix(i,j,k-1) &
                          +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i,j+1,k-1)+phiy(i,j,k-1)) &
                          *(onethird*b+c))/dxsq

             end do
          end do
       end do

    end if

  end subroutine stag_applyop_3d_work

end module stag_applyop_module
