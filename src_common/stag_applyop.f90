module stag_applyop_module

  use ml_layout_module
  use multifab_module
  use probin_common_module, only: visc_type
  use define_bc_module
  use multifab_physbc_module
  use convert_stag_module

  implicit none

  private

  public :: stag_applyop, stag_applyop_level

contains
  
  ! compute Lphi
  subroutine stag_applyop(mla,the_bc_tower,phi_fc,Lphi_fc,alpha_fc, &
                          beta_cc,beta_ed,gamma_cc,theta,dx)

    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: phi_fc(:,:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:,:)  ! face-centered
    type(multifab) , intent(inout) :: alpha_fc(:,:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc(:)    ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:,:)  ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma_cc(:)   ! cell-centered
    real(kind=dp_t), intent(in   ) :: theta,dx(:,:)

    ! local
    integer :: i,n,dm,nlevs
    type(multifab) :: alpha_fc_temp(mla%nlevel,mla%dim)

    dm = mla%dim
    nlevs = mla%nlevel

    ! multiply alpha_fc_temp by theta
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(alpha_fc_temp(n,i),mla%la(n),1,0,i)
          call multifab_copy_c(alpha_fc_temp(n,i),1,alpha_fc(n,i),1,1,0)
          call multifab_mult_mult_s_c(alpha_fc_temp(n,i),1,theta,1,0)
       end do
    end do

    do n=1,nlevs
       call stag_applyop_level(mla%la(n),the_bc_tower%bc_tower_array(n), &
                               phi_fc(n,:),Lphi_fc(n,:),alpha_fc_temp(n,:), &
                               beta_cc(n),beta_ed(n,:),gamma_cc(n),dx(n,:))
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(alpha_fc_temp(n,i))
       end do
    end do
       

  end subroutine stag_applyop

  ! compute Lphi
  subroutine stag_applyop_level(la,the_bc_level,phi_fc,Lphi_fc,alpha_fc,beta_cc, &
                                beta_ed,gamma_cc,dx,color_in)
    
    type(layout)   , intent(in   ) :: la
    type(bc_level) , intent(in   ) :: the_bc_level
    type(multifab) , intent(in   ) :: phi_fc(:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_fc(:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc     ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:)  ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma_cc    ! cell-centered
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ), optional :: color_in
    
    ! local
    integer :: i,dm,ng_p,ng_l,ng_a,ng_b,ng_e,ng_g
    integer :: lo(get_dim(la)), hi(get_dim(la))
    integer :: color
    
    real(kind=dp_t), pointer :: ppx(:,:,:,:)
    real(kind=dp_t), pointer :: ppy(:,:,:,:)
    real(kind=dp_t), pointer :: ppz(:,:,:,:)
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
    real(kind=dp_t), pointer :: kp(:,:,:,:)

    if (dx(1) .ne. dx(2)) then
       call bl_error("stag_applyop_2d requires the same dx in all directions")
    end if

    dm = get_dim(la)

    if (present(color_in)) then
       color = color_in
    else
       color = 0
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
       lpx => dataptr(Lphi_fc(1), i)
       lpy => dataptr(Lphi_fc(2), i)
       apx => dataptr(alpha_fc(1), i)
       apy => dataptr(alpha_fc(2), i)
       bp  => dataptr(beta_cc, i)
       kp  => dataptr(gamma_cc, i)
       lo = lwb(get_box(phi_fc(1), i))
       hi = upb(get_box(phi_fc(1), i))
       select case (dm)
       case (2)
          bnp => dataptr(beta_ed(1), i)
          call stag_applyop_2d(ppx(:,:,1,1),ppy(:,:,1,1),ng_p, &
                               lpx(:,:,1,1),lpy(:,:,1,1),ng_l, &
                               apx(:,:,1,1),apy(:,:,1,1),ng_a, &
                               bp(:,:,1,1),ng_b, &
                               bnp(:,:,1,1),ng_e, &
                               kp(:,:,1,1),ng_g, &
                               lo,hi,dx,color)
       case (3)
       ppz => dataptr(phi_fc(3), i)
       lpz => dataptr(Lphi_fc(3), i)
       apz => dataptr(alpha_fc(3), i)
       bp1 => dataptr(beta_ed(1), i)
       bp2 => dataptr(beta_ed(2), i)
       bp3 => dataptr(beta_ed(3), i)
       call stag_applyop_3d(ppx(:,:,:,1),ppy(:,:,:,1),ppz(:,:,:,1),ng_p, &
                            lpx(:,:,:,1),lpy(:,:,:,1),lpz(:,:,:,1),ng_l, &
                            apx(:,:,:,1),apy(:,:,:,1),apz(:,:,:,1),ng_a, &
                            bp(:,:,:,1),ng_b, &
                            bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                            kp(:,:,:,1),ng_g, &
                            lo,hi,dx,color)

       end select
    end do

    do i=1,dm
       ! set Lphi on physical domain boundaries to zero
       call multifab_physbc_domainvel(Lphi_fc(i),i,the_bc_level,dx)
    end do

  end subroutine stag_applyop_level

  subroutine stag_applyop_2d(phix,phiy,ng_p,Lpx,Lpy,ng_l, &
                             alphax,alphay,ng_a,beta,ng_b,beta_ed,ng_n, &
                             gamma,ng_g,lo,hi,dx,color)

    integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_l,ng_a,ng_b,ng_n,ng_g
    real(kind=dp_t), intent(in   ) ::    phix(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiy(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(inout) ::     Lpx(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpy(lo(1)-ng_l:,lo(2)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(lo(1)-ng_a:,lo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_ed(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(in   ) ::   gamma(lo(1)-ng_g:,lo(2)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: color
    
    ! local
    integer :: i,j

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds
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
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0

    if (visc_type .eq. -1) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     (beta(i,j)+beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))/dxsq) &
                     - ( phix(i+1,j)*beta(i,j) &
                     +phix(i-1,j)*beta(i-1,j) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) )/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     (beta(i,j)+beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))/dxsq) &
                     - ( phiy(i,j+1)*beta(i,j) &
                     +phiy(i,j-1)*beta(i,j-1) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) )/dxsq

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

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + 4.d0*b/dxsq) &
                     -(phix(i+1,j)+phix(i-1,j)+phix(i,j+1)+phix(i,j-1))*b/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 4.d0*b/dxsq) &
                     -(phiy(i,j+1)+phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j))*b/dxsq

             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))/dxsq) &
                     
                     -( 2.d0*phix(i+1,j)*beta(i,j) &
                     +2.d0*phix(i-1,j)*beta(i-1,j) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) &
                     
                     +phiy(i,j+1)*beta_ed(i,j+1) &
                     -phiy(i,j)*beta_ed(i,j) &
                     -phiy(i-1,j+1)*beta_ed(i,j+1) &
                     +phiy(i-1,j)*beta_ed(i,j) )/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))/dxsq) &
                     
                     -( 2.d0*phiy(i,j+1)*beta(i,j) &
                     +2.d0*phiy(i,j-1)*beta(i,j-1) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) &
                     
                     +phix(i+1,j)*beta_ed(i+1,j) &
                     -phix(i,j)*beta_ed(i,j) &
                     -phix(i+1,j-1)*beta_ed(i+1,j) &
                     +phix(i,j-1)*beta_ed(i,j) )/dxsq

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

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + 6.d0*b/dxsq) &
                     -(2.d0*phix(i+1,j)+2.d0*phix(i-1,j)+phix(i,j+1)+phix(i,j-1) &
                     +phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*b/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 6.d0*b/dxsq) &
                     -(2.d0*phiy(i,j+1)+2.d0*phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j) &
                     +phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*b/dxsq

             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do j=lo(2),hi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1)+1,offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i-1,j)+gamma(i-1,j) &
                     +beta_ed(i,j+1)+beta_ed(i,j) )/dxsq) &
                     
                     -( phix(i+1,j)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phix(i-1,j)*(fourthirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) &
                     
                     +phiy(i,j+1)*(beta_ed(i,j+1)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i,j)*(beta_ed(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i-1,j+1)*(beta_ed(i,j+1)-twothirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phiy(i-1,j)*(beta_ed(i,j)-twothirds*beta(i-1,j)+gamma(i-1,j)) )/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i,j-1)+gamma(i,j-1) &
                     +beta_ed(i+1,j)+beta_ed(i,j) )/dxsq) &
                     
                     -( phiy(i,j+1)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phiy(i,j-1)*(fourthirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) &
                     
                     +phix(i+1,j)*(beta_ed(i+1,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i,j)*(beta_ed(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i+1,j-1)*(beta_ed(i+1,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phix(i,j-1)*(beta_ed(i,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) )/dxsq

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

                Lpx(i,j) = phix(i,j)*(alphax(i,j)+(14.d0*b/3.d0+2.d0*c)/dxsq) &
                     -((phix(i+1,j)+phix(i-1,j))*(fourthirds*b+c) &
                     +(phix(i,j+1)+phix(i,j-1))*b &
                     +(phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*(onethird*b+c))/dxsq

             end do
          end do

       end if

       if (do_y) then

          do j=lo(2),hi(2)+1
             ioff = 0
             if ( offset .eq. 2 .and. mod(lo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=lo(1)+ioff,hi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j)+(14.d0*b/3.d0+2.d0*c)/dxsq) &
                     -((phiy(i,j+1)+phiy(i,j-1))*(fourthirds*b+c) &
                     +(phiy(i+1,j)+phiy(i-1,j))*b &
                     +(phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*(onethird*b+c))/dxsq

             end do
          end do

       end if

    end if

  end subroutine stag_applyop_2d

  subroutine stag_applyop_3d(phix,phiy,phiz,ng_p,Lpx,Lpy,Lpz,ng_l, &
                             alphax,alphay,alphaz,ng_a,beta,ng_b, &
                             beta_xy,beta_xz,beta_yz,ng_e, &
                             gamma,ng_g,lo,hi,dx,color)

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
    integer        , intent(in   ) :: color

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds
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
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0
    
    if (visc_type .eq. -1) then

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

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

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 6.d0*b/dxsq) &
                        -( phix(i+1,j,k)+phix(i-1,j,k) &
                        +phix(i,j+1,k)+phix(i,j-1,k) &
                        +phix(i,j,k+1)+phix(i,j,k-1))*b/dxsq

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

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 6.d0*b/dxsq) &
                        -( phiy(i+1,j,k)+phiy(i-1,j,k) &
                        +phiy(i,j+1,k)+phiy(i,j-1,k) &
                        +phiy(i,j,k+1)+phiy(i,j,k-1))*b/dxsq

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

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 6.d0*b/dxsq) &
                        -( phiz(i+1,j,k)+phiz(i-1,j,k) &
                        +phiz(i,j+1,k)+phiz(i,j-1,k) &
                        +phiz(i,j,k+1)+phiz(i,j,k-1))*b/dxsq

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

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

       end if

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),lo(3))

       if (do_x) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1)+1,offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 8.d0*b/dxsq) &
                        -( 2.d0*phix(i+1,j,k)+2.d0*phix(i-1,j,k) &
                        +phix(i,j+1,k)+phix(i,j-1,k) &
                        +phix(i,j,k+1)+phix(i,j,k-1) &
                        +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k) )*b/dxsq

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

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 8.d0*b/dxsq) &
                        -( 2.d0*phiy(i,j+1,k)+2.d0*phiy(i,j-1,k) &
                        +phiy(i+1,j,k)+phiy(i-1,j,k) &
                        +phiy(i,j,k+1)+phiy(i,j,k-1) &
                        +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k) )*b/dxsq

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

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 8.d0*b/dxsq) &
                        -( 2.d0*phiz(i,j,k+1)+2.d0*phiz(i,j,k-1) &
                        +phiz(i+1,j,k)+phiz(i-1,j,k) &
                        +phiz(i,j+1,k)+phiz(i,j-1,k) &
                        +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j,k-1)+phix(i,j,k-1) &
                        +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i,j+1,k-1)+phiy(i,j,k-1) )*b/dxsq

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

       end if

       if (do_y) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

       end if

       if (do_z) then

          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(lo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=lo(1)+ioff,hi(1),offset

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

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq) &
                        - ((phix(i+1,j,k)+phix(i-1,j,k))*(fourthirds*b+c) &
                        +(phix(i,j+1,k)+phix(i,j-1,k)+phix(i,j,k+1)+phix(i,j,k-1))*b &
                        +( phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k)) &
                        *(onethird*b+c))/dxsq

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

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k)+(20.d0*b/3.d0+2.d0*c)/dxsq) &
                        - ((phiy(i,j+1,k)+phiy(i,j-1,k))*(fourthirds*b+c) &
                        +(phiy(i+1,j,k)+phiy(i-1,j,k)+phiy(i,j,k+1)+phiy(i,j,k-1))*b &
                        +( phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k)) &
                        *(onethird*b+c))/dxsq

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

    end if

  end subroutine stag_applyop_3d

end module stag_applyop_module
