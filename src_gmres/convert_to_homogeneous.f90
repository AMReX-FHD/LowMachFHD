module convert_to_homogeneous_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use stag_applyop_module
  use div_and_grad_module

  implicit none

  private

  public :: convert_to_homogeneous
  
contains

  subroutine convert_to_homogeneous(mla,b_u,b_p,alpha_fc,beta,beta_ed,gamma, &
                                    theta,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: b_u(:,:)
    type(multifab) , intent(inout) :: b_p(:)
    type(multifab) , intent(inout) :: alpha_fc(:,:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: beta_ed(:,:)
    type(multifab) , intent(in   ) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: theta,dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm

    type(multifab) :: phi(mla%nlevel,mla%dim)
    type(multifab) :: Lphi(mla%nlevel,mla%dim)
    type(multifab) :: Dphi(mla%nlevel)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(Dphi(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge( phi(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(Lphi(n,i),mla%la(n),1,0,i)
       end do
    end do

  do n=1,nlevs
     do i=1,dm

        ! set phi = 0
        call setval(phi(n,i),0.d0,all=.true.)

        ! set values on physical boundaries
        call multifab_physbc_domainvel(phi(n,i),i,the_bc_tower%bc_tower_array(n),dx(n,:),.true.)

        ! fill periodic ghost cells
        call multifab_fill_boundary(phi(n,i))

        ! fill physical ghost cells
        call multifab_physbc_macvel(phi(n,i),i,the_bc_tower%bc_tower_array(n),dx(n,:),.true.)

     end do
  end do

  ! compute Lphi
  call stag_applyop(mla,the_bc_tower,phi,Lphi,alpha_fc,beta,beta_ed,gamma,theta,dx)

  ! subtract Lphi from b_u
  do n=1,nlevs
     do i=1,dm
        call multifab_sub_sub_c(b_u(n,i),1,Lphi(n,i),1,1,0)
     end do
  end do

  ! compute divergence of phi
  call compute_divu(mla,phi,Dphi,dx)

  ! add Dphi to b_p
  do n=1,nlevs
     call multifab_plus_plus_c(b_p(n),1,Dphi(n),1,1,0)
  end do

    do n=1,nlevs
       call destroy(Dphi(n))
       do i=1,dm
          call destroy(phi(n,i))
          call destroy(Lphi(n,i))
       end do
    end do

  end subroutine convert_to_homogeneous

end module convert_to_homogeneous_module
