module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use convert_m_to_umac_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use mk_baro_fluxdiv_module
  use diffusive_rhoc_fluxdiv_module
  use stochastic_rhoc_fluxdiv_module
  use bc_module
  use multifab_physbc_stag_module
  use probin_binarylm_module, only: barodiffusion_type
  use probin_common_module, only: rhobar, algorithm_type

  implicit none

  private

  public :: initial_projection

  ! special inhomogeneous boundary condition multifab
  ! vel_bc_n(nlevs,dm) are the normal velocities
  ! in 2D, vel_bc_t(nlevs,2) respresents
  !   1. y-velocity bc on x-faces (nodal)
  !   2. x-velocity bc on y-faces (nodal)
  ! in 3D, vel_bc_t(nlevs,6) represents
  !   1. y-velocity bc on x-faces (nodal in y and x)
  !   2. z-velocity bc on x-faces (nodal in z and x)
  !   3. x-velocity bc on y-faces (nodal in x and y)
  !   4. z-velocity bc on y-faces (nodal in z and y)
  !   5. x-velocity bc on z-faces (nodal in x and z)
  !   6. y-velocity bc on z-faces (nodal in y and z)
  type(multifab), allocatable, save :: vel_bc_n(:,:)
  type(multifab), allocatable, save :: vel_bc_t(:,:)

contains

  subroutine initial_projection(mla,mold,umac,sold,s_fc,prim,eta_ed,chi_fc,gp_fc, &
                                rhoc_fluxdiv,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: prim(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: chi_fc(:,:)
    type(multifab) , intent(in   ) :: gp_fc(:,:)
    type(multifab) , intent(inout) :: rhoc_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,dm,n,nlevs
    real(kind=dp_t) :: S_fac

    type(multifab) ::   mac_rhs(mla%nlevel)
    type(multifab) ::      divu(mla%nlevel)
    type(multifab) ::       phi(mla%nlevel)
    type(multifab) :: rhoinv_fc(mla%nlevel,mla%dim)

    real(kind=dp_t), allocatable :: weights(:)

    if (algorithm_type .eq. 0 .or. algorithm_type .eq. 1) then
       allocate(weights(1))
       weights(1) = 1.d0
    else if (algorithm_type .eq. 2) then
       allocate(weights(2))
       weights(1) = 1.d0
       weights(2) = 0.d0
    end if

    dm = mla%dim
    nlevs = mla%nlevel

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    call build_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_build(mac_rhs(n),mla%la(n),1,0)
       call multifab_build(divu(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(rhoinv_fc(n,i),mla%la(n),1    ,0,i)
       end do       
    end do

    do n=1,nlevs
       call setval(mac_rhs(n),0.d0)
       call setval(phi(n),0.d0)
    end do

    ! average sold to faces
    call average_cc_to_face(nlevs,sold,s_fc,1,scal_bc_comp,2,the_bc_tower%bc_tower_array)

    ! convert m^init to v^init in valid region
    call convert_m_to_umac(mla,s_fc,mold,umac,.true.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! build rhs = div(v^init) - S^0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,0.d0, &
                                   the_bc_tower%bc_tower_array)

    ! set rhoc_d_fluxdiv = div(rho*chi grad c)^0
    call diffusive_rhoc_fluxdiv(mla,rhoc_fluxdiv,1,prim,s_fc,chi_fc,dx, &
                                the_bc_tower%bc_tower_array,vel_bc_n)

    if (barodiffusion_type .gt. 0) then
       ! compute baro-diffusion flux divergence
       call mk_baro_fluxdiv(mla,rhoc_fluxdiv,1,s_fc,chi_fc,gp_fc,dx, &
                            the_bc_tower%bc_tower_array,vel_bc_n)
    end if

    ! set rhoc_s_fluxdiv = div(Psi^0)
    call stochastic_rhoc_fluxdiv(mla,the_bc_tower%bc_tower_array,rhoc_fluxdiv,s_fc, &
                                 chi_fc,dx,dt,vel_bc_n,weights)

    ! add diffusive, baro-diffusion, and stochastic rhoc fluxes to to mac_rhs
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,rhoc_fluxdiv(n),1,1,0)
    end do

    do n=1,nlevs
       do i=1,dm
          ! to deal with reservoirs
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! multiply mac_rhs by -S_fac
    ! now mac_rhs = -S^0
    do n=1,nlevs
       call multifab_mult_mult_s_c(mac_rhs(n),1,-S_fac,1,0)
    end do

    ! set divu = div(v^init)
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div(v^init) to mac_rhs
    ! now mac_rhs = div(v^init) - S^0
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,divu(n),1,1,0)
    end do

    ! compute (1/rho^0)
    do n=1,nlevs
       do i=1,dm
          call setval(rhoinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhoinv_fc(n,i),1,s_fc(n,i),1,1,0)
       end do
    end do

    ! solve div (1/rho^0) grad phi = div(v^init) - S^0
    ! solve to completion, i.e., use the 'full' solver
    call macproject(mla,phi,umac,rhoinv_fc,mac_rhs,dx,the_bc_tower,.true.)

    ! v^0 = v^init - (1/rho^0) grad phi
    call subtract_weighted_gradp(mla,umac,rhoinv_fc,phi,dx,the_bc_tower)

    ! fill ghost cells
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! compute mold by converting v^0 to m^0 in valid plus ghost region
    ! now mold has properly filled ghost cells
    call convert_m_to_umac(mla,s_fc,mold,umac,.false.)

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhoinv_fc(n,i))
       end do
    end do

  end subroutine initial_projection

  subroutine build_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs
    logical :: nodal_temp(3)

    dm = mla%dim
    nlevs = mla%nlevel

    allocate(vel_bc_n(nlevs,dm))
    if (dm .eq. 2) then
       allocate(vel_bc_t(nlevs,2))
    else if (dm .eq. 3) then
       allocate(vel_bc_t(nlevs,6))
    end if

    do n=1,nlevs
       ! boundary conditions
       do i=1,dm
          call multifab_build_edge(vel_bc_n(n,i),mla%la(n),1,0,i)
       end do
       if (dm .eq. 2) then
          ! y-velocity bc on x-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,1),mla%la(n),1,0)
          ! x-velocity bc on y-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,2),mla%la(n),1,0)
       else
          ! y-velocity bc on x-faces (nodal in y and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on x-faces (nodal in z and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on y-faces (nodal in x and y)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on y-faces (nodal in z and y)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on z-faces (nodal in x and z)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
          ! y-velocity bc on z-faces (nodal in y and z)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
       end if

       do i=1,dm
          call multifab_setval(vel_bc_n(n,i),0.d0,all=.true.)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_setval(vel_bc_t(n,i),0.d0,all=.true.)
       end do

    end do

  end subroutine build_bc_multifabs

  subroutine destroy_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,dm          
          call multifab_destroy(vel_bc_n(n,i))
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_destroy(vel_bc_t(n,i))
       end do
    end do

    deallocate(vel_bc_n,vel_bc_t)

  end subroutine destroy_bc_multifabs

end module initial_projection_module
