module advance_timestep_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use mk_advective_fluxdiv_module
  use mk_diffusive_fluxdiv_module
  use gmres_module
  use probin_lowmach_module, only: nscal
  use probin_common_module, only: fixed_dt, theta_fac

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,prim,chi,eta,kappa,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(in   ) :: chi(:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: kappa(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::         phi(mla%nlevel)

    type(multifab) ::      s_face(mla%nlevel,mla%dim)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::    chi_face(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv(mla%nlevel,mla%dim)

    type(multifab) :: eta_nodal(mla%nlevel)   ! averaged to nodes (2D only)
    type(multifab) ::  eta_edge(mla%nlevel,3) ! averaged to edges (3D only; xy/xz/yz edges)

    logical :: nodal_temp(mla%dim)

    integer :: i,dm,n,nlevs

    nlevs = mla%nlevel
    dm = mla%dim
    
    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),nscal,0)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1    ,0)
       call multifab_build(        phi(n),mla%la(n),1    ,1)
       do i=1,dm
          call multifab_build_edge(     s_face(n,i),mla%la(n),nscal,1,i)
          call multifab_build_edge(gmres_rhs_v(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(   chi_face(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_a_fluxdiv(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_d_fluxdiv(n,i),mla%la(n),1    ,0,i)
       end do
    end do

    ! nodal (in 2D) and edge-based (in 3D) eta
    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_build_nodal(eta_nodal(n),mla%la(n),1,0)
       end do
    else
       do n=1,nlevs
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(eta_edge(n,1),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(eta_edge(n,2),mla%la(n),1,0,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(eta_edge(n,3),mla%la(n),1,0,nodal_temp)
       end do
    end if

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
       call setval(phi(n)     ,0.d0,all=.true.)
       do i=1,dm
          call setval(m_a_fluxdiv(n,i),0.d0,all=.true.)
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Forward-Euler Scalar Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! compute prim from sold in valid region
    call convert_cons_to_prim(mla,sold,prim,.true.)

    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
    end do

    ! average sold to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,sold,s_face,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! average chi to faces
    call average_cc_to_face(nlevs,chi,chi_face,1,dm+2,1,the_bc_tower%bc_tower_array)

    ! compute advective flux divergence
    call mk_advective_s_fluxdiv(mla,umac,s_face,s_update,dx)

    ! compute del dot rho chi grad c
    call mk_diffusive_rhoc_fluxdiv(mla,s_update,2,prim,s_face,chi_face,dx, &
                                   the_bc_tower%bc_tower_array)

    ! snew = sold + dt * del dot (A+D)
    do n=1,nlevs
       call saxpy(snew(n),1.d0,sold(n),fixed_dt,s_update(n))
       call multifab_fill_boundary(snew(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Crank-Nicolson Velocity Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute eta on nodes (2D) or edges (3D)
    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta,eta_nodal,1,dm+2,1,the_bc_tower%bc_tower_array)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta,eta_edge,1,dm+2,1,the_bc_tower%bc_tower_array)
    end if

    ! build up the rhs_v - set rhs to rho^n v^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! compute diffusive flux divergence for m
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,eta_nodal,eta_edge, &
                                kappa,dx,the_bc_tower%bc_tower_array)

    ! multiply by dt/2 and add to rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv(n,i),1,fixed_dt/2.d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute advective flux divergence for m
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv,dx)

    ! multiply by dt and add to rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv(n,i),1,fixed_dt,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! initialize rhs_p to zero since subsequent subroutines will add to it
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! add del dot rho chi grad c to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_face,chi_face,dx, &
                                   the_bc_tower%bc_tower_array)

    ! multiply eta and kappa by dt/2
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,fixed_dt/2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,fixed_dt/2.d0,1,1)
    end do

    ! call gmres to compute umac^{n+1,*}
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,umac,phi,snew, &
               eta,kappa,theta_fac)

    ! restore eta and kappa by 2/dt
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0/fixed_dt,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0/fixed_dt,1,1)
    end do

    ! convert umac to mnew
    call convert_m_to_umac(mla,s_face,mnew,umac,.false.)

    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(mnew(n,i))
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute prim from snew in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)

    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
    end do

    ! average snew to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_face,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! reset s_update
    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
    end do
    
    ! compute advective flux divergence
    call mk_advective_s_fluxdiv(mla,umac,s_face,s_update,dx)

    ! compute del dot rho chi grad c
    call mk_diffusive_rhoc_fluxdiv(mla,s_update,2,prim,s_face,chi_face,dx, &
                                   the_bc_tower%bc_tower_array)

    ! snew = (1/2)*snew + (1/2)*sold + (1/2)*del dot (A+D)
    do n=1,nlevs
       call multifab_mult_mult_s_c(snew(n),1,0.5d0,nscal,0)
       call multifab_mult_mult_s_c(sold(n),1,0.5d0,nscal,0)
       call multifab_plus_plus_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_mult_mult_s_c(s_update(n),1,fixed_dt/2.d0,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
       call multifab_fill_boundary(snew(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Crank-Nicolson Velocity Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up the rhs_v - set rhs to rho^n v^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! add (dt/2) * old diffusive flux divergence to rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! multiply dt * old advective flux divergence by 1/2 and add to rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! reset m_a_fluxdiv
    do n=1,nlevs
       do i=1,dm
          call setval(m_a_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do

    ! compute advective flux divergence for m
    call mk_advective_m_fluxdiv(mla,umac,mnew,m_a_fluxdiv,dx)

    ! multiply by dt/2 and add to rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv(n,i),1,fixed_dt/2.d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! initialize rhs_p to zero since subsequent subroutines will add to it
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! add del dot rho chi grad c to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_face,chi_face,dx, &
                                   the_bc_tower%bc_tower_array)

    ! multiply eta and kappa by dt/2
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,fixed_dt/2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,fixed_dt/2.d0,1,1)
    end do

    ! call gmres to compute umac^{n+1,*}
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,umac,phi,snew, &
               eta,kappa,theta_fac)

    ! restore eta and kappa by 2/dt
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0/fixed_dt,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0/fixed_dt,1,1)
    end do

    ! convert umac to mnew
    call convert_m_to_umac(mla,s_face,mnew,umac,.false.)

    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(mnew(n,i))
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call destroy(s_update(n))
       call destroy(gmres_rhs_p(n))
       call destroy(phi(n))
       do i=1,dm
          call destroy(s_face(n,i))
          call destroy(gmres_rhs_v(n,i))
          call destroy(chi_face(n,i))
          call destroy(m_a_fluxdiv(n,i))
          call destroy(m_d_fluxdiv(n,i))
       end do
    end do

    if (dm .eq. 2) then
       do n=1,nlevs
          call multifab_destroy(eta_nodal(n))
       end do
    else
       do n=1,nlevs
          do i=1,3
             call multifab_destroy(eta_edge(n,i))
          end do
       end do
    end if

  end subroutine advance_timestep

end module advance_timestep_module
