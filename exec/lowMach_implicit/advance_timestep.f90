module advance_timestep_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use mk_advective_fluxdiv_module
  use probin_lowmach_module, only: nscal
  use probin_common_module, only: fixed_dt

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,chi,eta,kappa,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: chi(:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: kappa(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    type(multifab) :: s_update(mla%nlevel)
    type(multifab) ::   s_face(mla%nlevel,mla%dim)

    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    nlevs = mla%nlevel
    dm = mla%dim
    
    do n=1,nlevs
       call multifab_build(s_update(n),mla%la(n),nscal,0)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1,0)
       do i=1,dm
          call multifab_build_edge(s_face(n,i),mla%la(n),nscal,0,i)
          call multifab_build_edge(gmres_rhs_v(n,i),mla%la(n),1,0,i)
       end do
    end do

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Forward-Euler Scalar Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! average s to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,sold,s_face,i,dm+2,1,the_bc_level)
    end do

    ! compute del dot (-rho*v)
    call mk_advective_s_fluxdiv(mla,umac,s_face,s_update,dx)

    ! snew = sold + dt * del dot (-rho*v)
    do n=1,nlevs
       call saxpy(snew(n),1.d0,sold(n),fixed_dt,s_update(n))
    end do

    ! fill ghost cells
    do n=1,nlevs
       call multifab_fill_boundary(snew(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Crank-Nicolson Velocity Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! multiply eta and kappa by dt/2
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,fixed_dt/2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,fixed_dt/2.d0,1,1)
    end do

    ! build up the rhs - compute explicit diffusive term and multiply by -dt/2



    ! build up the rhs - add rho^n v^n
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! build up the rhs - add advective term
    call mk_advective_m_fluxdiv(mla,umac,mold,gmres_rhs_v,dx)

    ! call gmres
    

    ! restore eta and kappa by 2/dt
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0/fixed_dt,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0/fixed_dt,1,1)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Crank-Nicolson Velocity Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do n=1,nlevs
       call destroy(s_update(n))
       call destroy(gmres_rhs_p(n))
       do i=1,dm
          call destroy(s_face(n,i))
          call destroy(gmres_rhs_v(n,i))
       end do
    end do

  end subroutine advance_timestep

end module advance_timestep_module
