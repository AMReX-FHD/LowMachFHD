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

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,eta,chi,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: chi(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    type(multifab) :: s_update(mla%nlevel)
    type(multifab) :: s_face(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    nlevs = mla%nlevel
    dm = mla%dim
    
    do n=1,nlevs
       call multifab_build(s_update(n),mla%la(n),nscal,0)
       call setval(s_update(n),0.d0,all=.true.)
       do i=1,dm
          call multifab_build_edge(s_face(n,i),mla%la(n),nscal,0,i)
       end do
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Crank-Nicolson Velocity Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  end subroutine advance_timestep

end module advance_timestep_module
