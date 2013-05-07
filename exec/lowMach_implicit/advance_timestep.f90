module advance_timestep_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use convert_stag_module
  use convert_variables_module
  use mk_advective_fluxdiv_module
  use mk_diffusive_fluxdiv_module
  use mk_stochastic_fluxdiv_module
  use gmres_module
  use init_module
  use div_and_grad_module
  use probin_lowmach_module, only: nscal, rhobar, diff_coef, visc_coef
  use probin_common_module, only: fixed_dt

  implicit none

  private

  public :: advance_timestep

contains

  subroutine advance_timestep(mla,mold,mnew,umac,sold,snew,prim,pres,chi,eta,kappa, &
                              rhoc_d_fluxdiv,rhoc_s_fluxdiv,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: mold(:,:)
    type(multifab) , intent(inout) :: mnew(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: prim(:)
    type(multifab) , intent(inout) :: pres(:)
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: rhoc_d_fluxdiv(:)
    type(multifab) , intent(inout) :: rhoc_s_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab) ::    s_update(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::       dpres(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)

    type(multifab) ::        s_fc(mla%nlevel,mla%dim)
    type(multifab) :: gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) :: m_a_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) :: m_d_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) :: m_s_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::       dumac(mla%nlevel,mla%dim)
    type(multifab) ::    umac_old(mla%nlevel,mla%dim)
    type(multifab) ::       gradp(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: S_fac

    nlevs = mla%nlevel
    dm = mla%dim

    S_fac = (1.d0/rhobar(1) - 1.d0/rhobar(2))
    
    do n=1,nlevs
       call multifab_build(   s_update(n),mla%la(n),nscal,0)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1    ,0)
       call multifab_build(      dpres(n),mla%la(n),1    ,1)
       call multifab_build(       divu(n),mla%la(n),1    ,0)
       do i=1,dm
          call multifab_build_edge(       s_fc(n,i),mla%la(n),nscal,1,i)
          call multifab_build_edge(gmres_rhs_v(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_a_fluxdiv(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_d_fluxdiv(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(m_s_fluxdiv(n,i),mla%la(n),1    ,0,i)
          call multifab_build_edge(      dumac(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(   umac_old(n,i),mla%la(n),1    ,1,i)
          call multifab_build_edge(      gradp(n,i),mla%la(n),1    ,0,i)
       end do
    end do

    do n=1,nlevs
       call setval(s_update(n),0.d0,all=.true.)
       do i=1,dm
          call setval(m_a_fluxdiv(n,i),0.d0,all=.true.)
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
          call setval(m_s_fluxdiv(n,i),0.d0,all=.true.)
          call setval(      dumac(n,i),0.d0,all=.true.)
       end do
    end do

    ! make a temporary copy of umac at t^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_old(n,i),1,umac(n,i),1,1,1)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Forward-Euler Scalar Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! average sold = s^n to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,sold,s_fc,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! add A^n for s to s_update
    call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx)

    ! add D^n  for rho1 to s_update
    ! add St^n for rho1 to s_update
    do n=1,nlevs
       call multifab_plus_plus_c(s_update(n),2,rhoc_d_fluxdiv(n),1,1,0)
       call multifab_plus_plus_c(s_update(n),2,rhoc_s_fluxdiv(n),1,1,0)
    end do

    ! snew = s^{*,n+1} = s^n + dt * (A^n + D^n + St^n)
    do n=1,nlevs
       call saxpy(snew(n),1.d0,sold(n),fixed_dt,s_update(n))
       call multifab_fill_boundary(snew(n))
    end do

    ! compute prim^{*,n+1} from s^{*,n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Crank-Nicolson Velocity Predictor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold = m^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! compute s^{*,n+1} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! compute rho^{*,n+1} * v^n and temporarily store it in dumac
    call convert_m_to_umac(mla,s_fc,dumac,umac,.false.)

    do n=1,nlevs
       do i=1,dm

          ! subtract rho^{*,n+1} * v^n from gmres_rhs_v
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,dumac(n,i),1,1,0)

          ! multiply gmres_rhs_v by 1/dt
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/fixed_dt,1,0)

       end do
    end do

    ! compute grad pi^n
    call compute_gradp(mla,pres,gradp,dx)

    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! compute m_a_fluxdiv = A^n for m
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv,dx)

    ! add A^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_d_fluxdiv = D^n for m
    call mk_diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,kappa,dx, &
                                the_bc_tower%bc_tower_array)

    ! add D^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_s_fluxdiv = St^n for m
    call mk_stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,dx)

    ! add S^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! initialize rhs_p for gmres solve to zero since subsequent subroutines will add to it
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! add D^{*,n+1} to rhs_p
    call mk_diffusive_rhoc_fluxdiv(mla,gmres_rhs_p,1,prim,s_fc,chi,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add St^n to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_s_fluxdiv(n),1,1,0)
    end do

    ! reset s_update, then
    ! we save work by saving D^{*,n+1} + St^{*,n+1}
    ! this is used in Step 3 below
    do n=1,nlevs
       call multifab_setval_c(s_update(n),0.d0,1,1,all=.true.)
       call multifab_copy_c(s_update(n),2,gmres_rhs_p(n),1,1,0)
    end do

    ! multiply by -S_fac since gmres solves -div(u)=S
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! compute div(u)
    call compute_divu(mla,umac,divu,dx)

    ! add div(u) to gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,1)
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
          call multifab_setval(dpres(n)  ,0.d0,all=.true.)
       end do
    end do

    ! call gmres to compute v^{*,n+1}
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpres,snew, &
               eta,kappa,1.d0/fixed_dt)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,1)
    end do

    ! no need to compute p^{*,n+1} since we don't use it, 
    ! but keep dpres as initial guess for corrector
    ! compute v^{*,n+1} = v^n + dumac
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! convert v^{*,n+1} to m^{*,n+1}
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(mnew(n,i))
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! s_update already contains D^{*,n+1} + St^{*,n+1} from above
    ! add A^{*,n+1} for s to s_update
    call mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx)

    ! snew = s^{n+1} 
    !      = (1/2)*s^n + (1/2)*s^{*,n+1} + (dt/2)*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1})
    do n=1,nlevs
       call multifab_plus_plus_c(snew(n),1,sold(n),1,nscal,0)
       call multifab_mult_mult_s_c(snew(n),1,0.5d0,nscal,0)
       call multifab_mult_mult_s_c(s_update(n),1,fixed_dt/2.d0,nscal,0)
       call multifab_plus_plus_c(snew(n),1,s_update(n),1,nscal,0)
       call multifab_fill_boundary(snew(n))
    end do

    ! compute prim^{n+1} from s^{n+1} in valid region
    call convert_cons_to_prim(mla,snew,prim,.true.)
    do n=1,nlevs
       call multifab_fill_boundary(prim(n))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 4 - Crank-Nicolson Velocity Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set rhs to mold = m^n
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
       end do
    end do

    ! compute s^{n+1} to faces
    do i=1,nscal
       call average_cc_to_face(nlevs,snew,s_fc,i,dm+2,1,the_bc_tower%bc_tower_array)
    end do

    ! compute rho^{n+1} * v^n and temporarily store it in dumac
    call convert_m_to_umac(mla,s_fc,dumac,umac,.false.)

    do n=1,nlevs
       do i=1,dm

          ! subtract rho^{n+1} * v^n from gmres_rhs_v
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,dumac(n,i),1,1,0)

          ! multiply gmres_rhs_v by 1/dt
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/fixed_dt,1,0)

       end do
    end do

    ! gradp already contains grad pi^n
    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradp(n,i),1,1,0)
       end do
    end do

    ! m_a_fluxdiv already contains A^n for m
    ! add A^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! m_d_fluxdiv already contains D^n for m
    ! add D^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! m_s_fluxdiv already contains S^n for m
    ! add S^n for m to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! recompute chi, eta, and kappa
    call compute_chi(mla,chi,prim,dx)
    call compute_eta(mla,eta,prim,dx)
    call compute_kappa(mla,kappa,prim,dx)

    ! initialize rhs_p for gmres solve to zero since subsequent subroutines will add to it
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
    end do

    ! reset to zero since we only add to them
    do n=1,nlevs
       call setval(rhoc_d_fluxdiv(n),0.d0,all=.true.)
       call setval(rhoc_s_fluxdiv(n),0.d0,all=.true.)
    end do

    ! create D^{n+1}
    call mk_diffusive_rhoc_fluxdiv(mla,rhoc_d_fluxdiv,1,prim,s_fc,chi,dx, &
                                   the_bc_tower%bc_tower_array)

    ! add D^{n+1} to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_d_fluxdiv(n),1,1,0)
    end do

    ! fill the stochastic multifabs with a new set of random numbers
    call fill_stochastic(mla)

    ! create St^{n+1}
    call mk_stochastic_s_fluxdiv(mla,the_bc_tower%bc_tower_array,rhoc_s_fluxdiv,s_fc,chi,dx,1)

    ! add St^{n+1} to rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,rhoc_s_fluxdiv(n),1,1,0)
    end do    

    ! multiply by -S_fac since gmres solves -div(u)=S
    do n=1,nlevs
       call multifab_mult_mult_s_c(gmres_rhs_p(n),1,-S_fac,1,0)
    end do

    ! divu already contains div(u^n)
    ! add div(u) to gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,1)
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
          call multifab_setval(dpres(n)  ,0.d0,all=.true.)
       end do
    end do

    ! call gmres to compute v^{*,n+1}
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpres,snew, &
               eta,kappa,1.d0/fixed_dt)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,1)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,1)
    end do

    ! compute p^{n+1} = p^n + deltap
    ! compute v^{n+1} = v^n + dumac
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac(n,i),1,umac_old(n,i),1,1,0)
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! convert v^{n+1} to m^{n+1}
    call convert_m_to_umac(mla,s_fc,mnew,umac,.false.)
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(mnew(n,i))
       end do
    end do



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call multifab_destroy(s_update(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpres(n))
       call multifab_destroy(divu(n))
       do i=1,dm
          call multifab_destroy(s_fc(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(m_a_fluxdiv(n,i))
          call multifab_destroy(m_d_fluxdiv(n,i))
          call multifab_destroy(m_s_fluxdiv(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_old(n,i))
          call multifab_destroy(gradp(n,i))
       end do
    end do

  end subroutine advance_timestep

end module advance_timestep_module
