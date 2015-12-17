module advance_reaction_module

  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use bc_module
  use chemical_rates_module
  use advance_reaction_SSA_module
  use probin_reactdiff_module, only: nspecies, nreactions, reaction_type 

  implicit none

  private

  public :: advance_reaction

  ! here we use Mattingly's predictor-corrector with theta=0.5d0 (for rection_type=1).
  ! with these parameters this is actually equivalent to a traditional midpoint scheme.
  real(kind=dp_t), parameter :: theta = 0.5d0
  real(kind=dp_t), parameter :: alpha1 = 2.d0
  real(kind=dp_t), parameter :: alpha2 = 1.d0

contains

  ! this solves dn/dt = f(n) - g (note the minus sign for g)
  !  where f(n) are the chemical production rates (deterministic or stochastic)
  !  and g=ext_src_in is an optional, constant (in time) *deterministic* source term.
  ! to model stochastic particle production (sources) include g in the definition of f instead.

  subroutine advance_reaction(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src_in(:)

    ! local
    type(multifab) :: ext_src(mla%nlevel)
    type(multifab) :: rate1(mla%nlevel)  ! only used for reaction_type=0,1
    type(multifab) :: rate2(mla%nlevel)  ! only used for reaction_type=1

    integer :: nlevs, dm, n

    type(bl_prof_timer),save :: bpt

    real(kind=dp_t) :: mattingly_lin_comb_coef(1:2)  ! only used for reaction_type=1

    !!!!!!!!
    ! init !
    !!!!!!!!

    nlevs = mla%nlevel
    dm = mla%dim

    ! if there are no reactions to process, copy n_old to n_new,
    ! account for ext_src_in (if present) and return
    if(nreactions<1) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
       if (present(ext_src_in)) then
          call multifab_saxpy_3(n_new(n),-dt,ext_src_in(n))
       end if
       return
    end if
    
    call build(bpt,"advance_reaction")

    ! build
    do n=1,nlevs
      call multifab_build(ext_src(n),mla%la(n),nspecies,0)

      if (reaction_type .eq. 0) then  ! first-order tau-leaping or CLE
        call multifab_build(rate1(n),mla%la(n),nspecies,0)
      else if (reaction_type .eq. 1) then  ! second-order tau-leaping or CLE
        call multifab_build(rate1(n),mla%la(n),nspecies,0)
        call multifab_build(rate2(n),mla%la(n),nspecies,0)
      end if
    end do

    ! ext_src
    if (present(ext_src_in)) then
      do n=1,nlevs
        call multifab_copy_c(ext_src(n),1,ext_src_in(n),1,nspecies,0)
      end do
    else
      do n=1,nlevs
        call multifab_setval(ext_src(n),0.d0,all=.true.)
      end do
    end if

    !!!!!!!!!!!!!!!!!!
    ! advancing time !
    !!!!!!!!!!!!!!!!!!

    if (reaction_type .eq. 0) then  ! first-order tau-leaping or CLE 

      ! calculate rates
      ! rates could be deterministic or stochastic depending on use_Poisson_rng
      call chemical_rates(mla,n_old,rate1,dx,dt)

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,rate1(n))
        call multifab_saxpy_3(n_new(n),-dt,ext_src(n))  ! note the negative sign

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else if (reaction_type .eq. 1) then  ! second-order tau-leaping or CLE 

      ! calculate rates from a(n_old)
      call chemical_rates(mla,n_old,rate1,dx,theta*dt)

      ! predictor
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),theta*dt,rate1(n))
        call multifab_saxpy_3(n_new(n),-theta*dt,ext_src(n))  ! note the negative sign

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

      !!!!!!!!!!!!!!!
      ! second step !
      !!!!!!!!!!!!!!!

      mattingly_lin_comb_coef(1) = -alpha2
      mattingly_lin_comb_coef(2) = alpha1 

      ! calculate rates from 2*a(n_pred)-a(n_old)
      call chemical_rates(mla,n_old,rate2,dx,(1.d0-theta)*dt,n_new,mattingly_lin_comb_coef)

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),theta*dt,rate1(n))
        call multifab_saxpy_3(n_new(n),(1.d0-theta)*dt,rate2(n))
        call multifab_saxpy_3(n_new(n),-dt,ext_src(n))  ! note the negative sign
        ! also note that ext_src does not change in the time interval (t,t+dt) 

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else if (reaction_type .eq. 2) then  ! SSA

      call advance_reaction_SSA(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src)

    else

      call bl_error("advance_reaction: invalid reaction_type")

    end if

    !!!!!!!!!!!
    ! destroy !
    !!!!!!!!!!!

    do n=1,nlevs
      call multifab_destroy(ext_src(n))

      if (reaction_type .eq. 0) then  ! first-order tau-leaping or CLE 
        call multifab_destroy(rate1(n))
      else if (reaction_type .eq. 1) then  ! second-order tau-leaping or CLE 
        call multifab_destroy(rate1(n))
        call multifab_destroy(rate2(n))
      end if
    end do

    call destroy(bpt)

  end subroutine advance_reaction

end module advance_reaction_module
