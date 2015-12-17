module advance_reaction_SSA_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use compute_reaction_rates_module
  use probin_reactdiff_module, only: nspecies, nreactions, stoichiometric_factors, cross_section

  implicit none

  private

  public :: advance_reaction_SSA
  
contains

  ! advance_reaction solves dn/dt = f(n) - g (note the minus sign for g)
  !  where f(n) are the chemical production rates (deterministic or stochastic)
  !  and g=ext_src_in is an optional, constant (in time) *deterministic* source term.
  ! to model stochastic particle production (sources) include g in the definition of f instead.

  subroutine advance_reaction_SSA(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src_in(:)

    ! local
    integer         :: nlevs, dm, n, i
    real(kind=dp_t) :: dv

    type(mfiter)    :: mfi
    type(box)       :: tilebox
    integer         :: tlo(mla%dim), thi(mla%dim)
    integer         ::  lo(mla%dim),  hi(mla%dim)
    integer         :: ng_o, ng_n, ng_e

    real(kind=dp_t), pointer :: op(:,:,:,:) ! "o" for n_old
    real(kind=dp_t), pointer :: np(:,:,:,:) ! "n" for n_new
    real(kind=dp_t), pointer :: ep(:,:,:,:) ! "e" for ext_src

    type(multifab) :: ext_src(mla%nlevel)

    type(bl_prof_timer),save :: bpt

    ! this is a routine which is called by advance_reaction when reaction_type=2.
    ! note that the contribution of ext_src to n_new has not be implemented yet.
    ! when implementing it, remind the role and sign convention for ext_src in advance_reaction
    if (present(ext_src_in)) then
       call bl_error("advance_reaction_SSA does not use ext_src")
    end if

    nlevs = mla%nlevel
    dm = mla%dim
   
    ! no reaction case is already checked before this routine is called by advance_reaction 
    ! hence, the following condition should not hold
    if (nreactions .lt. 1) then
      call bl_error("advance_reaction: nreactions should be .ge. 1")
    end if  

    ! otherwise, complete the remaining part
    call build(bpt,"advance_reaction_SSA")

    do n=1,nlevs
       call multifab_build(ext_src(n),mla%la(n),nspecies,0)
    end do

    if (present(ext_src_in)) then
       do n=1,nlevs
          call multifab_copy_c(ext_src(n),1,ext_src_in(n),1,nspecies,0)
       end do
    else
       do n=1,nlevs
          call multifab_setval(ext_src(n),0.d0,all=.true.)
       end do
    end if

    dv = product(dx(1,1:dm))
    if (dm<3) dv = dv*cross_section


    !!!!! omp tiling

    ng_o = n_old(1)%ng
    ng_n = n_new(1)%ng
    ng_e = ext_src(1)%ng

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(op,np,ep,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi,n_old(n),tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          op => dataptr(n_old(n),i)
          np => dataptr(n_new(n),i)
          ep => dataptr(ext_src(n),i)

          lo = lwb(get_box(n_new(n),i))
          hi = upb(get_box(n_new(n),i))
          select case (dm)
          case (2)
             call advance_reaction_SSA_2d(op(:,:,1,:),ng_o,np(:,:,1,:),ng_n,ep(:,:,1,:),ng_e, &
                                      lo,hi,tlo,thi,dv,dt)
          case (3)
             call advance_reaction_SSA_3d(op(:,:,:,:),ng_o,np(:,:,:,:),ng_n,ep(:,:,:,:),ng_e, &
                                      lo,hi,tlo,thi,dv,dt)
          end select
       end do
    end do
    !$omp end parallel


    ! ensure ghost cells are consistent for n_new
    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies,             &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    do n=1,nlevs
       call multifab_destroy(ext_src(n))
    end do

    call destroy(bpt)

  end subroutine advance_reaction_SSA
  
  
  subroutine advance_reaction_SSA_2d(n_old,ng_o,n_new,ng_n,ext_src,ng_e,glo,ghi,tlo,thi,dv,dt)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_o, ng_n, ng_e
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: ext_src(glo(1)-ng_e:,glo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    
    ! local
    integer :: i, j

    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_SSA_cell(n_old(i,j,1:nspecies),n_new(i,j,1:nspecies), &
                                  ext_src(i,j,1:nspecies),dv,dt)
    end do
    end do

  end subroutine advance_reaction_SSA_2d


  subroutine advance_reaction_SSA_3d(n_old,ng_o,n_new,ng_n,ext_src,ng_e,glo,ghi,tlo,thi,dv,dt)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_o, ng_n, ng_e
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,glo(3)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: ext_src(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    
    ! local
    integer :: i, j, k

    do k=tlo(3),thi(3)
    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_SSA_cell(n_old(i,j,k,1:nspecies),n_new(i,j,k,1:nspecies), &
                                  ext_src(i,j,k,1:nspecies),dv,dt)
    end do
    end do
    end do

  end subroutine advance_reaction_SSA_3d


  subroutine advance_reaction_SSA_cell(n_old,n_new,ext_src,dv,dt)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: ext_src(:)
    real(kind=dp_t), intent(in   ) :: dv, dt

    real(kind=dp_t) :: avg_react_rate(1:nreactions)
    real(kind=dp_t) :: rTotal, rr, rSum, tau, t_local

    integer :: spec,reaction, which_reaction
    integer :: n_steps_SSA
    
    ! copy old state into new
    n_new = n_old

    t_local = 0.d0
    n_steps_SSA = 0

    EventLoop: do

       ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
       call compute_reaction_rates(n_new(1:nspecies),avg_react_rate,dv)

       ! compute reaction rates in units (# reactions) / (unit time)
       avg_react_rate = max(0.0d0,avg_react_rate*dv)

       ! sum the reaction rates
       rTotal = sum(avg_react_rate(1:nreactions))

       ! generate pseudorandom number in interval [0,1).
       call UniformRNG(rr)
       ! tau is how long until the next reaction occurs
       tau = -log(1-rr)/rTotal
       t_local = t_local + tau;

       if (t_local .gt. dt) exit EventLoop

       ! Select the next reaction according to relative rates
       call UniformRNG(rr)
       rr = rr*rTotal
       rSum = 0
       FindReaction: do reaction=1,nreactions
          rSum = rSum + avg_react_rate(reaction)
          which_reaction = reaction
          if( rSum >= rr ) exit FindReaction
       end do FindReaction

       ! update number densities for this reaction
       do spec=1,nspecies
          n_new(spec) = n_new(spec) + &
               (stoichiometric_factors(spec,2,which_reaction)-stoichiometric_factors(spec,1,which_reaction)) / dv
       end do
       
       n_steps_SSA = n_steps_SSA+1

    end do EventLoop

  end subroutine advance_reaction_SSA_cell

end module advance_reaction_SSA_module
