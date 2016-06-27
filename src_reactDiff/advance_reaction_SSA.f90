module advance_reaction_SSA_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use bl_rng_module
  use bl_random_module
  use compute_reaction_rates_module
  use probin_reactdiff_module, only: nspecies, nreactions, stoichiometric_factors, &
                                     cross_section, use_bl_rng

  implicit none

  private

  public :: advance_reaction_SSA
  
contains

  ! advance_reaction_SSA solves dn/dt = f(n) 
  !  where f(n) are the chemical production rates (deterministic or stochastic)
  subroutine advance_reaction_SSA(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:) 
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer         :: nlevs, dm, n, i
    real(kind=dp_t) :: dv

    type(mfiter)    :: mfi
    type(box)       :: tilebox
    integer         :: tlo(mla%dim), thi(mla%dim)
    integer         ::  lo(mla%dim),  hi(mla%dim)
    integer         :: ng_o, ng_n

    real(kind=dp_t), pointer :: op(:,:,:,:) ! "o" for n_old
    real(kind=dp_t), pointer :: np(:,:,:,:) ! "n" for n_new

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim

    ! no reaction case is already checked before this routine is called by advance_reaction 
    ! hence, the following condition should not hold but we do it for completeness
    if (nreactions .lt. 1) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
       return
    end if  

    ! otherwise, complete the remaining part
    call build(bpt,"advance_reaction_SSA")

    dv = product(dx(1,1:dm))*cross_section


    !!!!! omp tiling

    ng_o = n_old(1)%ng
    ng_n = n_new(1)%ng

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(op,np,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi,n_old(n),tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          op => dataptr(n_old(n),i)
          np => dataptr(n_new(n),i)

          lo = lwb(get_box(n_new(n),i))
          hi = upb(get_box(n_new(n),i))
          select case (dm)
          case (2)
             call advance_reaction_SSA_2d(op(:,:,1,:),ng_o,np(:,:,1,:),ng_n, &
                                      lo,hi,tlo,thi,dv,dt)
          case (3)
             call advance_reaction_SSA_3d(op(:,:,:,:),ng_o,np(:,:,:,:),ng_n, &
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

    call destroy(bpt)

  end subroutine advance_reaction_SSA
  
  
  subroutine advance_reaction_SSA_2d(n_old,ng_o,n_new,ng_n,glo,ghi,tlo,thi,dv,dt)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_o, ng_n
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    
    ! local
    integer :: i, j

    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_SSA_cell(n_old(i,j,1:nspecies),n_new(i,j,1:nspecies), &
                                      dv,dt)
    end do
    end do

  end subroutine advance_reaction_SSA_2d


  subroutine advance_reaction_SSA_3d(n_old,ng_o,n_new,ng_n,glo,ghi,tlo,thi,dv,dt)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_o, ng_n
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,glo(3)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    
    ! local
    integer :: i, j, k

    do k=tlo(3),thi(3)
    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_SSA_cell(n_old(i,j,k,1:nspecies),n_new(i,j,k,1:nspecies), &
                                      dv,dt)
    end do
    end do
    end do

  end subroutine advance_reaction_SSA_3d


  subroutine advance_reaction_SSA_cell(n_old,n_new,dv,dt)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
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
       if (use_bl_rng) then
          rr = bl_rng_get(rng_dist_uniform_real_reaction,rng_eng_reaction)
       else
          call UniformRNG(rr)
       end if
       ! tau is how long until the next reaction occurs
       tau = -log(1-rr)/rTotal
       t_local = t_local + tau;

       if (t_local .gt. dt) exit EventLoop

       ! Select the next reaction according to relative rates
       if (use_bl_rng) then
          rr = bl_rng_get(rng_dist_uniform_real_reaction,rng_eng_reaction)
       else
          call UniformRNG(rr)
       end if
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
