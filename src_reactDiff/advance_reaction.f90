module advance_reaction_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use compute_reaction_rates_module
  use probin_common_module, only: seed
  use probin_reactdiff_module, only: nspecies, nreactions, reaction_type, &
       stoichiometric_factors, use_Poisson_rng, cross_section

  implicit none

  private

  public :: advance_reaction, n_rejections
  
  ! keep track of the number of corrector steps that are null (rejected) in the second-order method
  integer*8 ::  n_rejections=0 ! Number of reaction rates set to zero in corrector stage

  ! Here we use Mattingly's predictor-corrector with theta=0.5d0, giving the parameters:
  real(kind=dp_t), parameter :: theta = 0.5d0
  real(kind=dp_t), parameter :: alpha1 = 2.0d0
  real(kind=dp_t), parameter :: alpha2 = 1.0d0

contains

  subroutine advance_reaction(mla,n_old,n_new,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,dm,i,ng_o,ng_n,lo(mla%dim),hi(mla%dim)
    real(kind=dp_t) :: dv

    real(kind=dp_t), pointer :: op(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer),save :: bpt
    
    nlevs = mla%nlevel
    dm = mla%dim

    ! There are no reactions to process!
    if(nreactions<1) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
       end do
       return
    end if   

    call build(bpt,"advance_reaction")

    ng_o = n_old(1)%ng
    ng_n = n_new(1)%ng

    dv = product(dx(1,1:dm))
    if (dm<3) dv = dv*cross_section

    ! if there are no reactions, copy old state into new and return
    if (nreactions .eq. 0) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
       call destroy(bpt)
       return
    end if

    do n=1,nlevs
       do i=1,nfabs(n_new(n))
          op => dataptr(n_new(n),i)
          np => dataptr(n_new(n),i)
          lo = lwb(get_box(n_new(n),i))
          hi = upb(get_box(n_new(n),i))
          select case (dm)
          case (2)
             call advance_reaction_2d(op(:,:,1,:),ng_o,np(:,:,1,:),ng_n,lo,hi,dv,dt)
          case (3)
             call advance_reaction_3d(op(:,:,:,:),ng_o,np(:,:,:,:),ng_n,lo,hi,dv,dt)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    call destroy(bpt)

  end subroutine advance_reaction
  
  
  subroutine advance_reaction_2d(n_old,ng_o,n_new,ng_n,lo,hi,dv,dt)

    integer        , intent(in   ) :: lo(:),hi(:),ng_o,ng_n
    real(kind=dp_t), intent(in   ) :: n_old(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    
    ! local
    integer :: i,j

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)       
       call advance_reaction_cell(n_old(i,j,1:nspecies), n_new(i,j,1:nspecies), dv, dt)
    end do
    end do

  end subroutine advance_reaction_2d

  subroutine advance_reaction_3d(n_old,ng_o,n_new,ng_n,lo,hi,dv,dt)

    integer        , intent(in   ) :: lo(:),hi(:),ng_o,ng_n
    real(kind=dp_t), intent(in   ) :: n_old(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:,:)
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_n:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    
    ! local
    integer :: i,j,k

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       call advance_reaction_cell(n_old(i,j,k,1:nspecies), n_new(i,j,k,1:nspecies), dv, dt)
    end do
    end do
    end do

  end subroutine advance_reaction_3d

  subroutine advance_reaction_cell(n_old,n_new,dv,dt)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dv,dt

    real(kind=dp_t) :: avg_reactions     (1:nreactions)
    real(kind=dp_t) :: avg_reactions_pred(1:nreactions)
    real(kind=dp_t) :: num_reactions     (1:nreactions)
    real(kind=dp_t) :: rTotal, rr, rSum, tau, t_local

    integer :: spec,reaction,which_reaction
    integer :: n_steps_SSA

    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       ! first-order tau-leaping or CLE

       ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
       call compute_reaction_rates(n_old(1:nspecies), avg_reactions, dv)
       !write(*,*) "PREDICTOR PROPENSITY=", real(avg_reactions)

       ! compute mean number of events over the time step
       if (reaction_type .eq. 1) avg_reactions = avg_reactions*theta ! Predictor step has length theta*dt
       avg_reactions = max(0.0d0, avg_reactions*dt*dv)

       do reaction=1,nreactions
          
          ! compute num_reactions(reaction)
          call sample_num_reactions(reaction)

          ! update number densities for this reaction
          n_new(1:nspecies) = n_old(1:nspecies) + num_reactions(reaction)/dv * &
             (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))

       end do

       if (reaction_type .eq. 1) then
          ! second-order tau-leaping or CLE corrector
          ! Mattingly predictor-corrector with theta=0.5d0

          ! save the mean reactions from the predictor
          avg_reactions_pred = avg_reactions

          ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
          call compute_reaction_rates(n_new(1:nspecies),avg_reactions, dv)
          !write(*,*) "CORRECTOR PROPENSITY=", real(avg_reactions); stop

          ! compute mean number of events over the time step
          avg_reactions = avg_reactions*dt*dv

          avg_reactions = (alpha1*avg_reactions_pred-alpha2*avg_reactions)*(1.d0-theta)

          do reaction=1,nreactions

             if (avg_reactions(reaction) .lt. 0.d0) then
                n_rejections = n_rejections + 1
                avg_reactions(reaction) = 0.d0
             end if

             ! compute num_reactions(reaction)
             call sample_num_reactions(reaction)

             ! update number densities for this reaction
             do spec=1,nspecies
                n_new(spec) = n_new(spec) + num_reactions(reaction)/dv * &
                  (stoichiometric_factors(spec,2,reaction)-stoichiometric_factors(spec,1,reaction))
             end do

          end do

       end if
       
    else if (reaction_type .eq. 2) then ! SSA

       t_local = 0.d0
       n_steps_SSA = 0

       EventLoop: do

          ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
          call compute_reaction_rates(n_old(1:nspecies), avg_reactions, dv)

          ! compute reaction rates in units (# reactions) / (unit time)
          avg_reactions = max(0.0d0, avg_reactions*dv)

          ! sum the reaction rates
          rTotal = sum(avg_reactions(1:nreactions))

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
             rSum = rSum + avg_reactions(reaction)
             which_reaction = reaction
             if( rSum >= rr ) exit FindReaction
          end do FindReaction

          ! update number densities for this reaction
          do spec=1,nspecies
             n_new(spec) = n_old(spec) + &
                  (stoichiometric_factors(spec,2,which_reaction)-stoichiometric_factors(spec,1,which_reaction)) / dv
          end do
          
          n_steps_SSA = n_steps_SSA+1

       end do EventLoop

    else
       call bl_error("advance_reaction: invalid reaction_type")
    end if

    contains

      ! compute num_reactions(:)
      subroutine sample_num_reactions(comp) ! Auxilliary routine (should be inlined by compiler)

        integer, intent(in) :: comp

        ! local
        integer :: tmp

        ! for each reaction, compute how many reactions will happen
        ! by sampling a Poisson (tau leaping) or Gaussian (CLE) number
        if (avg_reactions(comp) .gt. 0.d0) then
           select case(use_Poisson_rng)           
           case(1)
              ! Need a Poisson random number for tau leaping
              call PoissonRNG(number=tmp, mean=avg_reactions(comp))
              num_reactions(comp) = tmp ! Convert to real
           case(0)
              ! Need a Gaussian random number for CLE
              call NormalRNG(num_reactions(comp))
              num_reactions(comp) = avg_reactions(comp) + sqrt(avg_reactions(comp))*num_reactions(comp)
           case default ! Do deterministic chemistry   
              num_reactions(comp) = avg_reactions(comp)
           end select
        else
           num_reactions(comp) = 0
        end if

      end subroutine sample_num_reactions

  end subroutine advance_reaction_cell

end module advance_reaction_module
