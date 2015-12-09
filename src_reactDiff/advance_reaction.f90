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

  ! Here we use Mattingly's predictor-corrector with theta=0.5d0, giving these parameters
  ! With these parameters this is actually equivalent to a traditional midpoint scheme
  real(kind=dp_t), parameter :: theta = 0.5d0
  real(kind=dp_t), parameter :: alpha1 = 2.0d0
  real(kind=dp_t), parameter :: alpha2 = 1.0d0

contains

  ! This solves dn/dt = f(n) - g (note the minus sign for g)
  ! where f(n) are the chemical production rates (deterministic or stochastic)
  ! and g=ext_src_in is an optional, constant (in time) *deterministic* source term.
  ! To model stochastic particle production (sources) include g in the definition of f instead.
  ! If return_rates_in=0 (or not present), the code returns time-advanced n in the 'n_new' field.
  ! If return_rates_in=1, the code returns the deterministic rates.
  ! If return_rates_in=2, the code returns the stochastic rates sampled by using random numbers specified by use_Poisson_rng.
  ! Note that chemical *production* rate (per species) is
  ! the sum over reactions of chemical *reaction* rates (per reaction) times the stochiometric coefficients
  ! The optional argument n_pred_mattingly is for the second step of (unsplitting) Mattingly's scheme.
  ! For return_rates_in=2, if n_pred_mattingly is present, rates are sampled from the average reaction rates 2*a(n_pred_mattingly)-a(n_old).
  ! (cf. otherwisze, they are sampled from a(n_old).)

  subroutine advance_reaction(mla,n_old,n_new,dx,dt,the_bc_tower,return_rates_in,ext_src_in,n_pred_mattingly)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:) ! if return_rates_in=0 or not present, return new state
                                               ! if return_rates_in=1 or 2, return rates
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer, intent(in), optional  :: return_rates_in
    type(multifab) , intent(in   ), optional :: ext_src_in(:)
    type(multifab) , intent(in   ), optional :: n_pred_mattingly(:)

    ! local
    integer :: n,nlevs,dm,i,ng_o,ng_n,ng_e,ng_p
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t) :: dv

    real(kind=dp_t), pointer :: op(:,:,:,:) ! for n_old
    real(kind=dp_t), pointer :: np(:,:,:,:) ! for n_new
    real(kind=dp_t), pointer :: ep(:,:,:,:) ! for ext_src
    real(kind=dp_t), pointer :: pp(:,:,:,:) ! for n_pred_mattingly

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(multifab) :: ext_src(mla%nlevel)

    type(bl_prof_timer),save :: bpt
    
    integer :: return_rates ! Should we do an actual time step or compute the chemical rates?
    
    logical :: is_mattingly_2nd_step ! true only if return_rates_in=2 and present(n_pred_mattingly)

    nlevs = mla%nlevel
    dm = mla%dim
    
    return_rates=0
    if(present(return_rates_in)) return_rates=return_rates_in

    is_mattingly_2nd_step = .false.
    if(present(n_pred_mattingly)) then
       if(return_rates .eq. 2) then
          is_mattingly_2nd_step = .true.
       else
          call bl_error("advance_reaction: n_pred_mattingly can be given only with return_rates=2")
       end if
    end if

    ! check
    if(present(ext_src_in) .and. (return_rates .eq. 1 .or. return_rates .eq. 2)) then
       call bl_error("advance_reaction: return_rates>=1 with an external source not implemented")
    end if

    if(return_rates .lt. 0 .or. return_rates .gt. 2) then
       call bl_error("advance_reaction: invalid return_rates_in")
    end if

    ! There are no reactions to process!
    if(nreactions<1) then
       do n=1,nlevs
          if(return_rates .eq. 0) then
             call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng) ! make sure n_new contains the new state
          else   
             call setval(n_new(n),0.d0,all=.true.)
          end if   
       end do
       return
    end if   

    call build(bpt,"advance_reaction")

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

    ng_o = n_old(1)%ng
    ng_n = n_new(1)%ng
    ng_e = ext_src(1)%ng

    if (is_mattingly_2nd_step) then
       ng_p = n_pred_mattingly(1)%ng
    else
       ng_p = n_old(1)%ng ! won't be used
    end if

    dv = product(dx(1,1:dm))
    if (dm<3) dv = dv*cross_section

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(op,np,ep,pp,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi, n_old(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          op => dataptr(n_old(n),i)
          np => dataptr(n_new(n),i)
          ep => dataptr(ext_src(n),i)

          if (is_mattingly_2nd_step) then
             pp => dataptr(n_pred_mattingly(n),i)
          else
             pp => dataptr(n_old(n),i) ! won't be used
          end if

          lo = lwb(get_box(n_new(n),i))
          hi = upb(get_box(n_new(n),i))
          select case (dm)
          case (2)
             call advance_reaction_2d(op(:,:,1,:),ng_o,np(:,:,1,:),ng_n,ep(:,:,1,:),ng_e, &
                                      lo,hi,tlo,thi,dv,dt,return_rates,                   &
                                      is_mattingly_2nd_step,pp(:,:,1,:),ng_p)
          case (3)
             call advance_reaction_3d(op(:,:,:,:),ng_o,np(:,:,:,:),ng_n,ep(:,:,:,:),ng_e, &
                                      lo,hi,tlo,thi,dv,dt,return_rates,                   &
                                      is_mattingly_2nd_step,pp(:,:,:,:),ng_p)
          end select
       end do
    end do
    !$omp end parallel

    do n=1,nlevs
       call multifab_destroy(ext_src(n))
    end do

    if(return_rates .eq. 1 .or. return_rates .eq. 2) return ! We are done
    
    ! Ensure ghost cells are consistent or n_new
    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    call destroy(bpt)

  end subroutine advance_reaction
  
  
  subroutine advance_reaction_2d(n_old,ng_o,n_new,ng_n,ext_src,ng_e,glo,ghi,tlo,thi,dv,dt,return_rates, &
                                 is_mattingly_2nd_step,n_pred_mattingly,ng_p)

    integer        , intent(in   ) :: glo(:),ghi(:),tlo(:),thi(:),ng_o,ng_n,ng_e
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: ext_src(glo(1)-ng_e:,glo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    integer,         intent(in   ) :: return_rates
    logical,         intent(in   ) :: is_mattingly_2nd_step
    integer        , intent(in   ) :: ng_p
    real(kind=dp_t), intent(in   ) :: n_pred_mattingly(glo(1)-ng_p:,glo(2)-ng_p:,:)
    
    ! local
    integer :: i,j

    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_cell(n_old(i,j,1:nspecies), n_new(i,j,1:nspecies),           &
                                  ext_src(i,j,1:nspecies), dv, dt, return_rates,          &
                                  is_mattingly_2nd_step,n_pred_mattingly(i,j,1:nspecies))
    end do
    end do

  end subroutine advance_reaction_2d

  subroutine advance_reaction_3d(n_old,ng_o,n_new,ng_n,ext_src,ng_e,glo,ghi,tlo,thi,dv,dt,return_rates, &
                                 is_mattingly_2nd_step,n_pred_mattingly,ng_p)

    integer        , intent(in   ) :: glo(:),ghi(:),tlo(:),thi(:),ng_o,ng_n,ng_e
    real(kind=dp_t), intent(in   ) ::   n_old(glo(1)-ng_o:,glo(2)-ng_o:,glo(3)-ng_o:,:)
    real(kind=dp_t), intent(inout) ::   n_new(glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: ext_src(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    integer,         intent(in   ) :: return_rates
    logical,         intent(in   ) :: is_mattingly_2nd_step
    integer,         intent(in   ) :: ng_p
    real(kind=dp_t), intent(in   ) :: n_pred_mattingly(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:,:)
    
    ! local
    integer :: i,j,k

    do k=tlo(3),thi(3)
    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
       call advance_reaction_cell(n_old(i,j,k,1:nspecies), n_new(i,j,k,1:nspecies),         &
                                  ext_src(i,j,k,1:nspecies), dv, dt, return_rates,          &
                                  is_mattingly_2nd_step,n_pred_mattingly(i,j,k,1:nspecies))
    end do
    end do
    end do

  end subroutine advance_reaction_3d

  subroutine advance_reaction_cell(n_old,n_new,ext_src,dv,dt,return_rates, &
                                   is_mattingly_2nd_step,n_pred_mattingly)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: ext_src(:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    integer,         intent(in   ) :: return_rates
    logical,         intent(in   ) :: is_mattingly_2nd_step
    real(kind=dp_t), intent(in   ) :: n_pred_mattingly(:)

    real(kind=dp_t) :: avg_reactions     (1:nreactions)
    real(kind=dp_t) :: avg_reactions_pred(1:nreactions)
    real(kind=dp_t) :: num_reactions     (1:nreactions)
    real(kind=dp_t) :: rTotal, rr, rSum, tau, t_local

    integer :: spec,reaction,which_reaction
    integer :: n_steps_SSA
    
    !!!!!!!!!!!!!!!!!!!
    ! returning rates !
    !!!!!!!!!!!!!!!!!!!

    ! compute reaction rates only in units of (number density) / time

    if(return_rates .eq. 1) then  ! deterministic rates case
       call compute_reaction_rates(n_old(1:nspecies), avg_reactions, dv)

       n_new = 0.d0
       do reaction=1,nreactions
          n_new(1:nspecies) = n_new(1:nspecies) + avg_reactions(reaction) * &
             (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))
       end do

       return

    else if(return_rates .eq. 2) then  ! stochastic rates case

       if(is_mattingly_2nd_step) then
          ! calculate 2*a(n_pred)-a(n_old)
          ! note that the use of avg_reactions_pred is different from the time advancing case below.
          ! here, avg_reactions_pred saves a(n_pred), whereas it saves a(n_old) below. 
          call compute_reaction_rates(n_old(1:nspecies), avg_reactions, dv)
          call compute_reaction_rates(n_pred_mattingly(1:nspecies), avg_reactions_pred, dv)
          avg_reactions = 2.d0*avg_reactions_pred - avg_reactions
          ! As impletemented below for second-order tau-leaping, we could count n_rejections here.
       else
          ! calculate a(n_old)
          call compute_reaction_rates(n_old(1:nspecies), avg_reactions, dv)
       end if

       ! compute mean number of events over the time step
       avg_reactions = max(0.0d0, avg_reactions*dt*dv)

       ! sampling
       n_new = 0.d0
       do reaction=1,nreactions
          ! compute num_reactions(reaction)
          call sample_num_reactions(reaction)
          ! after loop, n_new will have stochastic rates (in units of (number density) / time)  
          n_new(1:nspecies) = n_new(1:nspecies) + num_reactions(reaction)/dv/dt * &
             (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))
       end do

       return

    end if

    !!!!!!!!!!!!!!!!!!
    ! advancing time !
    !!!!!!!!!!!!!!!!!!
 
    ! copy old state into new
    n_new = n_old

    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       ! first-order tau-leaping or CLE

       ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
       call compute_reaction_rates(n_new(1:nspecies), avg_reactions, dv)
       !write(*,*) "PREDICTOR PROPENSITY=", real(avg_reactions); stop
       ! save the mean reactions from the predictor
       avg_reactions_pred = avg_reactions

       ! compute mean number of events over the time step
       avg_reactions = max(0.0d0, avg_reactions*dt*dv)

       ! if reaction_type=1, we are doing a predictor/corrector, so scale the # of
       ! reactions by the predictor length (for midpoint, this is half the time interval)
       if (reaction_type .eq. 1) avg_reactions = avg_reactions*theta

       do reaction=1,nreactions
          
          ! compute num_reactions(reaction)
          call sample_num_reactions(reaction)

          ! update number densities for this reaction
          n_new(1:nspecies) = n_new(1:nspecies) + num_reactions(reaction)/dv * &
             (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))

       end do

       ! add contribution from external source
       if (reaction_type .eq. 1) then
          n_new(1:nspecies) = n_new(1:nspecies) - dt*theta*ext_src(1:nspecies)
       else
          n_new(1:nspecies) = n_new(1:nspecies) - dt*ext_src(1:nspecies)
       end if

       if (reaction_type .eq. 1) then
          ! second-order tau-leaping or CLE corrector
          ! Mattingly predictor-corrector with theta=0.5d0

          ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
          call compute_reaction_rates(n_new(1:nspecies), avg_reactions, dv)
          !write(*,*) "CORRECTOR PROPENSITY=", real(avg_reactions); stop
          
          ! Corrector rate is a linear combination of the two rates
          ! This is equivalent to a midpoint scheme
          avg_reactions = (alpha1*avg_reactions-alpha2*avg_reactions_pred)*(1.d0-theta)
          
          ! compute mean number of events over the time step
          avg_reactions = avg_reactions*dt*dv

          do reaction=1,nreactions

             if (avg_reactions(reaction) .lt. 0.d0) then
                n_rejections = n_rejections + 1
                avg_reactions(reaction) = 0.d0
             end if

             ! compute num_reactions(reaction)
             call sample_num_reactions(reaction)

             ! update number densities for this reaction
             n_new(1:nspecies) = n_new(1:nspecies) + num_reactions(reaction)/dv * &
                  (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))

          end do

          ! add contribution from external source
          n_new(1:nspecies) = n_new(1:nspecies) - dt*(alpha1-alpha2)*(1.d0-theta)*ext_src(1:nspecies)

       end if
       
    else if (reaction_type .eq. 2) then ! SSA

       t_local = 0.d0
       n_steps_SSA = 0

       EventLoop: do

          ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
          call compute_reaction_rates(n_new(1:nspecies), avg_reactions, dv)

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
             n_new(spec) = n_new(spec) + &
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
           case(-1)
               ! Do deterministic chemistry   
              num_reactions(comp) = avg_reactions(comp)
           case default    
              call bl_error("advance_reaction: invalid use_Poisson_rng")
           end select
        else
           num_reactions(comp) = 0
        end if

      end subroutine sample_num_reactions

  end subroutine advance_reaction_cell

end module advance_reaction_module
