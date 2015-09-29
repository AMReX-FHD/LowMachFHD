module advance_reaction_module

  use ml_layout_module
  use bc_module
  use define_bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use probin_common_module, only: seed
  use probin_reactdiff_module, only: nspecies, nreactions, reaction_type, &
       stoichiometric_factors

  implicit none

  private

  public :: advance_reaction, n_rejections
  
  ! Donev: It is important to keep track of the number of corrector steps that are null (rejected) in the second-order method
  ! This should be printed in the end (perhaps as a fraction of what number of steps were rejected)
  integer*8 ::  n_rejections=0 ! Number of reaction rates set to zero in corrector stage

  ! Donev: Made these compile-time constants since they are so simple
  ! Here we use Mattingly's predictor-corrector with theta=0.5d0, giving the parameters:
  real(kind=dp_t), parameter :: theta=0.5d0
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

    nlevs = mla%nlevel
    dm = mla%dim

    ng_o = n_old(1)%ng
    ng_n = n_new(1)%ng

    dv = product(dx(1,1:dm))

    ! if there are no reactions, copy old state into new and return
    if (nreactions .eq. 0) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
       end do
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

  end subroutine advance_reaction
  
  ! Donev: I extracted this local routine here outside of the loops
  ! This way code is not duplicated twice and it is clear that the code is purely local
  subroutine advance_reaction_cell(n_old,n_new,dv,dt)
    real(kind=dp_t), intent(in   ) :: n_old(n_species)
    real(kind=dp_t), intent(inout) :: n_new(n_species)
    real(kind=dp_t), intent(in   ) :: dv,dt

    real(kind=dp_t) :: avg_reactions     (1:nreactions)
    real(kind=dp_t) :: avg_reactions_pred(1:nreactions)

    integer :: num_reactions(1:nreactions)

    if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
       ! first-order tau-leaping

       ! compute reaction rates in terms of (reaction rate)/volume
       ! Donev: dv needs to be an argument to compute_reaction_rates
       ! This because things like n^2 need to be replaced by N*(N-1)/dv^2, where N=n*dV is number of molecules
       call compute_reaction_rates(n_old(i,j,:), avg_reactions)

       ! compute mean number of events over the time step
       avg_reactions = max(0.0_dp_t, avg_reactions*dt*dv) ! Donev: Moved test for negativity here

       do comp=1,nreactions
          call sample_num_reactions(comp)

          ! update number densities for this reaction
          ! Donev: Rewrote this to use array syntax since it reads nicer I think
          n_new(1:nspecies) = n_old(1:nspecies) + num_reactions(comp)/dv * &
             (stoichiometric_factors(1:nspecies,2,comp)-stoichiometric_factors(1:nspecies,1,comp))

       end do

       if (reaction_type .eq. 1) then
          ! second-order tau-leaping corrector
          ! Mattingly predictor-corrector with theta=0.5d0

          ! save the mean reactions from the predictor
          avg_reactions_pred = avg_reactions

          ! compute reaction rates in terms of (reaction rate)/volume
          call compute_reaction_rates(n_new(i,j,:),avg_reactions)

          ! compute mean number of events over the time step
          avg_reactions = avg_reactions*dt*dv

          avg_reactions = (alpha1*avg_reactions_pred-alpha2*avg_reactions)*(1.d0-theta)

          do comp=1,nreactions

             if (avg_reactions(comp) .lt. 0.d0) then
                n_rejections = n_rejections + 1
                avg_reactions(comp) = 0.d0
             end if

             call sample_num_reactions(comp)

             ! update number densities for this reaction
             do n=1,nspecies
                n_new(n) = n_old(n) + num_reactions(comp)/dv * &
                  (stoichiometric_factors(1:nspecies,2,comp)-stoichiometric_factors(1:nspecies,1,comp))
             end do

          end do

       end if
       
    else if (reaction_type .eq. 2) then ! SSA
    
       call bl_error("advance_reaction: reaction_type=2 (SSA) not supported yet")

    else
       call bl_error("advance_reaction: invalid reaction_type")
    end if

  contains
    
    subroutine sample_num_reactions(comp) ! Auxilliary routine (should be inlined by compiler)
      integer, intent(in) :: comp
      ! for each reaction, compute how many reactions will happen
      ! by sampling a Poisson or Gaussian number
      if (avg_reactions(comp) .gt. 0.d0) then
         ! Donev: Either do tau leaping or CLE:
         if(use_Poisson_rng) then
             call PoissonNumber(number=num_reactions(comp), mean=avg_reactions(comp))
         else
             call NormalRNG(num_reactions(comp))
             num_reactions(comp) = avg_reactions(comp) + sqrt(avg_reactions(comp))*num_reactions(comp)
         end if    
      else
         num_reactions(comp) = 0
      end if
    end subroutine           
    
  end subroutine
  
  subroutine advance_reaction_2d(n_old,ng_o,n_new,ng_n,lo,hi,dv,dt)

    integer        , intent(in   ) :: lo(:),hi(:),ng_o,ng_n
    real(kind=dp_t), intent(in   ) :: n_old(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    
    ! local
    integer :: i,j,comp,n

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
    integer :: i,j,k,comp,n

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       call advance_reaction_cell(n_old(i,j,k,1:nspecies), n_new(i,j,k,1:nspecies), dv, dt)
    end do
    end do
    end do

  end subroutine advance_reaction_3d

  ! Donev: This routine should be something in a separate .f90 file that can be overwritten in an exec module
  ! That is, instead of listing all sorts of options in the code, each run can override this routine if needed
  ! We will provide a default routine based on the Law of Mass Action (LMA) here as a default but applications can change the chemistry laws
  ! Donev: Will write the default routine later
  subroutine compute_reaction_rates(n_in,reaction_rates)
    
    real(kind=dp_t), intent(in   ) :: n_in(:)
    real(kind=dp_t), intent(inout) :: reaction_rates(:)

    ! reaction 1: n1 + n2 -> n3
    reaction_rates(1) = 0.0001d0*n_in(1)*n_in(2)

    ! reaction 2: n3 -> n1 + n2
    reaction_rates(2) = 0.d0
    
  end subroutine compute_reaction_rates

end module advance_reaction_module
