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

  public :: advance_reaction

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! can delete this once reactions are working
    do n=1,nlevs
       call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  subroutine advance_reaction_2d(n_old,ng_o,n_new,ng_n,lo,hi,dv,dt)

    integer        , intent(in   ) :: lo(:),hi(:),ng_o,ng_n
    real(kind=dp_t), intent(in   ) :: n_old(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real(kind=dp_t), intent(in   ) :: dv,dt
    
    ! local
    integer :: i,j,comp,n

    real(kind=dp_t) :: reaction_rates(1:nreactions)
    real(kind=dp_t) :: avg_reactions(1:nreactions)

    integer :: num_reactions(1:nreactions)

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       
       if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
          ! first-order tau-leaping

          ! compute reaction rates in terms of reactions/volume
          call compute_reaction_rates(n_old(i,j,:),reaction_rates(:))

          ! compute mean number of events over the time step
          avg_reactions(:) = reaction_rates(:)*dt*dv

          do comp=1,nreactions

             ! for each reaction, compute how many reactions will happen
             ! by sampling a Poisson number
             if (reaction_rates(comp) .gt. 0.d0) then
                call PoissonNumber(num_reactions(comp),avg_reactions(comp))
             else
                num_reactions(comp) = 0
             end if

             ! update number densities for this reaction
             do n=1,nspecies
                n_new(i,j,n) = n_old(i,j,n) + &
                     num_reactions(comp)*stoichiometric_factors(n,comp) / dv
             end do

          end do

          if (reaction_type .eq. 1) then
             ! second-order tau-leaping corrector
             call bl_error("advance_reaction: reaction_type=1 not supported yet")
             
          end if

       else if (reaction_type .eq. 2 .or. reaction_type .eq. 3) then
          ! first-order CLE
          call bl_error("advance_reaction: reaction_type=2 not supported yet")

          if (reaction_type .eq. 3) then
             ! second-order CLE corrector
             call bl_error("advance_reaction: reaction_type=3 not supported yet")

          end if

       else if (reaction_type .eq. 4) then
          ! SSA
          call bl_error("advance_reaction: reaction_type=4 not supported yet")

       else
          call bl_error("advance_reaction: invalid reaction_type")
       end if

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

    real(kind=dp_t) :: reaction_rates(1:nreactions)
    real(kind=dp_t) :: avg_reactions(1:nreactions)

    integer :: num_reactions(1:nreactions)

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       
       if (reaction_type .eq. 0 .or. reaction_type .eq. 1) then
          ! first-order tau-leaping

          ! compute reaction rates
          reaction_rates(:) = 0.d0

          ! compute mean number of events
          avg_reactions(:) = reaction_rates(:)*dt

          do comp=1,nreactions

             ! for each reaction, compute how many reactions will happen
             ! by sampling a Poisson number
             if (reaction_rates(comp) .gt. 0.d0) then
                call PoissonNumber(num_reactions(comp),avg_reactions(comp))
             else
                num_reactions(comp) = 0
             end if

             ! update number densities for this reaction
             do n=1,nspecies
                n_new(i,j,k,n) = n_old(i,j,k,n) + &
                     num_reactions(comp)*stoichiometric_factors(n,comp) / dv
             end do

          end do

          if (reaction_type .eq. 1) then
             ! second-order tau-leaping corrector
             call bl_error("advance_reaction: reaction_type=1 not supported yet")
             
          end if

       else if (reaction_type .eq. 2 .or. reaction_type .eq. 3) then
          ! first-order CLE
          call bl_error("advance_reaction: reaction_type=2 not supported yet")

          if (reaction_type .eq. 3) then
             ! second-order CLE corrector
             call bl_error("advance_reaction: reaction_type=3 not supported yet")

          end if

       else if (reaction_type .eq. 4) then
          ! SSA
          call bl_error("advance_reaction: reaction_type=4 not supported yet")

       else
          call bl_error("advance_reaction: invalid reaction_type")
       end if

    end do
    end do
    end do

  end subroutine advance_reaction_3d

  subroutine compute_reaction_rates(n_in,reaction_rates)
    
    real(kind=dp_t), intent(in   ) :: n_in(:)
    real(kind=dp_t), intent(inout) :: reaction_rates(:)

    ! reaction 1: n1 + n2 -> n3
    reaction_rates(1) = 0.0001d0*n_in(1)*n_in(2)

    ! reaction 2: n3 -> n1 + n2
    reaction_rates(2) = 0.d0
    
  end subroutine compute_reaction_rates

end module advance_reaction_module
