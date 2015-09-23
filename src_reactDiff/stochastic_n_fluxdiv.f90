module stochastic_n_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_fill_random_module
  use multifab_physbc_module
  use average_to_faces_module
  use div_and_grad_module
  use probin_common_module, only: variance_coef_mass, initial_variance
  use probin_reactdiff_module, only: nspecies

  implicit none

  private

  public :: stochastic_n_fluxdiv, fill_mass_stochastic, &
       init_mass_stochastic, destroy_mass_stochastic, &
       add_n_fluctuations

  ! stochastic fluxes for mass densities are face-centered
  type(multifab), allocatable, save :: stoch_W_fc(:,:,:)

  integer, save :: n_rngs ! how many random number stages
  
contains
  
  subroutine stochastic_n_fluxdiv(mla,n_cc,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                  the_bc_tower,increment_in)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: n_cc(:)
    type(multifab) , intent(in   )  :: diff_coef_face(:,:)
    type(multifab) , intent(inout)  :: stoch_fluxdiv(:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    real(kind=dp_t), intent(in   )  :: dt
    type(bc_tower) , intent(in   )  :: the_bc_tower
    logical  , intent(in), optional :: increment_in

    integer :: i,dm,n,nlevs
    real(kind=dp_t)  :: variance

    type(multifab) :: flux(mla%nlevel,mla%dim)

    logical :: increment_div

    type(bl_prof_timer), save :: bpt

    dm = mla%dim
    nlevs = mla%nlevel

    call build(bpt,"stochastic_n_fluxdiv")

    increment_div = .false.
    if (present(increment_in)) increment_div = increment_in

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! average n to faces, store in "flux" so as to avoid an extra multifab
    ! alternatively one could multiply diff_coeff_face with this number
    call average_to_faces(mla,n_cc,flux,1,1,nspecies)

    ! assumble fluxes on faces, sqrt(2*D_k*n_k / (dt*dV)) * random_normal
    variance = sqrt(2.d0*variance_coef_mass/(product(dx(1,1:dm))*dt))        
    call assemble_stoch_n_fluxes(mla,n_cc,diff_coef_face,flux)
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s(flux(n,i), variance, 0)
       end do
    end do      

    ! take flux divergence
    call compute_div(mla,flux,stoch_fluxdiv,dx,1,1,nspecies,increment_div)

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine stochastic_n_fluxdiv

  subroutine assemble_stoch_n_fluxes(mla,n_cc,diff_coef_face,flux)

    ! note: n averaged to faces is stored in "flux" on entry to this subroutine

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: n_cc(:)
    type(multifab) , intent(in   )  :: diff_coef_face(:,:)
    type(multifab) , intent(inout)  :: flux(:,:) ! On input this contains the face-averaged number densities

    integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
    integer :: ng_n,ng_d,ng_f,ng_s

    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: dx(:,:,:,:)
    real(kind=dp_t), pointer :: dy(:,:,:,:)
    real(kind=dp_t), pointer :: dz(:,:,:,:)
    real(kind=dp_t), pointer :: fx(:,:,:,:)
    real(kind=dp_t), pointer :: fy(:,:,:,:)
    real(kind=dp_t), pointer :: fz(:,:,:,:)
    real(kind=dp_t), pointer :: sx(:,:,:,:)
    real(kind=dp_t), pointer :: sy(:,:,:,:)
    real(kind=dp_t), pointer :: sz(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_cc(1)%ng
    ng_d = diff_coef_face(1,1)%ng
    ng_f = flux(1,1)%ng
    ng_s = stoch_W_fc(1,1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(n_cc(n))
          np => dataptr(n_cc(n),i)
          dx => dataptr(diff_coef_face(n,1),i)
          dy => dataptr(diff_coef_face(n,2),i)
          fx => dataptr(flux(n,1),i)
          fy => dataptr(flux(n,2),i)
          sx => dataptr(stoch_W_fc(n,1,1),i)
          sy => dataptr(stoch_W_fc(n,2,1),i)
          lo = lwb(get_box(n_cc(n),i))
          hi = upb(get_box(n_cc(n),i))
          select case (dm)
          case (2)
             call assemble_stoch_n_fluxes_2d(np(:,:,1,:),ng_n, &
                                             dx(:,:,1,:),dy(:,:,1,:),ng_d, &
                                             fx(:,:,1,:),fy(:,:,1,:),ng_f, &
                                             sx(:,:,1,:),sy(:,:,1,:),ng_s, lo,hi)
          case (3)
             dz => dataptr(diff_coef_face(n,3),i)
             fz => dataptr(flux(n,3),i)
             sz => dataptr(stoch_W_fc(n,3,1),i)
             call assemble_stoch_n_fluxes_3d(np(:,:,:,:),ng_n, &
                                             dx(:,:,:,:),dy(:,:,:,:),dz(:,:,:,:),ng_d, &
                                             fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),ng_f, &
                                             sx(:,:,:,:),sy(:,:,:,:),sz(:,:,:,:),ng_s, lo,hi)
          end select
       end do
    end do

    ! sync the fluxes at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_internal_sync(flux(n,i))
          call multifab_fill_boundary(flux(n,i))  
       end do
    end do

  end subroutine assemble_stoch_n_fluxes

  subroutine assemble_stoch_n_fluxes_2d(n_cc,ng_n,coefx,coefy,ng_d,fluxx,fluxy,ng_f, &
                                        stochx,stochy,ng_s,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_d,ng_f,ng_s
    real(kind=dp_t), intent(in   ) ::   n_cc(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t), intent(in   ) ::  coefx(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefy(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: stochx(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochy(lo(1)-ng_s:,lo(2)-ng_s:,:)

    integer :: i,j

    ! x-fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fluxx(i,j,1:nspecies) = &
               sqrt(coefx(i,j,1:nspecies)*fluxx(i,j,1:nspecies))*stochx(i,j,1:nspecies)
       end do
    end do

    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fluxy(i,j,1:nspecies) = &
               sqrt(coefy(i,j,1:nspecies)*fluxy(i,j,1:nspecies))*stochy(i,j,1:nspecies)
       end do
    end do

  end subroutine assemble_stoch_n_fluxes_2d

  subroutine assemble_stoch_n_fluxes_3d(n_cc,ng_n,coefx,coefy,coefz,ng_d,fluxx,fluxy,fluxz,ng_f, &
                                        stochx,stochy,stochz,ng_s,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_d,ng_f,ng_s
    real(kind=dp_t), intent(in   ) ::   n_cc(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real(kind=dp_t), intent(in   ) ::  coefx(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefy(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefz(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: stochx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    integer :: i,j,k

    ! x-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fluxx(i,j,k,1:nspecies) = &
                  sqrt(coefx(i,j,k,1:nspecies)*fluxx(i,j,k,1:nspecies))*stochx(i,j,k,1:nspecies)
          end do
       end do
    end do

    ! y-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fluxy(i,j,k,1:nspecies) = &
                  sqrt(coefy(i,j,k,1:nspecies)*fluxy(i,j,k,1:nspecies))*stochy(i,j,k,1:nspecies)
          end do
       end do
    end do

    ! z-fluxes
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fluxz(i,j,k,1:nspecies) = &
                  sqrt(coefz(i,j,k,1:nspecies)*fluxz(i,j,k,1:nspecies))*stochz(i,j,k,1:nspecies)
          end do
       end do
    end do

  end subroutine assemble_stoch_n_fluxes_3d

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_mass_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,i,dm,comp
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_mass_stochastic")

    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(stoch_W_fc(mla%nlevel, mla%dim, n_rngs))

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             ! we need one face-centered flux for each concentration
             call multifab_build_edge(stoch_W_fc(n,i,comp),mla%la(n),nspecies,0,i)
          end do
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

    call destroy(bpt)

  end subroutine init_mass_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_mass_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm,comp
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"destroy_mass_stochastic")

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             call multifab_destroy(stoch_W_fc(n,i,comp))
          end do
       end do
    end do
    
    deallocate(stoch_W_fc)

    call destroy(bpt)

  end subroutine destroy_mass_stochastic

  subroutine fill_mass_stochastic(mla,the_bc_level)
  
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! Local variables
    integer :: dm,nlevs,box,i,rng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_mass_stochastic")

    nlevs = mla%nlevel
    dm    = mla%dim    
    
    ! generate and store the stochastic flux (random numbers)
    do rng=1, n_rngs
       do i = 1,dm
          call multifab_fill_random(stoch_W_fc(:,i,rng))
       end do
    end do

    ! apply boundary conditions to stochastic fluxes
    call stoch_mass_bc(mla,the_bc_level)
  
    call destroy(bpt)

  end subroutine fill_mass_stochastic

  subroutine stoch_mass_bc(mla,the_bc_level)
    
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local
    integer :: n,nlevs,dm,idim,comp,i,ng_f
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stoch_mass_bc")

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_f = stoch_W_fc(1,1,1)%ng

    do n=1,nlevs
       do idim=1,dm
          do comp=1,n_rngs
             do i=1,nfabs(stoch_W_fc(n,idim,comp))
                fp => dataptr(stoch_W_fc(n,idim,comp),i)
                lo = lwb(get_box(stoch_W_fc(n,idim,comp),i))
                hi = upb(get_box(stoch_W_fc(n,idim,comp),i))
                select case (dm)
                case (2)
                   call stoch_mass_bc_2d(fp(:,:,1,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                case (3)
                   call stoch_mass_bc_3d(fp(:,:,:,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                end select
             end do
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine stoch_mass_bc

  subroutine stoch_mass_bc_2d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),:)
         end if

      else

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_2d

  subroutine stoch_mass_bc_3d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
         end if

      else if (idim .eq. 2) then

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:)
         end if

      else

         if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = 0.d0
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = 0.d0
         end if

         if (phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:)
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_3d

  subroutine add_n_fluctuations(mla,n_init,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,dm,comp,n_cell
    real(kind=dp_t) :: dn_sum

    type(multifab) :: n_temp(mla%nlevel)

    nlevs = mla%nlevel
    dm = mla%dim

    ! the number of ghost cells must match variance_mfab input to multifab_fill_random
    ! the values in the ghost cells do not have to be added to n_init since we
    ! fill ghost cells for n_init afterwards
    do n=1,nlevs
       call multifab_build(n_temp(n),mla%la(n),nspecies,n_init(n)%ng)
    end do

    n_cell = multifab_volume(n_temp(1)) / nspecies

    ! create a multifab full of random numbers
    do n=1,nlevs
       call multifab_fill_random(n_temp(n:n), &
                                 variance_mfab=n_init, &
                                 variance=initial_variance*variance_coef_mass/product(dx(n,1:dm)))

       ! Make sure this sums to zero
       do comp=1, nspecies
          dn_sum = multifab_sum_c(n_temp(n),comp,1) / dble(n_cell)
          call multifab_sub_sub_s_c(n_temp(n),comp,dn_sum,1,0)
       end do   
          
    end do

    do n=1,nlevs
       call multifab_plus_plus_c(n_init(n),1,n_temp(n),1,nspecies,0)
       call multifab_fill_boundary(n_init(n))
       call multifab_physbc(n_init(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    do n=1,nlevs
       call multifab_destroy(n_temp(n))
    end do    

  end subroutine add_n_fluctuations
   
end module stochastic_n_fluxdiv_module
