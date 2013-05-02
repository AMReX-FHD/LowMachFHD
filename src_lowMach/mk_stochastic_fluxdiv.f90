module mk_stochastic_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use analysis_module
  use probin_lowmach_module, only: nscal, rhobar, visc_coef, diff_coef, variance_coef, &
                                   conc_scal, stoch_stress_form, filtering_width

  implicit none

  private

  public :: mk_stochastic_s_fluxdiv, mk_stochastic_m_fluxdiv, &
       fill_stochastic, init_stochastic, destroy_stochastic

  ! Stochastic fluxes for momentum are generated on:
  ! -cell-centered grid for diagonal components
  ! -node-centered (2D) or edge-centered (3D) grid for off-diagonal components
  type(multifab), allocatable, save :: mflux_cc(:), mflux_nd(:), mflux_ed(:,:)

  ! Stochastic fluxes for scalars are face-centered
  type(multifab), allocatable, save :: sflux_fc(:,:)
  
  logical   , save :: warn_bad_c=.false. ! Should we issue warnings about c<0 or c>1
  real(dp_t), save :: kT=1.d0
  real(dp_t), save :: mol_mass(2)=1.d0
  
contains

  ! Note that here we *increment* stoch_s_force so it must be initialized externally!
  subroutine mk_stochastic_s_fluxdiv(mla,the_bc_level,stoch_s_force,s_face,chi,dx,dt)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_s_force(:)
    type(multifab) , intent(in   ) :: s_face(:,:)
    type(multifab) , intent(in   ) :: chi(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: dt

    ! if chi varies in space, average chi to faces
    if (diff_coef < 0) then

    end if

    ! create temporary multifab with scaled random numbers

    ! if chi varies in space, need to multiply pointwise by sqrt(chi)
    if (diff_coef < 0) then

    end if

    ! multiply pointwise by
    ! rho * mu_c^-1 k_b T = rho * c * (1-c)
    

    ! apply boundary conditions

    ! sync up random numbers at boundaries and ghost cells

    ! add divergence to stoch_s_force



  end subroutine mk_stochastic_s_fluxdiv

  ! Note that here we *increment* stoch_m_force so it must be initialized externally!
  subroutine mk_stochastic_m_fluxdiv(mla,the_bc_level,stoch_m_force,eta,dx,dt)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(in   ) :: eta(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: dt


    ! if eta varies in space, average eta to nodes (2D) or edges (3D)
    if (visc_coef < 0) then

    end if

    ! create temporary multifab with scaled random numbers

    ! if eta varies in space, need to multiply pointwise by sqrt(eta)
    if (visc_coef < 0) then

    end if

    ! apply boundary conditions

    ! sync up random numbers at boundaries and ghost cells

    ! add divergence to stoch_m_force



  end subroutine mk_stochastic_m_fluxdiv

  ! fill the stochastic multifabs with random numbers
  subroutine fill_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm

    nlevs = mla%nlevel
    dm = mla%dim

    ! Diagonal components of stochastic stress tensor for momentum
    select case(stoch_stress_form)
       case (0) ! Non-symmetric
          call multifab_fill_random(mflux_cc)
       case default ! Symmetric
          call multifab_fill_random(mflux_cc, variance=2.d0)
       end select

      ! Off-diagonal components of stochastic stress tensor for momentum
      if (dm .eq. 2) then
         ! in 2D, we need 2 random fluxes at each node
         select case(stoch_stress_form)
         case(0) ! Non-symmetric
            call multifab_fill_random(mflux_nd)
         case default ! Symmetric
            call multifab_fill_random(mflux_nd, comp=1)
            do n=1,nlevs
               call multifab_copy_c(mflux_nd(n),2,mflux_nd(n),1)
            end do
         end select
      else if (dm .eq. 3) then
         ! in 3D, we need 2 random fluxes at each edge
         select case(stoch_stress_form)
         case(0) ! Non-symmetric
            call multifab_fill_random(mflux_ed(:,1))
            call multifab_fill_random(mflux_ed(:,2))
            call multifab_fill_random(mflux_ed(:,3))
         case default ! Symmetric
            call multifab_fill_random(mflux_ed(:,1), comp=1)
            call multifab_fill_random(mflux_ed(:,2), comp=1)
            call multifab_fill_random(mflux_ed(:,3), comp=1)
            do n = 1, nlevs
               call multifab_copy_c(mflux_ed(n,1),2,mflux_ed(n,1),1)
               call multifab_copy_c(mflux_ed(n,2),2,mflux_ed(n,2),1)
               call multifab_copy_c(mflux_ed(n,3),2,mflux_ed(n,3),1)
            end do   
         end select             
      end if
            
      if(nscal>1) then ! Stochastic diffusive flux
         do i=1,dm
            call multifab_fill_random(sflux_fc(:,i))
         end do
      end if    

  contains

    ! fill a multifab with random numbers
    subroutine multifab_fill_random(mfab, comp, variance, variance_mfab)
    type(multifab), intent(inout)           :: mfab(:)
    integer       , intent(in   ), optional :: comp ! Only one component
    real(dp_t)    , intent(in   ), optional :: variance
    type(multifab), intent(in   ), optional :: variance_mfab(:)

    integer :: n,box

    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)
  
    !--------------------------------------
    do n=1,size(mfab)
       do box = 1, nfabs(mfab(n))
          if(present(comp)) then
             fp => dataptr(mfab(n),box,comp,1)
          else
             fp => dataptr(mfab(n),box)
          end if
          
          call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers
          
          if(present(variance_mfab)) then ! Must have same distribution
             fpvar => dataptr(variance_mfab(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if
          if(present(variance)) then
             fp=sqrt(variance)*fp
          end if
       end do
    end do

    end subroutine multifab_fill_random

  end subroutine fill_stochastic

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm
    logical :: nodal_temp(mla%dim)
    
    nlevs = mla%nlevel
    dm = mla%dim

    allocate(sflux_fc(mla%nlevel,mla%dim))
    allocate(mflux_cc(mla%nlevel))
    allocate(mflux_nd(mla%nlevel))
    allocate(mflux_ed(mla%nlevel,3))

    do n=1,nlevs
       if(nscal>1) then
          do i=1,dm
             ! we need one face-centered flux for each concentration
            call multifab_build_edge(sflux_fc(n,i),mla%la(n),nscal-1,filtering_width,i)
         end do
       end if
       ! we need dm cell-centered fluxes for momentum
       call multifab_build(mflux_cc(n),mla%la(n),dm,max(1,filtering_width))
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          nodal_temp = .true.
          call multifab_build(mflux_nd(n),mla%la(n),2,filtering_width,nodal_temp)
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(mflux_ed(n,1),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed(n,2),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed(n,3),mla%la(n),2,filtering_width,nodal_temp)
       end if
    end do

  end subroutine init_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm
    
    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       if(nscal>1) then
          do i=1,dm
             call multifab_destroy(sflux_fc(n,i))
         end do
       end if
       call multifab_destroy(mflux_cc(n))
       if (dm .eq. 2) then
          call multifab_destroy(mflux_nd(n))
       else if (dm .eq. 3) then
          call multifab_destroy(mflux_ed(n,1))
          call multifab_destroy(mflux_ed(n,2))
          call multifab_destroy(mflux_ed(n,3))
       end if
    end do
    
    deallocate(mflux_cc,mflux_nd,mflux_ed,sflux_fc)

  end subroutine destroy_stochastic

end module mk_stochastic_fluxdiv_module
