module fluid_charge_module

 
  use ml_layout_module
  use convert_stag_module
  use div_and_grad_module
  use define_bc_module
  use bc_module
  use mass_flux_utilities_module
  use matvec_mul_module
  use compute_mixture_properties_module
  use multifab_physbc_module
  use zero_edgeval_module
  use probin_common_module, only: molmass, k_B, rhobar, nspecies, prob_lo, prob_hi
  use probin_charged_module, only: charge_per_mass, dpdt_factor, &
                                   dielectric_const, dielectric_type, &
                                   E_ext_type, E_ext_value, zero_eps_on_wall_type, &
                                   zero_eps_on_wall_left_end, zero_eps_on_wall_right_start

  implicit none

  private

  public :: dot_with_z, dot_with_z_face, compute_charge_coef, &
            enforce_charge_neutrality, implicit_potential_coef, modify_S, &
            compute_permittivity, compute_Lorentz_force, compute_E_ext, &
            zero_eps_on_wall
  
contains

  ! mfdotz = mf dot z
  subroutine dot_with_z(mla,mf,mfdotz,abs_z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mf(:)
    type(multifab) , intent(inout) :: mfdotz(:)
    logical        , intent(in), optional :: abs_z
      ! Dot with abs(z) to estimate norms for relative errors (default=false)

    ! local variables
    integer :: n,nlevs,comp
    real(kind=dp_t) :: z_temp(nspecies)
    
    z_temp=charge_per_mass(1:nspecies)
    if(present(abs_z)) then
      if(abs_z) z_temp=abs(z_temp)
    end if  

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_setval(mfdotz(n),0.d0,all=.true.)
       do comp=1,nspecies
          call multifab_saxpy_3_cc(mfdotz(n),1,z_temp(comp),mf(n),comp,1,all=.true.)
       end do
    end do

  end subroutine dot_with_z

  ! mfdotz = mf dot z
  subroutine dot_with_z_face(mla,mf,mfdotz,abs_z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: mf(:,:)
    type(multifab) , intent(inout) :: mfdotz(:,:)
    logical        , intent(in), optional :: abs_z
      ! Dot with abs(z) to estimate norms for relative errors (default=false)

    ! local variables
    integer :: n,nlevs,i,dm,comp
    real(kind=dp_t) :: z_temp(nspecies)

    z_temp=charge_per_mass(1:nspecies)
    if(present(abs_z)) then
      if(abs_z) z_temp=abs(z_temp)
    end if  

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_setval(mfdotz(n,i),0.d0,all=.true.)
          do comp=1,nspecies
             call multifab_saxpy_3_cc(mfdotz(n,i),1,z_temp(comp),mf(n,i),comp,1,all=.true.)
          end do
       end do
    end do

  end subroutine dot_with_z_face

  ! compute cell-centered mass diffusion coefficients due to charge fluid
  ! charge_coef_i = rho_i z/(n k_B T)
  subroutine compute_charge_coef(mla,rho,Temp,charge_coef)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(in   ) :: rho(:)
    type(multifab ), intent(in   ) :: Temp(:)
    type(multifab ), intent(inout) :: charge_coef(:)

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_3,ng_5
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    real(kind=dp_t), pointer :: dp5(:,:,:,:)

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = rho(1)%ng
    ng_3 = Temp(1)%ng
    ng_5 = charge_coef(1)%ng

    do n=1,nlevs
       do i=1,nfabs(charge_coef(n))
          dp1 => dataptr(rho(n),i)
          dp3 => dataptr(Temp(n),i)
          dp5 => dataptr(charge_coef(n),i)
          lo = lwb(get_box(charge_coef(n),i))
          hi = upb(get_box(charge_coef(n),i))
          select case (dm)
          case (2)
             call compute_charge_coef_2d(dp1(:,:,1,:),ng_1, &
                                         dp3(:,:,1,1),ng_3, &
                                         dp5(:,:,1,:),ng_5, lo,hi)
          case (3)
             call compute_charge_coef_3d(dp1(:,:,:,:),ng_1, &
                                         dp3(:,:,:,1),ng_3, &
                                         dp5(:,:,:,:),ng_5, lo,hi)
          end select
       end do
    end do

  contains

    subroutine compute_charge_coef_2d(rho,ng_1,Temp,ng_3, &
                                      charge_coef,ng_5,lo,hi)
      
      integer         :: lo(:),hi(:),ng_1,ng_3,ng_5
      real(kind=dp_t) ::         rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
      real(kind=dp_t) ::        Temp(lo(1)-ng_3:,lo(2)-ng_3:)
      real(kind=dp_t) :: charge_coef(lo(1)-ng_5:,lo(2)-ng_5:,:)

      ! local variables
      integer :: i,j,comp
      real(kind=dp_t) :: n

      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1

         n = 0.d0
         do comp=1,nspecies
            n = n + rho(i,j,comp)/molmass(comp)
         end do
            
         do comp=1,nspecies
            charge_coef(i,j,comp) = &
                 rho(i,j,comp)*charge_per_mass(comp)/(n*k_B*Temp(i,j))
         end do

      end do
      end do

    end subroutine compute_charge_coef_2d

    subroutine compute_charge_coef_3d(rho,ng_1,Temp,ng_3, &
                                      charge_coef,ng_5,lo,hi)
      
      integer         :: lo(:),hi(:),ng_1,ng_3,ng_5
      real(kind=dp_t) ::         rho(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
      real(kind=dp_t) ::        Temp(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
      real(kind=dp_t) :: charge_coef(lo(1)-ng_5:,lo(2)-ng_5:,lo(3)-ng_5:,:)

      ! local variables
      integer :: i,j,k,comp
      real(kind=dp_t) :: n

      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1

         n = 0.d0
         do comp=1,nspecies
            n = n + rho(i,j,k,comp)/molmass(comp)
         end do
            
         do comp=1,nspecies
            charge_coef(i,j,k,comp) = &
                 rho(i,j,k,comp)*charge_per_mass(comp)/(n*k_B*Temp(i,j,k))
         end do

      end do
      end do
      end do

    end subroutine compute_charge_coef_3d

  end subroutine compute_charge_coef

  subroutine enforce_charge_neutrality(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    type(multifab) :: charge(mla%nlevel)

    integer :: i,dm,nlevs,n,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t) :: charge_temp

    ! pointers into multifab
    real(kind=dp_t), pointer :: rp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_r = rho(1)%ng

    if (nlevs .ne. 1) then
       call bl_error("enforce_charge_neutrality only works with nlevs=1")
    end if

    do n=1,nlevs
       call multifab_build(charge(n),mla%la(n),1,0)
    end do

    ! compute total charge in each cell
    call dot_with_z(mla,rho,charge)

    ! integrate charge over domain
    charge_temp = multifab_sum_c(charge(1),1,1)

    if (parallel_IOProcessor()) then
       print*,'enforce_charge_neutrality: charge before',charge_temp
    end if
    
    ! modify the (0,0,0) cell
    ! for positively charged domains, pick the positive species with the largest rho 
    ! and subtract density.  Pick the negative species with the largest rho_i and 
    ! add density (and vice versa)
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          rp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call enforce_charge_neutrality_2d(rp(:,:,1,:),ng_r,lo,hi,charge_temp)
          case (3)
             call enforce_charge_neutrality_3d(rp(:,:,:,:),ng_r,lo,hi,charge_temp)
          end select
       end do
    end do

    ! compute total charge in each cell
    call dot_with_z(mla,rho,charge)

    ! integrate charge over domain
    charge_temp = multifab_sum_c(charge(1),1,1)

    if (parallel_IOProcessor()) then
       print*,'enforce_charge_neutrality: charge after',charge_temp
    end if

    do n=1,nlevs
       call multifab_destroy(charge(n))
    end do

  contains

    subroutine enforce_charge_neutrality_2d(rho,ng_r,lo,hi,net_charge)
      
      integer         :: lo(:),hi(:),ng_r
      real(kind=dp_t) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
      real(kind=dp_t) :: net_charge

      ! local variables
      integer :: i,j,comp,negative_comp,positive_comp
      logical :: is_positive(nspecies),is_negative(nspecies)
      real(kind=dp_t) :: rho_temp

      is_positive = .false.
      is_negative = .false.
      do comp=1,nspecies
         if (charge_per_mass(comp) .gt. 0.d0) then
            is_positive(comp) = .true.
         end if
         if (charge_per_mass(comp) .lt. 0.d0) then
            is_negative(comp) = .true.
         end if
      end do

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         if (i .eq. 0 .and. j .eq. 0) then

            ! find the negatively charged species with largest rho_i
            rho_temp = -1.d0
            negative_comp = -1
            do comp=1,nspecies
               if (is_negative(comp)) then
                  if (rho(i,j,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,comp)
                     negative_comp = comp
                  end if
               end if
            end do

            ! find the positively charged species with largest rho_i
            rho_temp = -1.d0
            positive_comp = -1
            do comp=1,nspecies
               if (is_positive(comp)) then
                  if (rho(i,j,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,comp)
                     positive_comp = comp
                  end if
               end if
            end do

            if (net_charge .lt. 0.d0) then
               rho(i,j,negative_comp) = rho(i,j,negative_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(negative_comp))
               rho(i,j,positive_comp) = rho(i,j,positive_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(positive_comp))
            else
               rho(i,j,positive_comp) = rho(i,j,positive_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(positive_comp))
               rho(i,j,negative_comp) = rho(i,j,negative_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(negative_comp))
            end if

         end if

      end do
      end do

    end subroutine enforce_charge_neutrality_2d

    subroutine enforce_charge_neutrality_3d(rho,ng_r,lo,hi,net_charge)
      
      integer         :: lo(:),hi(:),ng_r
      real(kind=dp_t) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
      real(kind=dp_t) :: net_charge

      ! local variables
      integer :: i,j,k,comp,negative_comp,positive_comp
      logical :: is_positive(nspecies),is_negative(nspecies)
      real(kind=dp_t) :: rho_temp

      is_positive = .false.
      is_negative = .false.
      do comp=1,nspecies
         if (charge_per_mass(comp) .gt. 0.d0) then
            is_positive(comp) = .true.
         end if
         if (charge_per_mass(comp) .lt. 0.d0) then
            is_negative(comp) = .true.
         end if
      end do

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         if (i .eq. 0 .and. j .eq. 0 .and. k .eq. 0) then

            ! find the negatively charged species with largest rho_i
            rho_temp = -1.d0
            negative_comp = -1
            do comp=1,nspecies
               if (is_negative(comp)) then
                  if (rho(i,j,k,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,k,comp)
                     negative_comp = comp
                  end if
               end if
            end do

            ! find the positively charged species with largest rho_i
            rho_temp = -1.d0
            positive_comp = -1
            do comp=1,nspecies
               if (is_positive(comp)) then
                  if (rho(i,j,k,comp) .gt. rho_temp) then
                     rho_temp = rho(i,j,k,comp)
                     positive_comp = comp
                  end if
               end if
            end do

            if (net_charge .lt. 0.d0) then
               rho(i,j,k,negative_comp) = rho(i,j,k,negative_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(negative_comp))
               rho(i,j,k,positive_comp) = rho(i,j,k,positive_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(positive_comp))
            else
               rho(i,j,k,positive_comp) = rho(i,j,k,positive_comp) &
                    - abs(0.5d0*net_charge/charge_per_mass(positive_comp))
               rho(i,j,k,negative_comp) = rho(i,j,k,negative_comp) &
                    + abs(0.5d0*net_charge/charge_per_mass(negative_comp))
            end if

         end if

      end do
      end do
      end do

    end subroutine enforce_charge_neutrality_3d

  end subroutine enforce_charge_neutrality


  ! compute face-centered A_\Phi (an nspecies vector)
  ! compute vector rho W z / (n k_B T) on cell centers and average to faces
  ! compute tensor rho W chi on cell centers and average to faces
  ! multiply them together and store the resulting vector in A_Phi
  subroutine implicit_potential_coef(mla,rho,Temp,A_Phi,the_bc_tower,rhoWchi_face_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(in   ) :: rho(:)
    type(multifab ), intent(in   ) :: Temp(:)
    type(multifab ), intent(inout) :: A_Phi(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab ), intent(in), optional :: rhoWchi_face_in(:,:)
    
    ! local
    type(multifab) :: rhotot_temp (mla%nlevel)
    type(multifab) :: charge_coef (mla%nlevel)
    type(multifab) :: molarconc   (mla%nlevel)
    type(multifab) :: molmtot     (mla%nlevel)
    type(multifab) :: D_therm     (mla%nlevel)
    type(multifab) :: Hessian     (mla%nlevel)
    type(multifab) :: D_bar       (mla%nlevel)
    type(multifab) :: rhoWchi     (mla%nlevel)
    type(multifab) :: rhoWchi_face(mla%nlevel,mla%dim)

    integer :: n,nlevs,i,dm

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(rhotot_temp(n),  mla%la(n), 1,           rho(n)%ng)
       call multifab_build(charge_coef(n),  mla%la(n), nspecies,    1)
       call multifab_build(molarconc(n),    mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(molmtot(n),      mla%la(n), 1,           rho(n)%ng)
       call multifab_build(D_therm(n),      mla%la(n), nspecies,    rho(n)%ng)
       call multifab_build(Hessian(n),      mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(D_bar(n),        mla%la(n), nspecies**2, rho(n)%ng)
       call multifab_build(rhoWchi(n),      mla%la(n), nspecies**2, rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(rhoWchi_face(n,i),  mla%la(n), nspecies**2, 0, i)
       end do
    end do

    call compute_rhotot(mla,rho,rhotot_temp,ghost_cells_in=.true.)

    ! compute rho W z / (n k_B T) on cell centers
    call compute_charge_coef(mla,rho,Temp,charge_coef)

    ! average charge_coef to faces (store in A_Phi)
    call average_cc_to_face(nlevs,charge_coef,A_Phi,1,c_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array,.true.)

    ! rhoWchi on faces
    if (present(rhoWchi_face_in)) then
       ! copy rhoWchi to faces
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(rhoWchi_face(n,i),1,rhoWchi_face_in(n,i),1,nspecies**2,0)
          end do
       end do
    else
       ! compute rhoWchi on cell centers
       call compute_molconc_molmtot(mla,rho,rhotot_temp,molarconc,molmtot)
       call compute_mixture_properties(mla,rho,rhotot_temp,D_bar,D_therm,Hessian)
       call compute_rhoWchi(mla,rho,rhotot_temp,molarconc,rhoWchi,D_bar)

       ! average rhoWchi to faces
       call average_cc_to_face(nlevs, rhoWchi, rhoWchi_face, 1, tran_bc_comp, &
                               nspecies**2, the_bc_tower%bc_tower_array, .false.) 
    end if

    ! multiply A_Phi by rhoWchi_face
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, A_Phi(n,i), rhoWchi_face(n,i), nspecies)
       end do
    end do    

    ! deallocate memory
    do n=1,nlevs
       call multifab_destroy(rhotot_temp(n))
       call multifab_destroy(charge_coef(n))
       call multifab_destroy(molarconc(n))
       call multifab_destroy(molmtot(n))
       call multifab_destroy(rhoWchi(n))
       call multifab_destroy(D_bar(n))
       call multifab_destroy(D_therm(n))
       call multifab_destroy(Hessian(n))
       do i=1,dm
          call multifab_destroy(rhoWchi_face(n,i))
       end do
    end do

  end subroutine implicit_potential_coef

  subroutine modify_S(mla,rho,S_inc,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(in   ) :: rho(:)
    type(multifab ), intent(inout) :: S_inc(:)
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    
    real(kind=dp_t) :: S_avg

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = rho(1)%ng
    ng_2 = S_inc(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(S_inc(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call modify_S_2d(dp1(:,:,1,:),ng_1,dp2(:,:,1,1),ng_2,lo,hi,dt)
          case (3)
             call modify_S_3d(dp1(:,:,:,:),ng_1,dp2(:,:,:,1),ng_2,lo,hi,dt)
          end select
       end do
    end do

    ! subtract off average so S sums to zero
    S_avg = multifab_sum_c(S_inc(1),1,1) / multifab_volume(S_inc(1))
    call multifab_sub_sub_s_c(S_inc(1),1,S_avg,1,0)

  contains

    subroutine modify_S_2d(rho,ng_1,S_inc,ng_2,lo,hi,dt)
      
      integer         :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t) :: dt
      real(kind=dp_t) ::   rho(lo(1)-ng_1:,lo(2)-ng_1:,:)
      real(kind=dp_t) :: S_inc(lo(1)-ng_2:,lo(2)-ng_2:)

      ! local variables
      integer :: i,j,comp
      real(kind=dp_t) :: rhorhobar

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         rhorhobar = 0.d0
         do comp=1,nspecies
            rhorhobar = rhorhobar + rho(i,j,comp)/rhobar(comp)
         end do

         S_inc(i,j) = S_inc(i,j) - (dpdt_factor/dt)*(rhorhobar - 1.d0)/rhorhobar

      end do
      end do

    end subroutine modify_S_2d

    subroutine modify_S_3d(rho,ng_1,S_inc,ng_2,lo,hi,dt)
      
      integer         :: lo(:),hi(:),ng_1,ng_2
      real(kind=dp_t) :: dt
      real(kind=dp_t) ::   rho(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:,:)
      real(kind=dp_t) :: S_inc(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)

      ! local variables
      integer :: i,j,k,comp
      real(kind=dp_t) :: rhorhobar

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

         rhorhobar = 0.d0
         do comp=1,nspecies
            rhorhobar = rhorhobar + rho(i,j,k,comp)/rhobar(comp)
         end do

         S_inc(i,j,k) = S_inc(i,j,k) - (dpdt_factor/dt)*(rhorhobar - 1.d0)/rhorhobar

      end do
      end do
      end do

    end subroutine modify_S_3d

  end subroutine modify_S

  subroutine compute_permittivity(mla,permittivity,rho,rhotot,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: permittivity(:)
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2,ng_3
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    if (dielectric_type .eq. 0) then
       return
    end if
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = permittivity(1)%ng
    ng_2 = rho(1)%ng
    ng_3 = rhotot(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp1 => dataptr(permittivity(n),i)
          dp2 => dataptr(rho(n),i)
          dp3 => dataptr(rhotot(n),i)
          lo = lwb(get_box(permittivity(n),i))
          hi = upb(get_box(permittivity(n),i))
          select case (dm)
          case (2)
             call compute_permittivity_2d(dp1(:,:,1,1),ng_1,dp2(:,:,1,:),ng_2, &
                                          dp3(:,:,1,1),ng_3,lo,hi)
          case (3)
             call compute_permittivity_3d(dp1(:,:,:,1),ng_1,dp2(:,:,:,:),ng_2, &
                                          dp3(:,:,:,1),ng_3,lo,hi)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_fill_boundary(permittivity(n))
       call multifab_physbc(permittivity(n),1,c_bc_comp,1,the_bc_tower%bc_tower_array(n))
    end do

  contains

    subroutine compute_permittivity_2d(permittivity,ng_1,rho,ng_2,rhotot,ng_3,lo,hi)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: permittivity(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) ::          rho(lo(1)-ng_2:,lo(2)-ng_2:,:)
      real(kind=dp_t) ::       rhotot(lo(1)-ng_3:,lo(2)-ng_3:)
      
      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            call compute_permittivity_local(permittivity(i,j),rho(i,j,:),rhotot(i,j))

         end do
      end do

    end subroutine compute_permittivity_2d

    subroutine compute_permittivity_3d(permittivity,ng_1,rho,ng_2,rhotot,ng_3,lo,hi)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: permittivity(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) ::          rho(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:)
      real(kind=dp_t) ::       rhotot(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
      
      ! local
      integer :: i,j,k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               call compute_permittivity_local(permittivity(i,j,k),rho(i,j,k,:),rhotot(i,j,k))


            end do
         end do
      end do
      
    end subroutine compute_permittivity_3d

    subroutine compute_permittivity_local(perm,rho,rhotot)

      real(kind=dp_t) :: perm, rho(:), rhotot

      ! local
      real(kind=dp_t) :: c1

      if (dielectric_type .eq. 1) then
         c1 = rho(1)/rhotot
         perm = (1.d0+c1)*dielectric_const
      end if

    end subroutine compute_permittivity_local

  end subroutine compute_permittivity

  subroutine compute_Lorentz_force(mla,Lorentz_force,grad_Epot,permittivity,charge, &
                                   dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: Lorentz_force(:,:)
    type(multifab) , intent(in   ) :: grad_Epot(:,:)
    type(multifab) , intent(in   ) :: permittivity(:)
    type(multifab) , intent(in   ) :: charge(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: i,n,dm,nlevs
    integer :: ng_1,ng_2,ng_3
    integer :: lo(mla%dim),hi(mla%dim)

    type(multifab) :: temp_cc(mla%nlevel)

    ! pointers into multifabs
    real(kind=dp_t), pointer :: dp1x(:,:,:,:)
    real(kind=dp_t), pointer :: dp1y(:,:,:,:)
    real(kind=dp_t), pointer :: dp1z(:,:,:,:)
    real(kind=dp_t), pointer :: dp2x(:,:,:,:)
    real(kind=dp_t), pointer :: dp2y(:,:,:,:)
    real(kind=dp_t), pointer :: dp2z(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    logical :: use_qE_Lorentz
    
    dm = mla%dim
    nlevs = mla%nlevel

    use_qE_Lorentz = .false.

    if (use_qE_Lorentz) then

       call average_cc_to_face(nlevs,charge,Lorentz_force,1,c_bc_comp,1, &
                               the_bc_tower%bc_tower_array)

       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_s(Lorentz_force(n,i),-1.d0,0)
             call multifab_mult_mult_c(Lorentz_force(n,i),1,grad_Epot(n,i),1,1,0)
          end do
       end do

       return

    end if

    ! temporary multifab to hold div(-eps*E)
    do n=1,nlevs
       call multifab_build(temp_cc(n),mla%la(n),1,1)
    end do

    ! Lorentz force = E div (eps*E) - (1/2) E^2 grad(eps)

    ! discretize E div(eps*E)
    ! start by averaging epsilon to faces
    ! store this in Lorentz_force temporarily
    if (dielectric_type .eq. 0) then
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(Lorentz_force(n,i),dielectric_const,all=.true.)
          end do
       end do
    else
       call average_cc_to_face(nlevs,permittivity,Lorentz_force,1,c_bc_comp,1, &
                               the_bc_tower%bc_tower_array) 

      if (zero_eps_on_wall_type .gt. 0) then
          ! set beta to set to zero on certain boundary faces
          call zero_eps_on_wall(mla,Lorentz_force,dx)
       end if

    end if

    ! multiply by grad_Epot = -E to get eps*E on faces
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(Lorentz_force(n,i),1,grad_Epot(n,i),1,1,0)
       end do
    end do

    ! take divergence of -eps*E and put it in cc_temp
    call compute_div(mla,Lorentz_force,temp_cc,dx,1,1,1)

    ! fill ghost cells for div(-eps*E)
    do n=1,nlevs
       call multifab_fill_boundary(temp_cc(n))
       call multifab_physbc(temp_cc(n),1,c_bc_comp,1,the_bc_tower%bc_tower_array(n))
    end do

    ! average div(-eps*E) to faces, store in Lorentz_force
    call average_cc_to_face(nlevs,temp_cc,Lorentz_force,1,c_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! multiply by -E to get E*div(eps*E) on faces
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(Lorentz_force(n,i),1,grad_Epot(n,i),1,1,0)
       end do
    end do       

    if (dielectric_type .ne. 0) then

       ! spatially-varying permittivity
       ! add -(1/2) E^2 grad(eps)

       ng_1 = Lorentz_force(1,1)%ng
       ng_2 = grad_Epot(1,1)%ng
       ng_3 = permittivity(1)%ng

       do n=1,nlevs
          do i=1,nfabs(Lorentz_force(n,1))
             dp1x => dataptr(Lorentz_force(n,1),i)
             dp1y => dataptr(Lorentz_force(n,2),i)
             dp2x => dataptr(grad_Epot(n,1),i)
             dp2y => dataptr(grad_Epot(n,2),i)
             dp3  => dataptr(permittivity(n),i)
             lo = lwb(get_box(Lorentz_force(n,1),i))
             hi = upb(get_box(Lorentz_force(n,1),i))
             select case (dm)
             case (2)
                call Lorentz_force_gradeps_2d(dp1x(:,:,1,1),dp1y(:,:,1,1),ng_1, &
                                              dp2x(:,:,1,1),dp2y(:,:,1,1),ng_2, &
                                              dp3(:,:,1,1),ng_3,lo,hi,dx(n,:))
             case (3)
                dp1z => dataptr(Lorentz_force(n,3),i)
                dp2z => dataptr(grad_Epot(n,3),i)
                call Lorentz_force_gradeps_3d(dp1x(:,:,:,1),dp1y(:,:,:,1),dp1z(:,:,:,1),ng_1, &
                                              dp2x(:,:,:,1),dp2y(:,:,:,1),dp2z(:,:,:,1),ng_2, &
                                              dp3(:,:,:,1),ng_3,lo,hi,dx(n,:))
             end select
          end do
       end do

    end if

    ! set force on walls to be zero since normal velocity is zero
    do n=1,nlevs
       call zero_edgeval_walls(Lorentz_force(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    do n=1,nlevs
       call multifab_destroy(temp_cc(n))
    end do

  contains

    subroutine Lorentz_force_gradeps_2d(forcex,forcey,ng_1,Ex,Ey,ng_2,perm,ng_3,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:)
      real(kind=dp_t) ::     Ex(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::     Ey(lo(1)-ng_2:,lo(2)-ng_2:)
      real(kind=dp_t) ::   perm(lo(1)-ng_3:,lo(2)-ng_3:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j

      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         forcex(i,j) = forcex(i,j) - 0.5d0 * ( Ex(i,j)**2 + (0.25d0*(Ey(i,j)+Ey(i,j+1)+Ey(i-1,j)+Ey(i-1,j+1)))*2 ) &
              * (perm(i,j)-perm(i-1,j)) / dx(1)
      end do
      end do

      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
         forcey(i,j) = forcey(i,j) - 0.5d0 * ( Ey(i,j)**2 + (0.25d0*(Ex(i,j)+Ex(i+1,j)+Ex(i,j-1)+Ex(i+1,j-1)))**2 ) &
              * (perm(i,j)-perm(i,j-1)) / dx(1)
      end do
      end do

    end subroutine Lorentz_force_gradeps_2d

    subroutine Lorentz_force_gradeps_3d(forcex,forcey,forcez,ng_1,Ex,Ey,Ez,ng_2,perm,ng_3,lo,hi,dx)

      integer         :: lo(:),hi(:),ng_1,ng_2,ng_3
      real(kind=dp_t) :: forcex(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcey(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) :: forcez(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)
      real(kind=dp_t) ::     Ex(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) ::     Ey(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) ::     Ez(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:)
      real(kind=dp_t) ::   perm(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)
      real(kind=dp_t) :: dx(:)

      ! local variables
      integer :: i,j,k

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
         forcex(i,j,k) = forcex(i,j,k) - 0.5d0 * ( Ex(i,j,k)**2 &
              + (0.25d0*(Ey(i,j,k)+Ey(i,j+1,k)+Ey(i-1,j,k)+Ey(i-1,j+1,k)))*2 &
              + (0.25d0*(Ez(i,j,k)+Ez(i,j,k+1)+Ez(i-1,j,k)+Ez(i-1,j,k+1)))*2 ) &
              * (perm(i,j,k)-perm(i-1,j,k)) / dx(1)
      end do
      end do
      end do

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
         forcey(i,j,k) = forcey(i,j,k) - 0.5d0 * ( Ey(i,j,k)**2 &
              + (0.25d0*(Ex(i,j,k)+Ex(i+1,j,k)+Ex(i,j-1,k)+Ex(i+1,j-1,k)))**2 &
              + (0.25d0*(Ez(i,j,k)+Ez(i,j,k+1)+Ez(i,j-1,k)+Ez(i,j-1,k+1)))**2 ) &
              * (perm(i,j,k)-perm(i,j-1,k)) / dx(2)
      end do
      end do
      end do

      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
         forcez(i,j,k) = forcez(i,j,k) - 0.5d0 * ( Ez(i,j,k)**2 &
              + (0.25d0*(Ex(i,j,k)+Ex(i+1,j,k)+Ex(i,j,k-1)+Ex(i+1,j,k-1)))**2 &
              + (0.25d0*(Ey(i,j,k)+Ey(i,j+1,k)+Ey(i,j,k-1)+Ey(i,j+1,k-1)))**2 ) &
              * (perm(i,j,k)-perm(i,j,k-1)) / dx(3)
      end do
      end do
      end do

    end subroutine Lorentz_force_gradeps_3d

  end subroutine compute_Lorentz_force

  subroutine compute_E_ext(mla,E_ext)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: E_ext(:,:)

    integer :: n, nlevs, dm, i

    nlevs = mla%nlevel
    dm = mla%dim

    if (E_ext_type .eq. 1) then

       do n=1,nlevs
          do i=1,dm
             call multifab_setval(E_ext(n,i),E_ext_value(i),all=.true.)
          end do
       end do

    end if

  end subroutine compute_E_ext

  subroutine zero_eps_on_wall(mla,beta,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: beta(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: i,dm,nlevs,n,ng_b
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: bpx(:,:,:,:)
    real(kind=dp_t), pointer :: bpy(:,:,:,:)
    real(kind=dp_t), pointer :: bpz(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_b = beta(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(beta(n,1))
          lo = lwb(get_box(beta(n,1),i))
          hi = upb(get_box(beta(n,1),i))
          bpy => dataptr(beta(n,2),i)
          select case (dm)
          case (2)
             call zero_eps_on_wall_2d(bpy(:,:,1,1),ng_b,lo,hi,dx(n,:))
          case (3)
             call bl_error("zero_eps_on_wall_3d not written yet")
          end select
       end do
    end do

  end subroutine zero_eps_on_wall

  subroutine zero_eps_on_wall_2d(betay,ng_b,lo,hi,dx)
      
    integer         :: lo(:),hi(:),ng_b
    real(kind=dp_t) :: betay(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t) :: dx(:)

    integer :: i,j

    real(kind=dp_t) :: x,y,Lx

    Lx = prob_hi(1) - prob_lo(1)    

    if (zero_eps_on_wall_type .eq. 1) then

       ! y-faces
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dx(2)*j
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1)*(i+0.5d0)

             if (y .eq. 0.d0 .and. (x .le. zero_eps_on_wall_left_end*Lx .or. x .ge. zero_eps_on_wall_right_start*Lx) ) then
                betay(i,j) = 0.d0
             end if

          end do
       end do

    end if

  end subroutine zero_eps_on_wall_2d

end module fluid_charge_module
