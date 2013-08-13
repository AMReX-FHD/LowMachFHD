module init_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use convert_stag_module
  use bc_module
  use probin_lowmach_module, only: rhobar, diff_coef, visc_coef, &
                                   smoothing_width, c_init, material_properties
  use probin_common_module , only: prob_lo, prob_hi, prob_type

  implicit none

  private

  public :: init, compute_eta, compute_chi, compute_kappa

contains

  ! Important note: For periodic boundaries, the init routines should fill out
  ! *both* sides of the domain with values, even though this is duplicate information
  ! We ensure the two sides are bitwise identical, but low and/or high sides may win

  subroutine init(m,s,p,dx,mla,time)

    type(multifab) , intent(inout) :: m(:,:),s(:),p(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t), pointer :: mxp(:,:,:,:), myp(:,:,:,:), mzp(:,:,:,:)
    real(kind=dp_t), pointer :: sop(:,:,:,:), pp(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: nlevs,n,i,ng_m,ng_s,ng_p,dm

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_m = m(1,1)%ng
    ng_s = s(1)%ng
    ng_p = p(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(s(n))
          mxp => dataptr(m(n,1),i)
          myp => dataptr(m(n,2),i)
          sop => dataptr(s(n),i)
          pp  => dataptr(p(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call init_2d(mxp(:,:,1,1), myp(:,:,1,1), &
                          sop(:,:,1,:), pp(:,:,1,1), &
                          lo, hi, ng_m, ng_s, ng_p, dx(n,:), time)
          case (3)
             mzp => dataptr(m(n,3),i)
             call init_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), &
                          sop(:,:,:,:), pp(:,:,:,1), &
                          lo, hi, ng_m, ng_s, ng_p, dx(n,:), time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(m(n,i))
       end do

    enddo

  end subroutine init

  subroutine init_2d(mx,my,s,p,lo,hi,ng_m,ng_s,ng_p,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s, ng_p
    real(kind=dp_t), intent(inout) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  p(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y,y1,y2,r
    real(kind=dp_t) :: one_third_domain1,one_third_domain2

    select case (prob_type)
    case (0)

       ! constant density
       mx = 0.d0
       my = 0.d0

       p = 0.d0

       ! set c to c_init(1)
       s(:,:,2) = c_init(1)

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! compute rho with eos
             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))

             ! compute rho*c
             s(i,j,2) = s(i,j,1)*s(i,j,2)

          end do
       end do

    case (1)

       ! lo density spherical bubble

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! set c using tanh smoothing
             s(i,j,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                  (1.d0 + tanh((r-2.5d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))

             ! compute rho using eos
             s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))

             ! compute rho*c
             s(i,j,2) = s(i,j,1)*s(i,j,2)
          enddo
       enddo

    case (2)

       ! bilayer interface (stripe)

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
       one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)

       do j=lo(2),hi(2)
          y1 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain1)
          y2 =(prob_lo(2) + dx(2)*(dble(j)+0.5d0) - one_third_domain2)
        
          do i=lo(1),hi(1)
             ! tanh smoothing
             if(abs(smoothing_width)>epsilon(1.0d0)) then
                s(i,j,2) = c_init(1)+ 0.5d0*(c_init(2)-c_init(1))*&
                   (tanh(y1/(smoothing_width*dx(2))) - tanh(y2/(smoothing_width*dx(2))))
                s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))
                s(i,j,2) = s(i,j,1)*s(i,j,2)
             else
                ! Try to initialize exactly as we do in the HDMD simulations,
                ! with finite-volume averaging of sharp interface
                if((y1<-0.5d0*dx(2)).or.(y2>0.5d0*dx(2))) then
                   s(i,j,2) = 0
                   s(i,j,1) = rhobar(2)
                else if((y1>0.5d0*dx(2)).and.(y2<-0.5d0*dx(2))) then
                   s(i,j,2) = rhobar(1)
                   s(i,j,1) = rhobar(1)
                else if(y1 <= 0.5d0*dx(2)) then
                   s(i,j,2) = (max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(1)
                   s(i,j,1) = s(i,j,2) + (1.0d0-max(0.0d0,min(0.5d0+y1/dx(2),1.0d0)))*rhobar(2)
                else 
                   s(i,j,2) = (1.0d0-max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(1)
                   s(i,j,1) = s(i,j,2) + (max(0.0d0,min(0.5d0+y2/dx(2),1.0d0)))*rhobar(2)
                end if   
             end if  
          enddo
       enddo

    case (3)

       ! one fluid on top of another
       ! could be used for diffusive mixing or Rayleigh-Taylor,
       ! depending on where the heavier fluid is and if gravity is on

       mx = 0.d0
       my = 0.d0

       p = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          do i=lo(1),hi(1)

             if (y .lt. y1) then
                s(i,j,2) = c_init(1)
                s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))             
             else
                s(i,j,2) = c_init(2)
                s(i,j,1) = 1.0d0/(s(i,j,2)/rhobar(1)+(1.0d0-s(i,j,2))/rhobar(2))   
             end if
             s(i,j,2) = s(i,j,1)*s(i,j,2)
             
          end do
       end do

    case default

       call bl_error("init_2d: invalid prob_type")

    end select

  end subroutine init_2d

  subroutine init_3d(mx,my,mz,s,p,lo,hi,ng_m,ng_s,ng_p,dx,time)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m, ng_s, ng_p
    real(kind=dp_t), intent(inout) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(inout) ::  p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z,r

    select case (prob_type)
    case (0)

       ! constant density
       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       ! set c to c_init(1)
       s(:,:,:,2) = c_init(1)

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                ! compute rho with eos
                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))

                ! compute rho*c
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)

             end do
          end do
       end do

    case (1)

       ! lo density spherical bubble

       mx = 0.d0
       my = 0.d0
       mz = 0.d0

       p = 0.d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))
                
                r = sqrt (x**2 + y**2 + z**2)

                ! set c using tanh smoothing
                s(i,j,k,2) = c_init(1) + 0.5d0*(c_init(2)-c_init(1))* &
                     (1.d0 + tanh((r-2.0d0*smoothing_width*dx(1))/(smoothing_width*dx(1))))
                
                ! compute rho using eos
                s(i,j,k,1) = 1.0d0/(s(i,j,k,2)/rhobar(1)+(1.0d0-s(i,j,k,2))/rhobar(2))

                ! compute rho*c
                s(i,j,k,2) = s(i,j,k,1)*s(i,j,k,2)
             enddo
          enddo
       end do

    case default

       call bl_error("init_3d: invalid prob_type")
          
    end select

  end subroutine init_3d

  subroutine compute_chi(mla,chi,chi_fc,prim,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: chi(:)
    type(multifab) , intent(inout) :: chi_fc(:,:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    integer :: nlevs,dm,i,n,ng_c,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_c = chi(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_c) then
       call bl_error('ng_p must be no less than ng_c')
    end if

    do n=1,nlevs
       do i=1,nfabs(chi(n))
          cp => dataptr(chi(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(chi(n), i))
          hi = upb(get_box(chi(n), i))
          select case (dm)
          case (2)
             call compute_chi_2d(cp(:,:,1,1),ng_c,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_chi_3d(cp(:,:,:,1),ng_c,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

    call average_cc_to_face(nlevs,chi,chi_fc,1,tran_bc_comp,1,the_bc_level)

  end subroutine compute_chi

  subroutine compute_chi_2d(chi,ng_c,prim,ng_p,lo,hi,dx)

    ! compute chi in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j

    select case (prob_type)
    case (2)

       ! simple test function
       ! chi = chi0*(1 + a*c + b*c^2)
       do j=lo(2)-ng_c,hi(2)+ng_c
          do i=lo(1)-ng_c,hi(1)+ng_c
             conc = max(min(prim(i,j,2), 1.d0), 0.d0)
             chi(i,j) = abs(diff_coef)*(1.d0 + material_properties(1,1)*conc + &
                                            material_properties(2,1)*conc**2)
          end do
       end do

    case default

       chi = diff_coef

    end select

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(chi,ng_c,prim,ng_p,lo,hi,dx)

    ! compute chi in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_p
    real(kind=dp_t), intent(inout) ::  chi(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j,k

    select case (prob_type)
    case (2)

       ! simple test function
       do k=lo(3)-ng_c,hi(3)+ng_c
          do j=lo(2)-ng_c,hi(2)+ng_c
             do i=lo(1)-ng_c,hi(1)+ng_c
                conc = max(min(prim(i,j,k,2), 1.0_dp_t), 0.d0)
                chi(i,j,k) = abs(diff_coef)*(1.d0 + material_properties(1,1)*conc + &
                                                    material_properties(2,1)*conc**2)
             end do
          end do
       end do

       ! simple test function
       ! chi = chi0*(1 + a*c + b*c^2)

    case default

       chi = diff_coef

    end select

  end subroutine compute_chi_3d

  subroutine compute_eta(mla,eta,eta_ed,prim,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    integer :: nlevs,dm,i,n,ng_e,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_e = eta(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_e) then
       call bl_error('ng_p must be no less than ng_e')
    end if

    do n=1,nlevs
       do i=1,nfabs(eta(n))
          cp => dataptr(eta(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(eta(n), i))
          hi = upb(get_box(eta(n), i))
          select case (dm)
          case (2)
             call compute_eta_2d(cp(:,:,1,1),ng_e,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_eta_3d(cp(:,:,:,1),ng_e,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

    if (dm .eq. 2) then
       call average_cc_to_node(nlevs,eta,eta_ed(:,1),1,tran_bc_comp,1,the_bc_level)
    else if (dm .eq. 3) then
       call average_cc_to_edge(nlevs,eta,eta_ed,1,tran_bc_comp,1,the_bc_level)
    end if

  end subroutine compute_eta

  subroutine compute_eta_2d(eta,ng_e,prim,ng_p,lo,hi,dx)

    ! compute eta in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j

    select case (prob_type)
    case (2)

       ! simple test function
       ! eta = eta0*(1+a*c) / (1+b*c)
       do j=lo(2)-ng_e,hi(2)+ng_e
       do i=lo(1)-ng_e,hi(1)+ng_e
          conc = max(min(prim(i,j,2), 1.d0), 0.d0)
          eta(i,j) = abs(visc_coef)*(1.d0 + material_properties(1,2)*conc) / &
                                    (1.d0 + material_properties(2,2)*conc)
       end do
       end do

    case default

       eta = visc_coef

    end select

  end subroutine compute_eta_2d

  subroutine compute_eta_3d(eta,ng_e,prim,ng_p,lo,hi,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng_e, ng_p
    real(kind=dp_t), intent(inout) ::  eta(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: conc
    integer :: i,j,k

    select case (prob_type)
    case (2)

       ! simple test function
       do k=lo(3)-ng_e,hi(3)+ng_e
       do j=lo(2)-ng_e,hi(2)+ng_e
       do i=lo(1)-ng_e,hi(1)+ng_e
          conc = max(min(prim(i,j,k,2), 1.d0), 0.d0)
          eta(i,j,k) = abs(visc_coef)*(1.d0 + material_properties(1,2)*conc) / &
                                      (1.d0 + material_properties(2,2)*conc)

       end do
       end do
       end do

    case default

       eta = visc_coef

    end select

  end subroutine compute_eta_3d

  subroutine compute_kappa(mla,kappa,prim,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(in   ) :: prim(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    integer :: nlevs,dm,i,n,ng_k,ng_p
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    ng_k = kappa(1)%ng
    ng_p = prim(1)%ng

    nlevs = mla%nlevel
    dm = mla%dim

    if(ng_p<ng_k) then
       call bl_error('ng_p must be no less than ng_k')
    end if

    do n=1,nlevs
       do i=1,nfabs(kappa(n))
          cp => dataptr(kappa(n), i)
          pp => dataptr(prim(n), i)
          lo = lwb(get_box(kappa(n), i))
          hi = upb(get_box(kappa(n), i))
          select case (dm)
          case (2)
             call compute_kappa_2d(cp(:,:,1,1),ng_k,pp(:,:,1,:),ng_p,lo,hi,dx(n,:))
          case (3)
             call compute_kappa_3d(cp(:,:,:,1),ng_k,pp(:,:,:,:),ng_p,lo,hi,dx(n,:))
          end select
       end do
    end do

  end subroutine compute_kappa

  subroutine compute_kappa_2d(kappa,ng_k,prim,ng_p,lo,hi,dx)

    ! compute kappa in valid AND ghost regions
    ! the ghost cells in prim have already been filled properly

    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:)
    real(kind=dp_t), intent(in   ) ::  prim(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    select case (prob_type)
    case default
       kappa = 1.d0
    end select

  end subroutine compute_kappa_2d

  subroutine compute_kappa_3d(kappa,ng_k,prim,ng_p,lo,hi,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng_k, ng_p
    real(kind=dp_t), intent(inout) :: kappa(lo(1)-ng_k:,lo(2)-ng_k:,lo(3)-ng_k:)
    real(kind=dp_t), intent(in   ) ::  prim(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    select case (prob_type)
    case default
       kappa = 1.d0
    end select

  end subroutine compute_kappa_3d

end module init_module
