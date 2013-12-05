module external_force_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use multifab_fill_ghost_module
  use probin_common_module
  use probin_multispecies_module

  implicit none

  private

  public :: external_source

contains

  subroutine external_source(mla,rho,fluxdiv,molmtot,Dbar,mass,prob_lo,prob_hi,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(:)
    type(multifab) , intent(in   ) :: molmtot(:)
    real(kind=dp_t), intent(in   ) :: Dbar(:,:)
    real(kind=dp_t), intent(in   ) :: mass(:)
    real(kind=dp_t), intent(in   ) :: prob_lo(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(rho(1)%dim)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local variables
    integer :: i,n,ng,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointer for rho,fluxdiv,molmtot
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)

    dm    = mla%dim         ! dimensionality
    ng    = fluxdiv(1)%ng   ! number of ghost cells
    nlevs = mla%nlevel      ! number of levels

    do n=1,nlevs
       do i = 1, nfabs(rho(n))
          dp  => dataptr(rho(n),i)
          dp1 => dataptr(fluxdiv(n),i)
          dp2 => dataptr(molmtot(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call external_source_2d(dp(:,:,1,:),dp1(:,:,1,:),dp2(:,:,1,1),Dbar(:,:),mass,ng,lo,hi,prob_lo,prob_hi,dx(n,:),time)
          case (3)
             call external_source_3d(dp(:,:,:,:),dp1(:,:,:,:),dp2(:,:,:,1),Dbar(:,:),mass,ng,lo,hi,prob_lo,prob_hi,dx(n,:),time)
          end select
       end do
    enddo

  end subroutine external_source

  subroutine external_source_2d(rho,fluxdiv,molmtot,Dbar,mass,ng,lo,hi,prob_lo,prob_hi,dx,time)

    integer          :: lo(:),hi(:),ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t)  :: fluxdiv(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:)
    real(kind=dp_t)  :: prob_lo(2),prob_hi(2)
    real(kind=dp_t)  :: dx(:),time
    real(kind=dp_t)  :: Dbar(:,:),mass(:) 

    ! local variables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,r,L(2),r_temp

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
   
    !print*, lbound(fluxdiv), ubound(fluxdiv) 
    select case(init_type)
    case(4)
     
     ! for specific box, now start loops over alloted cells    
     do j=lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
        do i=lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
 
             r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2)
               
             ! for m1 = m2 = m and m1 != m2 != m  (for 2-species)
             r_temp = (dexp(-(beta*time) - (r**2*(1.0d0 + time))/(4.0d0*Dbar(1,2)*time))*(16.0d0*beta*Dbar(1,2)**2*&
                      dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI*time + alpha*(r**2 + 4.0d0*beta*Dbar(1,2)*time)))/(64.0d0*&
                      Dbar(1,2)**3*M_PI**2*time**2)
  
             fluxdiv(i,j,1) = fluxdiv(i,j,1) - r_temp 
             fluxdiv(i,j,2) = fluxdiv(i,j,2) + r_temp

        enddo
     enddo
   
    end select 

  end subroutine external_source_2d

  subroutine external_source_3d(rho,fluxdiv,molmtot,Dbar,mass,ng,lo,hi,prob_lo,prob_hi,dx,time)

    integer          :: lo(:),hi(:),ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t)  :: fluxdiv(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t)  :: molmtot(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t)  :: prob_lo(3),prob_hi(3)
    real(kind=dp_t)  :: dx(:),time
    real(kind=dp_t)  :: Dbar(:,:),mass(:) 

    ! local variables
    integer          :: i,j,k
    real(kind=dp_t)  :: x,y,z,r,L(3),r_temp

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    select case(init_type)
    case(4)
    
     ! for specific box, now start loops over alloted cells    
     do k=lo(3), hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
        do j=lo(2), hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
           do i=lo(1), hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0           
   
              r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2)
              
              ! for m1 = m2 = m and m1 != m2 != m (for 2-species)
              r_temp = (dexp(-(beta*time) - (r**2*(1.0d0 + time))/(4.0d0*Dbar(1,2)*time))*&
                       (0.02244839026564582d0*beta*Dbar(1,2)**3.0d0*dexp(r**2/(4.0d0*Dbar(1,2)))*(Dbar(1,2)*&
                       time)**3.5d0 + alpha*Dbar(1,2)**1.5d0*(0.0001259825563796855d0*r**2*(Dbar(1,2)*time)**2.5d0 +& 
                       0.000503930225518742d0*beta*(Dbar(1,2)*time)**3.5d0)))/(Dbar(1,2)**3.0d0*(Dbar(1,2)*time)**5.0d0) 
 
              fluxdiv(i,j,k,1) = fluxdiv(i,j,k,1) - r_temp 
              fluxdiv(i,j,k,2) = fluxdiv(i,j,k,2) + r_temp 
 
           enddo
        enddo
     enddo

    end select 
 
  end subroutine external_source_3d

end module external_force_module 
