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

    dm    = mla%dim     ! dimensionality
    ng    = fluxdiv(1)%ng   ! number of ghost cells
    nlevs = mla%nlevel  ! number of levels

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
    real(kind=dp_t)  :: x,y,r,tau,L(2)

    tau = 1.0d0
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

             ! for m1 \eq m2 \eq m
             
             fluxdiv(i,j,1) = fluxdiv(i,j,1) - dexp(-r**2/(4.0d0*Dbar(1,2)*time) - time)/(4.0d0*Dbar(1,2)*M_PI*time) 
             fluxdiv(i,j,2) = fluxdiv(i,j,2) + dexp(-r**2/(4.0d0*Dbar(1,2)*time) - time)/(4.0d0*Dbar(1,2)*M_PI*time) 
             !fluxdiv(i,j,1) = fluxdiv(i,j,1) - (dexp(-(r**2*(1.0d0 + time))/(4.0d0*Dbar(1,2)*time) - time/tau)*&
             !                 (4.0d0*Dbar(1,2)*(1.0d0 + 4.0d0*Dbar(1,2)*dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI)*time +&
             !                 r**2*tau))/(64.0d0*Dbar(1,2)**3*M_PI**2*time**2*tau) 
           ! 
           !  fluxdiv(i,j,2) = fluxdiv(i,j,2) + (dexp(-(r**2*(1.0d0 + time))/(4.0d0*Dbar(1,2)*time) - time/tau)*&
           !                   (4.0d0*Dbar(1,2)*(1.0d0 + 4.0d0*Dbar(1,2)*dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI)*time +&
           !                   r**2*tau))/(64.0d0*Dbar(1,2)**3*M_PI**2*time**2*tau)
 
             ! for m1 \neq m2 \neq m
            !fluxdiv(i,j,1) = fluxdiv(i,j,1) + (dexp((-2.0d0*time)/tau - (r**2*(2.0d0 + time))/(4.0d0*Dbar(1,2)*time))*&
            !                  ((mass(1) - mass(2))*tau*(-2.0d0*r**2 + 4.0d0*Dbar(1,2)*time - r**2*time - 8.0d0*Dbar(1,2)*&
            !                  dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI*(r**2 - 2.0d0*Dbar(1,2)*time)) - 4.0d0*Dbar(1,2)*dexp(r**2/&
            !                  (4.0d0*Dbar(1,2)*time) + time/tau)*M_PI*time*(molmtot(i,j)*(1.0d0 + 4.0d0*Dbar(1,2)*dexp(r**2/&
            !                  (4.0d0*Dbar(1,2)))*M_PI)*(-(r**2*tau) + 4.0d0*Dbar(1,2)*time*(tau + time)) + mass(2)*tau*(r**2 -&
            !                  4.0d0*Dbar(1,2)*time + r**2*time + 4.0d0*Dbar(1,2)*dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI*(r**2 -&
            !                  4.0d0*Dbar(1,2)*time)))))/(256.0d0*Dbar(1,2)**4*molmtot(i,j)*M_PI**3*tau*time**4)

            !fluxdiv(i,j,2) = fluxdiv(i,j,2) + (dexp((-2.0d0*time)/tau - (r**2*(2.0d0 + time))/(4.0d0*Dbar(1,2)*time))*&
            !                 ((mass(1) - mass(2))*tau*(-2.0d0*r**2 + 4.0d0*Dbar(1,2)*time - r**2*time - 8.0d0*Dbar(1,2)*&
            !                 dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI*(r**2 - 2.0d0*Dbar(1,2)*time)) + 4.0d0*Dbar(1,2)*dexp(r**2/&
            !                 (4.0d0*Dbar(1,2)*time) + time/tau)*M_PI*time*(molmtot(i,j)*(1.0d0 + 4.0d0*Dbar(1,2)*dexp(r**2/&
            !                 (4.0d0*Dbar(1,2)))*M_PI)*(-(r**2*tau) + 4.0d0*Dbar(1,2)*time*(tau + time)) + mass(1)*tau*(r**2 -&
            !                 4.0d0*Dbar(1,2)*time + r**2*time + 4.0d0*Dbar(1,2)*dexp(r**2/(4.0d0*Dbar(1,2)))*M_PI*(r**2 -&
            !                 4.0d0*Dbar(1,2)*time)))))/(256.0d0*Dbar(1,2)**4*molmtot(i,j)*M_PI**3*tau*time**4)
           
            !print*, fluxdiv(i,j,1), fluxdiv(i,j,2)
 
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
    real(kind=dp_t)  :: x,y,z,rsq,tau,L(3)

    tau = 1.0d0
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
   
              rsq = (x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2
              fluxdiv(i,j,k,1) = fluxdiv(i,j,k,1) + (dexp(-rsq*(1.0d0+2.0d0*time)/(4.0d0*Dbar_in(1)*time)-&
                               time/tau)*(-0.0224484d0*Dbar_in(1)**21.0d0*(Dbar_in(1)*time)**3.5d0*dexp(3.0d0*&
                               rsq/(4.0d0*Dbar_in(1))) + 5.65621d0*10.0d0**(-6.0d0)*Dbar_in(1)**17.0d0*rsq*&
                               (Dbar_in(1)*time)**3.5d0*tau*dexp(rsq/(4.0d0*Dbar_in(1))) + 6.34864d0*10.0d0**&
                               (-8.0d0)*Dbar_in(1)**14.5d0*rsq*sqrt(rsq)*(Dbar_in(1)*time)**3.5d0*tau + &
                               Dbar_in(1)**19.5d0*dexp(rsq/(2.0d0*Dbar_in(1)))*(-0.00100786d0*(Dbar_in(1)*&
                               time)**3.5d0 - 0.000125983d0*rsq*(Dbar_in(1)*time)**2.5d0*tau)+ Dbar_in(1)**&
                               18.0d0*dexp(rsq/(4.0d0*Dbar_in(1)))*(-0.0000113124d0*(Dbar_in(1)*time)**3.5d0-&
                               2.82811d0*10.0d0**(-6.0d0)*rsq*(Dbar_in(1)*time)**2.5d0*tau)))/(Dbar_in(1)**&
                               19.5d0*(0.0224484d0 + Dbar_in(1)**1.5d0*dexp(rsq/(4.0d0*Dbar_in(1))))*(Dbar_in(1)*&
                               time)**5.0d0*tau)

              fluxdiv(i,j,k,2) = fluxdiv(i,j,k,2) + (dexp(-rsq*(1.0d0+2.0d0*time)/(4.0d0*Dbar_in(1)*time)-&
                               time/tau)*(0.0224484d0*Dbar_in(1)**21.0d0*(Dbar_in(1)*time)**3.5d0*dexp(3.0d0*&
                               rsq/(4.0d0*Dbar_in(1))) + Dbar_in(1)**14.5d0*rsq*sqrt(rsq)*(-6.34864d0*10.0d0**&
                               (-8.0d0)*(Dbar_in(1)*time)**3.5d0 + 2.82811d0*10.0d0**(-6.0d0)*dexp(rsq/(4.0d0*&
                               Dbar_in(1)*time)+time/tau)*(Dbar_in(1)*time)**5.0d0)*tau + Dbar_in(1)**17.0d0*rsq*&
                               (-5.65621d0*10.0d0**(-6.0d0)*(Dbar_in(1)*time)**3.5d0*dexp(rsq/(4.0d0*Dbar_in(1)))+&
                               0.000251965d0*(Dbar_in(1)*time)**5.0d0*dexp(rsq*(1.0d0+time)/(4.0d0*Dbar_in(1)*time)+&
                               time/tau))*tau + Dbar_in(1)**18.0d0*dexp(rsq/(4.0d0*Dbar_in(1)))*(0.0000113124d0*&
                               (Dbar_in(1)*time)**3.5d0 + 2.82811d0*10.0d0**(-6.0d0)*rsq*(Dbar_in(1)*time)**2.5d0*&
                               tau) + Dbar_in(1)**19.5d0*dexp(rsq/(2.0d0*Dbar_in(1)))*(0.00100786d0*(Dbar_in(1)*&
                               time)**3.5d0 + 0.000125983d0*rsq*(Dbar_in(1)*time)**2.5d0*tau)))/(Dbar_in(1)**&
                               19.5d0*(0.0224484d0 + Dbar_in(1)**1.5d0*dexp(rsq/(4.0d0*Dbar_in(1))))*(Dbar_in(1)*&
                               time)**5.0d0*tau)
 
           enddo
        enddo
     enddo

    end select 
 
  end subroutine external_source_3d

end module external_force_module 
