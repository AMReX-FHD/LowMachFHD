module correction_flux_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use probin_multispecies_module

  implicit none

  private

  public :: correction_flux

contains

  subroutine correction_flux(mla,rho,rho_tot,flux,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: rho_tot(:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: n,i,dm,nlevs,ng_p,ng_g
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointer for rho, rho_tot, flux_x, flux_y and flux_z 
    real(kind=dp_t), pointer :: dp(:,:,:,:)      ! for rho
    real(kind=dp_t), pointer :: dp1(:,:,:,:)     ! for rho_tot
    real(kind=dp_t), pointer :: flux_x(:,:,:,:)
    real(kind=dp_t), pointer :: flux_y(:,:,:,:)
    real(kind=dp_t), pointer :: flux_z(:,:,:,:)

    dm    = mla%dim       ! dimensionality
    nlevs = mla%nlevel    ! number of levels 
    ng_p  = rho(1)%ng     ! number of ghost cells for rho
    ng_g  = flux(1,1)%ng  ! number of ghost cells for flux

    ! loop over all boxes   
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp  => dataptr(rho(n), i)
          dp1 => dataptr(rho_tot(n), i)
          flux_x => dataptr(flux(n,1), i)
          flux_y => dataptr(flux(n,2), i)
          lo = lwb(get_box(rho(n), i))
          hi = upb(get_box(rho(n), i))
          
          select case (dm)
          case (2)
          call correction_flux_2d(dp(:,:,1,:),dp1(:,:,1,1),ng_p,flux_x(:,:,1,:),flux_y(:,:,1,:),&
                                  ng_g,lo,hi)
          case (3)
          flux_z => dataptr(flux(n,3), i)
          call correction_flux_3d(dp(:,:,:,:),dp1(:,:,:,1),ng_p,flux_x(:,:,:,:),flux_y(:,:,:,:),&
                                  flux_z(:,:,:,:),ng_g,lo,hi)
          end select
       enddo
    enddo

  contains
    
    subroutine correction_flux_2d(rho,rho_tot,ng_p,flux_x,flux_y,ng_g,lo,hi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_p:,lo(2)-ng_p:,:)
      real(kind=dp_t), intent(in   ) :: rho_tot(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) :: flux_x(lo(1)-ng_g:,lo(2)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_y(lo(1)-ng_g:,lo(2)-ng_g:,:)

      ! local
      integer         :: i,j,n
      real(kind=dp_t) :: sumx,sumy,corr

      ! x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
               
            ! free the data
            sumx = 0.d0
            corr = 0.d0
              
            ! sum the x-fluxes upto nspecies-1 
            do n=1, nspecies-1
               sumx = sumx + flux_x(i,j,n)
            enddo
              
            ! caculate corr and print error if not zero 
            corr = flux_x(i,j,nspecies) + sumx
            if(corr .gt. fraction_tolerance) print*, "Sum of x-flux=", corr
 
            if(corr .gt. rho_tot(i,j)*1e-8) then
               write(*,*) "Error: sum of x-fluxes greater than rho_tot*1e-8"
               write(*,*) "sum is",corr
            endif
              
            ! correct x-flux for last species  
            flux_x(i,j,nspecies) = -sumx             

         enddo
      enddo
   
      ! y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

            ! free the data
            sumy  = 0.d0
            corr = 0.d0
              
            ! sum the y-fluxes upto nspecies-1 
            do n=1, nspecies-1
               sumy = sumy + flux_y(i,j,n)
            enddo
              
            ! caculate corr and print error if not zero 
            corr = flux_y(i,j,nspecies) + sumy
            if(corr .gt. fraction_tolerance) print*, "Sum of y-flux=", corr
            
            if(corr .gt. rho_tot(i,j)*1e-8) then
               write(*,*) "Error: sum of y-fluxes greater than rho_tot*1e-8"
               write(*,*) "sum is",corr         
            endif
              
            ! correct y-flux for last species  
            flux_y(i,j,nspecies) = -sumy             
 
         enddo
      enddo

    end subroutine correction_flux_2d

    subroutine correction_flux_3d(rho,rho_tot,ng_p,flux_x,flux_y,flux_z,ng_g,lo,hi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
      real(kind=dp_t), intent(in   ) :: rho_tot(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) :: flux_x(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_y(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_z(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)

      ! local
      integer         :: i,j,k,n
      real(kind=dp_t) :: sumx,sumy,sumz,corr
      
      ! x-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
   
               ! free the data
               sumx = 0.d0
               corr = 0.d0
              
               ! sum the x-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumx = sumx + flux_x(i,j,k,n)
               enddo
              
               ! caculate corr and print error if not zero 
               corr = flux_x(i,j,k,nspecies) + sumx 
              
               if(corr .gt. fraction_tolerance) print*, "Sum of x-flux=", corr
               if(corr .gt. rho_tot(i,j,k)*1e-8) then
                  write(*,*) "Error: sum of x-fluxes greater than rho_tot*1e-8"             
                  write(*,*) "sum is",corr
               endif
              
               ! correct x-flux for last species  
               if(abs(flux_x(i,j,k,nspecies)).gt.abs(sumx) .or. abs(flux_x(i,j,k,nspecies)) & 
                  .lt.abs(sumx)) then 
                  flux_x(i,j,k,nspecies) = -sumx             
               endif

            enddo
         enddo
      enddo

      ! y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
              
               ! free the data
               sumy  = 0.d0
               corr = 0.d0
              
               ! sum the y-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumy = sumy + flux_y(i,j,k,n)
               enddo
              
               ! caculate corr and print error if not zero 
               corr = flux_y(i,j,k,nspecies) + sumy 
               
               if(corr .gt. fraction_tolerance) print*, "Sum of y-flux=", corr
               if(corr .gt. rho_tot(i,j,k)*1e-8) then
                  write(*,*) "Error: sum of y-fluxes greater than rho_tot*1e-8"             
                  write(*,*) "sum is",corr
               endif
              
               ! correct y-flux for last species  
               if(abs(flux_y(i,j,k,nspecies)).gt.abs(sumy) .or. abs(flux_y(i,j,k,nspecies)) & 
                  .lt.abs(sumy)) then 
                  flux_y(i,j,k,nspecies) = -sumy             
               endif
 
            enddo
         enddo
      enddo

      ! z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
            
               ! free the data
               sumz  = 0.d0
               corr = 0.d0
              
               ! sum the z-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumz = sumz + flux_z(i,j,k,n)
               enddo
              
               ! caculate corr and print error if not zero 
               corr = flux_z(i,j,k,nspecies) + sumz 
               
               if(corr .gt. fraction_tolerance) print*, "Sum of z-flux=", corr
               if(corr .gt. rho_tot(i,j,k)*1e-8) then
                  write(*,*) "Error: sum of z-fluxes greater than rho_tot*1e-8"             
                  write(*,*) "sum is",corr
               endif
              
               ! correct z-flux for last species  
               if(abs(flux_z(i,j,k,nspecies)).gt.abs(sumz) .or. abs(flux_z(i,j,k,nspecies)) & 
                  .lt.abs(sumz)) then 
                  flux_z(i,j,k,nspecies) = -sumz             
               endif
 
            enddo
         enddo
      enddo

    end subroutine correction_flux_3d

  end subroutine correction_flux

end module correction_flux_module
