module init_solution_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use probin_module
  use bc_module
  use multifab_physbc_inhomogeneous_module

  use BoxLibRNGs

  implicit none

  private
  public :: init_solution

contains

  subroutine init_solution(mla,x_u,x_p,dx,time,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: nlevs,n,i,ng_u,ng_p,dm

    nlevs = mla%nlevel
    dm = mla%dim

    ng_u = x_u(1,1)%ng
    ng_p = x_p(1)%ng

    do n=1,nlevs
       do i = 1, nfabs(x_u(n,1))
          ump => dataptr(x_u(n,1),i)
          vmp => dataptr(x_u(n,2),i)
          pp  => dataptr(x_p(n),i)
          lo =  lwb(get_box(x_u(n,1),i))
          hi =  upb(get_box(x_u(n,1),i))
          select case (dm)
          case (2)
            call init_solution_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_u, &
                              pp(:,:,1,1), ng_p, lo, hi, dx(n,:), time)             
            
          case (3)
            wmp => dataptr(x_u(n,3),i)
            call init_solution_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                              pp(:,:,:,1), ng_p, lo, hi, dx(n,:), time)
            
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(x_u(n,i))
          ! It is also necessary to make sure the boundary values are consistent with the BCs:
          call multifab_physbc_domainvel_inhomogeneous(x_u(n,i),1,i,1,the_bc_level(n),dx(n,:))
       end do

    enddo

  end subroutine init_solution

  subroutine init_solution_2d(x_u,x_v,ng_u,x_p,ng_p,lo,hi,dx,time)

    integer, intent(in) :: lo(:), hi(:), ng_u, ng_p
    real(kind=dp_t), intent(inout) :: x_u(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: x_v(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: x_p(lo(1)-ng_p:,lo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y

    ! staggered x-velocity
    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + i*dx(1)
          x_u(i,j) = solution_val_2d(1,x,y)
       end do
    end do

    ! staggered y-velocity
    do j=lo(2),hi(2)+1
       y = prob_lo(2) + j*dx(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (i+0.5d0)*dx(1)
          x_v(i,j) = solution_val_2d(2,x,y)
       end do
    end do

    ! pressure
    do j=lo(2),hi(2)
       y = prob_lo(2) + (j+0.5d0)*dx(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (i+0.5d0)*dx(1)
          x_p(i,j) = solution_val_2d(3,x,y)
       end do
    end do

  contains
    
    function solution_val_2d(comp,x,y) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: velx,vely,Lx,Ly,freq,pfreq,vcst,ucst,ufac,y1,y2,visc_coef
      real(kind=dp_t) :: mean_value_u, mean_value_v, mean_value_p, random(1)

      visc_coef = coeff_mag(2)/coeff_mag(1) 

      select case ( abs(prob_sol) )
      case(1,2) ! Taylor-vortex to provide init sol and boudary condition

         velx = 1.d0
         vely = 1.d0
         Lx = prob_hi(1)-prob_lo(1)
         Ly = prob_hi(2)-prob_lo(2)

         if (comp .eq. 1) then
            ! staggered x-velocity
            val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
                 *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
                 *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
         else if (comp .eq. 3) then
            ! pressure
            val = -exp(-16.d0*M_PI**2*visc_coef*time/(Lx*Ly))* &
                 ( cos(4.d0*M_PI*(x-velx*time)/Lx) + cos(4.d0*M_PI*(y-vely*time)/Ly) )
         end if
         
      case(3) !  density=1 & smooth vis, with exact solution, strain form

         if (comp .eq. 1) then
            ! staggered x-velocity
            val = exp(-8.0d0*M_PI**2*time)*sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = exp(-8.0d0*M_PI**2*time)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*y)
         else if (comp .eq. 3) then
            ! pressure
            val = exp(-16.d0*M_PI**2*time)*(cos(4.d0*M_PI*x)+cos(4.d0*M_PI*y))
         end if

      case(4, 5, 20) ! 4: spherical bubble, 5: stripe, 20: Driven cavity 

         val = 0.d0

      case(6) ! see case (11) in varden-staggered/exact_solution

         ! incompressible modes with BCs       
         freq  = mode_coefs(1)*2.0*M_PI ! frequency in the periodic, "mod(prob_dir+1,DIM)" direction

         select case(mode_coefs(1))
         case(1)
            select case(mode_coefs(2))
               ! frequency in non-periodic, "prob_dir" direction
            case(1)
               pfreq  = 1.8583684196603974620d0 
            case(2)
               pfreq  = 5.4245857755826005234d0
            case(3)
               pfreq  = 8.8049901634565292604d0
            case default
               call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_sol=6 in init_solution_2d')
            end select

         case(2)
            select case(mode_coefs(2))
               ! frequency in non-periodic, "prob_dir" direction
            case(1)
               pfreq  = 2.1789521056066303662d0
            case(2)
               pfreq  = 5.7874274827134140679d0
            case(3)
               pfreq  = 9.0932821164949254163d0
            case default
               call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_sol=6 in init_solution_2d')
            end select

         case default
            call bl_error('ERROR: Only mode_coefs(1)<=2 supported for prob_sol=6 in init_solution_2d')
         end select

         ucst = visc_coef*(freq*freq+pfreq*pfreq) 
         ! since this is an exact solution for the time-dependant Stokes equations, 
         ! we multiply by a small prefactor so that the nonlinearity is negligable
         ufac = (1E-8)*dexp(-ucst*time)

         if (prob_dir .eq. 1) then

            if (comp .eq. 1) then
               val = ufac*sin(freq*y)/cos(pfreq)*(cos(pfreq*x)*cosh(freq) &
                    -cosh(freq*x)*cos(pfreq)) 
            else if (comp .eq. 2) then
               val = ufac*cos(freq*y)/sin(pfreq)*(sin(pfreq*x)*sinh(freq) &
                    -sinh(freq*x)*sin(pfreq)) 
            end if

         else if (prob_dir .eq. 2) then

            if (comp .eq. 1) then
               val = ufac*cos(freq*x)/sin(pfreq)*(sin(pfreq*y)*sinh(freq) &
                    -sinh(freq*y)*sin(pfreq)) 
            else if (comp .eq. 2) then
               val = ufac*sin(freq*x)/cos(pfreq)*(cos(pfreq*y)*cosh(freq) &
                    -cosh(freq*y)*cos(pfreq)) 
            end if

         else
            call bl_error('init_solution.f90: invalid prob_dir for prob_sol=6')
         end if

         ! We no longer maintain pressure, but here it is for the record:
         if(.false.) then
            if (comp .eq. 3) then
               ! this may be wrong--check pfac!?!
               val = -ufac*sinh(freq*x)*sin(freq*y)/freq
            end if
         end if

         call bl_error('solution_val_2d: pressure not defined for prob_sol=6')

      case (31,32) ! for testing discretization accuracy 
         if (comp .eq. 1) then
            ! staggered x-velocity
            val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
         else if (comp .eq. 3) then
            ! pressure
            val = cos(4.d0*M_PI*x)+cos(4.d0*M_PI*y)
         end if

      case(100) ! random initial solution

         call UniformRNGs(random,size(random))
         if (.false. .and. (comp .eq. 3) ) then ! Testing of scaling problems
            ! scale pressure by eta/dx
            val = random(1) * coeff_mag(2)/dx(1)            
         else
            val = random(1)
         end if

      case default
         call bl_error('ERROR: bad or unimplemented choice for prob_sol in solution_val_2d')
      end select

    end function solution_val_2d

  end subroutine init_solution_2d

  subroutine init_solution_3d(x_u,x_v,x_w,ng_u,x_p,ng_p,lo,hi,dx,time)

    integer, intent(in) :: lo(:), hi(:), ng_u, ng_p
    real(kind=dp_t), intent(inout) :: x_u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: x_v(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: x_w(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: x_p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    ! staggered x-velocity
    do k=lo(3),hi(3)
       z = prob_lo(3) + (k+0.5d0)*dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + i*dx(1)
             x_u(i,j,k) = solution_val_3d(1,x,y,z)
          end do
       end do
    end do

    ! staggered y-velocity
    do k=lo(3),hi(3)
       z = prob_lo(3) + (k+0.5d0)*dx(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + j*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (i+0.5d0)*dx(1)
             x_v(i,j,k) = solution_val_3d(2,x,y,z)
          end do
       end do
    end do

    ! staggered w-velocity
    do k=lo(3),hi(3)+1
       z = prob_lo(3) + k*dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (i+0.5d0)*dx(1)
             x_w(i,j,k) = solution_val_3d(3,x,y,z)
          enddo
       enddo
    enddo

    ! pressure
    do k=lo(3),hi(3)
       z = prob_lo(3) + (k+0.5d0)*dx(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (i+0.5d0)*dx(1)
             x_p(i,j,k) = solution_val_3d(4,x,y,z)
          end do
       end do
    end do
    
  contains

    function solution_val_3d(comp,x,y,z) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y,z
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: velx,vely,Lx,Ly,freq,pfreq,vcst,ucst,y1,y2,visc_coef
      real(kind=dp_t) :: mean_value_u, mean_value_v, mean_value_p
      real(kind=dp_t) :: ufac,vfac,pfac,hx,hy,hz,random(1)

      visc_coef = coeff_mag(2)/coeff_mag(1) 

      select case ( abs(prob_sol) )
      case(1) ! Taylor-vortex, quasi-2D

         velx = 1.d0
         vely = 1.d0
         Lx = prob_hi(1)-prob_lo(1)
         Ly = prob_hi(2)-prob_lo(2)

         if (comp .eq. 1) then
            ! staggered x-velocity
            val = velx - 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
                 *cos(2.d0*M_PI*(x-velx*time)/Lx)*sin(2.d0*M_PI*(y-vely*time)/Ly)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = vely + 2.d0*exp(-8.d0*M_PI**2*visc_coef*time/(Lx*Ly)) &
                 *sin(2.d0*M_PI*(x-velx*time)/Lx)*cos(2.d0*M_PI*(y-vely*time)/Ly)
         else if (comp .eq. 3) then
            ! staggered w-velocity
            val = 0.d0
         else if (comp .eq. 4) then
            ! pressure
            val = -exp(-16.d0*M_PI**2*visc_coef*time/(Lx*Ly))* &
                 ( cos(4.d0*M_PI*(x-velx*time)/Lx) + cos(4.d0*M_PI*(y-vely*time)/Ly) )
         end if

      case(2)    ! ABC flow    

         freq  = 2.d0*M_PI
         ufac = dexp(-freq*freq*visc_coef*time)
         pfac = dexp(-2.0d0*freq*freq*visc_coef*time)
         hx = ABC_coefs(1)
         hy = ABC_coefs(2)
         hz = ABC_coefs(3)

         if (comp .eq. 1) then
            val = 1.0d0 + ufac*(hz*cos(freq*(y-time))+hx*sin(freq*(z-time)))
         else if (comp .eq. 2) then
            val = 1.0d0 + ufac*(hy*sin(freq*(x-time))+hx*cos(freq*(z-time)))
         else if (comp .eq. 3) then
            val = 1.0d0 + ufac*(hy*cos(freq*(x-time))+hz*sin(freq*(y-time)))
         else if (comp .eq. 4) then
            val = -pfac*(hx*hz*cos(freq*(y-time))*sin(freq*(z-time))+hx*hy*sin(freq*(x-time))*cos(freq*(z-time))&
                 +hy*hz*cos(freq*(x-time))*sin(freq*(y-time)))
         end if

      case(3) ! smooth density & vis, with exact solution, strain form
         
         if (comp .eq. 1) then
            ! staggered x-velocity
            val = exp(-8.0d0*M_PI**2*time)*sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = exp(-8.0d0*M_PI**2*time)*cos(2.d0*M_PI*x)*cos(2.d0*M_PI*y)
         else if (comp .eq. 3) then
            ! staggered w-velocity
            val = 0.d0
         else if (comp .eq. 4) then
            ! pressure
            val = exp(-16.d0*M_PI**2*time)*(cos(4.d0*M_PI*x)+cos(4.d0*M_PI*y))
         end if

      case(4,5,20) ! 4: spherical bubble, 5: stripe, 20: Driven cavity

         val = 0.d0

      case(6)
         ! incompressible modes with BCs
         ! Only quasi-2D flow supported for now

         freq  = mode_coefs(1)*2.0*M_PI ! frequency in the periodic, "mod(prob_dir+1,DIM)" direction

         select case(mode_coefs(1))
         case(1)
            select case(mode_coefs(2))
               ! frequency in non-periodic, "prob_dir" direction
            case(1)
               pfreq  = 1.8583684196603974620d0 
            case(2)
               pfreq  = 5.4245857755826005234d0
            case(3)
               pfreq  = 8.8049901634565292604d0
            case default
               call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_sol=6 in init_solution_2d')
            end select
         case(2)
            select case(mode_coefs(2))
               ! frequency in non-periodic, "prob_dir" direction
            case(1)
               pfreq  = 2.1789521056066303662d0
            case(2)
               pfreq  = 5.7874274827134140679d0
            case(3)
               pfreq  = 9.0932821164949254163d0
            case default
               call bl_error('ERROR: Only mode_coefs(2)<=3 supported for prob_sol=6 in init_solution_2d')
            end select
         case default
            call bl_error('ERROR: Only mode_coefs(1)<=2 supported for prob_sol=6 in init_solution_2d')
         end select

         ucst = visc_coef * (freq*freq + pfreq*pfreq)
         ! since this is an exact solution for the time-dependant Stokes equations, 
         ! we multiply by a small prefactor so that the nonlinearity is negligable
         ufac = (1E-8)*dexp(-ucst*time)

         if (prob_dir .eq. 1) then

            if (comp .eq. 1) then
               val = ufac*sin(freq*y)/cos(pfreq)*(cos(pfreq*x)*cosh(freq) &
                    -cosh(freq*x)*cos(pfreq)) 
            else if (comp .eq. 2) then
               val = ufac*cos(freq*y)/sin(pfreq)*(sin(pfreq*x)*sinh(freq) &
                    -sinh(freq*x)*sin(pfreq)) 
            else if (comp .eq. 3) then
               val = 0.d0
            end if

         else if (prob_dir .eq. 2) then

            if (comp .eq. 1) then
               val = 0.d0
            else if (comp .eq. 2) then
               val = ufac*sin(freq*z)/cos(pfreq)*(cos(pfreq*y)*cosh(freq) &
                    -cosh(freq*y)*cos(pfreq))
            else if (comp .eq. 3) then
               val = ufac*cos(freq*z)/sin(pfreq)*(sin(pfreq*y)*sinh(freq) &
                    -sinh(freq*y)*sin(pfreq)) 
            end if

         else if (prob_dir .eq. 3) then

            if (comp .eq. 1) then
               val = ufac*cos(freq*x)/sin(pfreq)*(sin(pfreq*z)*sinh(freq) &
                    -sinh(freq*z)*sin(pfreq)) 
            else if (comp .eq. 2) then
               val = 0.d0
            else if (comp .eq. 3) then
               val = ufac*sin(freq*x)/cos(pfreq)*(cos(pfreq*z)*cosh(freq) &
                    -cosh(freq*z)*cos(pfreq))
            end if

         else
            call bl_error('init_solutions.f90: invalid prob_dir for prob_sol=6')
         end if

         call bl_error('solution_val_3d: pressure not defined for prob_sol=6')

      case (31,32) ! for testing discretization accuracy 
         if (comp .eq. 1) then
            ! staggered x-velocity
            val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
         else if (comp .eq. 2) then
            ! staggered y-velocity
            val = cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)
         else if (comp .eq. 3) then
            ! staggered z-velocity
            val = sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*cos(2.d0*M_PI*z)
         else if (comp .eq. 4) then
            ! pressure
            val = cos(4.d0*M_PI*x)+cos(4.d0*M_PI*y)+cos(4.d0*M_PI*z)
         end if

      case(100) ! random initial solution

         call UniformRNGs(random,size(random))
         val = random(1)

      case default
         call bl_error('ERROR: bad or unimplemented choice for prob_sol in solution_val_3d')
      end select

    end function solution_val_3d

  end subroutine init_solution_3d

end module init_solution_module
