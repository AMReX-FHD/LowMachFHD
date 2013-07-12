module init_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use BoxLibRNGs
  use multifab_physbc_stag_module
  use bc_module

  implicit none

  private

  public :: init_solution, init_rhs, init_mat

contains

  subroutine init_solution(mla,x_u,x_p,dx,time,the_bc_level,vel_bc_n,vel_bc_t)

    use probin_module       , only: prob_dir, prob_sol, coeff_mag, ABC_coefs, mode_coefs
    use probin_common_module, only: prob_lo, prob_hi

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(in   ) :: vel_bc_n(:,:)
    type(multifab) , intent(in   ) :: vel_bc_t(:,:)

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
          call multifab_physbc_domainvel(x_u(n,i),i,the_bc_level(n),dx(n,:),vel_bc_n(n,:))
       end do

    enddo

  contains

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

    end subroutine init_solution_3d

    function solution_val_2d(comp,x,y) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: velx,vely,Lx,Ly,freq,pfreq,vcst,ucst,ufac,y1,y2,visc_coef
      real(kind=dp_t) :: mean_value_u, mean_value_v, mean_value_p

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

         call UniformRNG(val)

      case default
         call bl_error('ERROR: bad or unimplemented choice for prob_sol in solution_val_2d')
      end select

    end function solution_val_2d

    function solution_val_3d(comp,x,y,z) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y,z
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: velx,vely,Lx,Ly,freq,pfreq,vcst,ucst,y1,y2,visc_coef
      real(kind=dp_t) :: mean_value_u, mean_value_v, mean_value_p
      real(kind=dp_t) :: ufac,vfac,pfac,hx,hy,hz

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
            val = -pfac*( hx*hz*cos(freq*(y-time))*sin(freq*(z-time)) &
                         +hx*hy*sin(freq*(x-time))*cos(freq*(z-time))&
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

         call UniformRNG(val)

      case default
         call bl_error('ERROR: bad or unimplemented choice for prob_sol in solution_val_3d')
      end select

    end function solution_val_3d

  end subroutine init_solution

  subroutine init_rhs(mla,b_u,b_p,dx,time,the_bc_level)

    use probin_module       , only: prob_coeff, prob_sol, coeff_mag, ABC_coefs, theta_fac
    use probin_common_module, only: fixed_dt, prob_lo, prob_hi, visc_type

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: b_u(:,:)
    type(multifab) , intent(inout) :: b_p(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: uxp(:,:,:,:), uyp(:,:,:,:), uzp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: n,nlevs,i,dm,ng_u,ng_p

    nlevs = mla%nlevel
    dm = mla%dim

    ng_u = b_u(1,1)%ng
    ng_p = b_p(1)%ng

    if (ng_u .ne. 0) then
       call bl_error("init_rhs expects b_u to have 0 ghost cells")
    end if

    if (ng_p .ne. 0) then
       call bl_error("init_rhs expects b_p to have 0 ghost cells")
    end if

    do n=1,nlevs
       do i=1,nfabs(b_u(n,1))
          uxp => dataptr(b_u(n,1),i)
          uyp => dataptr(b_u(n,2),i)
          pp  => dataptr(b_p(n),i)
          lo = lwb(get_box(b_u(n,1),i))
          hi = upb(get_box(b_u(n,1),i))
          select case (dm)
          case (2)
             call init_rhs_2d(uxp(:,:,1,1), uyp(:,:,1,1), ng_u, &
                              pp(:,:,1,1), ng_p, &
                              lo, hi, dx(n,:), time)

          case (3)
             uzp => dataptr(b_u(n,3),i)
             call init_rhs_3d(uxp(:,:,:,1), uyp(:,:,:,1), uzp(:,:,:,1), ng_u, &
                              pp(:,:,:,1), ng_p, &
                              lo, hi, dx(n,:), time)
          end select
       end do

       ! For periodic boundaries, ensure the low and high side are consistent:
       ! Note: multifab_internal_sync compares the box number of the two boxes 
       ! with overlapping values and the data on the box with lower number wins. 
       do i=1,dm
          call multifab_internal_sync(b_u(n,i))
          ! This will set the forcing on Dirichlet boundaries to zero:
          call multifab_physbc_domainvel(b_u(n,i),i,the_bc_level(n),dx(n,:))
       end do

    enddo

  contains

    subroutine init_rhs_2d(b_ux,b_uy,ng_u,b_p,ng_p,lo,hi,dx,time)

      use bl_constants_module

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p
      real(kind=dp_t), intent(inout) :: b_ux(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(inout) :: b_uy(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(inout) ::  b_p(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      ! local 
      integer :: i,j
      real(kind=dp_t) :: x,y

      ! staggered x-component
      do j=lo(2),hi(2)
         y = prob_lo(2) + (j+0.5d0)*dx(2)
         do i=lo(1),hi(1)+1
            x = prob_lo(1) + i*dx(1)  
            b_ux(i,j) = rhs_val_2d(1,x,y)
         end do
      end do

      ! staggered y-component
      do j=lo(2),hi(2)+1
         y = prob_lo(2) + j*dx(2)
         do i=lo(1),hi(1)
            x = prob_lo(1) + (i+0.5d0)*dx(1)
            b_uy(i,j) = rhs_val_2d(2,x,y)
         end do
      end do

      ! cell-centered pressure
      do j=lo(2),hi(2)
         y = prob_lo(2) + (j+0.5d0)*dx(2)
         do i=lo(1),hi(1)
            x = prob_lo(1) + (i+0.5d0)*dx(1)
            b_p(i,j) = rhs_val_2d(3,x,y)
         end do
      end do

    end subroutine init_rhs_2d

    subroutine init_rhs_3d(b_ux,b_uy,b_uz,ng_u,b_p,ng_p,lo,hi,dx,time)

      use bl_constants_module

      integer        , intent(in   ) :: lo(:), hi(:), ng_u, ng_p
      real(kind=dp_t), intent(inout) :: b_ux(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) :: b_uy(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) :: b_uz(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) ::  b_p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      ! local 
      real(kind=dp_t) :: x, y, z
      integer :: i, j, k

      ! staggered x-component
      do k=lo(3),hi(3)
         z = prob_lo(3) + (k+0.5d0)*dx(3)
         do j=lo(2),hi(2)
            y = prob_lo(2) + (j+0.5d0)*dx(2)
            do i=lo(1),hi(1)+1
               x = prob_lo(1) + i*dx(1)
               b_ux(i,j,k) = rhs_val_3d(1,x,y,z)
            end do
         end do
      end do

      ! staggered y-component
      do k=lo(3),hi(3)
         z = prob_lo(3) + (k+0.5d0)*dx(3)
         do j=lo(2),hi(2)+1
            y = prob_lo(2) + j*dx(2)
            do i=lo(1),hi(1)
               x = prob_lo(1) + (i+0.5d0)*dx(1)
               b_uy(i,j,k) = rhs_val_3d(2,x,y,z)
            end do
         end do
      end do

      ! staggered z-component
      do k=lo(3),hi(3)+1
         z = prob_lo(3) + k*dx(3)
         do j=lo(2),hi(2)
            y = prob_lo(2) + (j+0.5d0)*dx(2)
            do i=lo(1),hi(1)
               x = prob_lo(1) + (i+0.5d0)*dx(1)
               b_uz(i,j,k) = rhs_val_3d(3,x,y,z)
            end do
         end do
      end do

      ! cell-centered pressure
      do k=lo(3),hi(3)
         z = prob_lo(3) + (k+0.5d0)*dx(3)
         do j=lo(2),hi(2)
            y = prob_lo(2) + (j+0.5d0)*dx(2)
            do i=lo(1),hi(1)
               x = prob_lo(1) + (i+0.5d0)*dx(1)
               b_p(i,j,k) = rhs_val_3d(4,x,y,z)
            end do
         end do
      end do

    end subroutine init_rhs_3d

    function rhs_val_2d(comp,x,y) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: visc_coef, L, freq

      visc_coef=coeff_mag(2)/coeff_mag(1)

      if( (abs(prob_sol)==1) .and. (abs(prob_coeff)==1) .and. (theta_fac==0.0d0) ) then
         ! Constant-coefficient Taylor vortex

         L = prob_hi(1)-prob_lo(1)
         freq  = 2.d0*M_PI

         select case(abs(visc_type))      ! which visc_type is using?
         case(1,2,3)      ! as Taylor vortex is divergence free, visc_type=1, or 2 or 3 no difference

            if (comp .eq. 1) then
               ! staggered x-component
               val=-4.d0*M_PI*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*(-4.d0*M_PI*visc_coef*cos(freq*(-x+time)/L)*&
                    sin(freq*(-y+time)/L)+exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(4.d0*M_PI*(-x+time)/L)*L)/L**2     
               ! if only viscous term is concerned, uncomment the following code 
               !val=-16.d0*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*cos(freq*(x-time)/L)*M_PI**2*sin(freq*(y-time)/L)/L**2 
            else if (comp .eq. 2) then
               ! staggered y-component
               val=-4.d0*M_PI*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*(4.d0*M_PI*visc_coef*sin(freq*(-x+time)/L)*&
                    cos(freq*(-y+time)/L)+exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(4.d0*M_PI*(-y+time)/L)*L)/L**2
               ! if only viscous term is concerned, uncomment the following code 
               !val=16.d0*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(freq*(x-time)/L)*M_PI**2*cos(freq*(y-time)/L)/L**2
            else if (comp .eq. 3) then
               ! cell-centered pressure
               val = 0.d0
            end if

         case default
            call bl_error('ERROR:rhs_val_2d, bad or unimplemented choice for visc_type')
         end select

      elseif ( (abs(prob_sol)==2) .or. (abs(prob_sol)==4) .or. (abs(prob_sol)==6) .or. (abs(prob_sol)==7) ) then
         ! RHS is 0 

         val = 0.d0 

      elseif ( (abs(prob_sol)==3) .and. (abs(prob_coeff)==5) .and. (theta_fac==0.0d0) ) then  
         ! density=1, smooth viscosity flow

         freq = 2.d0*M_PI
         if (abs(visc_type)==2) then      ! visc_type has to be -2 

            if (comp .eq. 1) then
               ! staggered x-component
               val = -4.0d0*M_PI*exp(-8.0d0*M_PI**2*time)*(-2.0d0*sin(freq*x)*M_PI*cos(freq*x)+&
                    2.0d0*sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2+exp(-8.0d0*M_PI**2*time)*sin(4.0d0*M_PI*x))
            else if (comp .eq. 2) then
               ! staggered y-component
               val = 4.0d0*M_PI*exp(-8.0d0*M_PI**2*time)*(2.0d0*cos(freq*x)**2*cos(freq*y)*M_PI*&
                    sin(freq*y)-exp(-8.0d0*M_PI**2*time)*sin(4.0d0*M_PI*y))
            else if (comp .eq. 3) then
               ! cell-centered pressure
               val = 0.d0
            end if

         else 
            call bl_error('ERROR:rhs_val_2d, bad choice for visc_type')
         end if

      elseif ( (abs(prob_sol)==31) .and. ((prob_coeff==5) .or. (prob_coeff==1)) ) then 
         ! smoothly variable viscosity or vis=1, divergence not free solution, for checking discretization accuracy

         freq  = 2.d0*M_PI

         select case(abs(visc_type))  
         case(1)

            if (comp .eq. 1) then
               ! staggered x-component
               if (prob_coeff .eq. 5) then 
                  val = -(-2.d0*sin(freq*x)*(3.d0*cos(freq*x)-4.d0*cos(freq*x)*&
                       cos(freq*y)**2+4.d0*sin(freq*y))*M_PI**2) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)     ! last line is dp/dx
               else   
                  ! by default prob_coeff=1
                  val = 8.d0*sin(freq*x)*M_PI**2*sin(freq*y) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)     ! last line is dp/dx
               end if
            else if (comp .eq. 2) then
               ! staggered y-component  ! for testing steady state problem
               if (prob_coeff .eq. 5) then 
                  val = 2.d0*M_PI**2*(-1.d0+cos(freq*y)**2+3.d0*cos(freq*x)**2-&
                       4.d0*cos(freq*x)**2*cos(freq*y)**2+4.d0*cos(freq*x)*sin(freq*y)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)     ! last line is dp/dy
               else   
                  ! by default prob_coeff=1
                  val = 8.d0*cos(freq*x)*M_PI**2*sin(freq*y) &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)     ! last line is dp/dy
               end if
            end if

         case(2)

            if (comp .eq. 1) then
               ! staggered x-component
               if (prob_coeff .eq. 5) then 
                  val = -(-2.d0*sin(freq*x)*(5.d0*cos(freq*x)-6.d0*cos(freq*x)*cos(freq*y)**2+&
                       6.d0*sin(freq*y)+2.d0*cos(freq*x)*cos(freq*y)*sin(freq*y)+2.d0*cos(freq*y))*M_PI**2) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)     ! last line is dp/dx
               else  
                  ! by default prob_coeff=1
                  val = 4.d0*sin(freq*x)*M_PI**2*(3.d0*sin(freq*y)+cos(freq*y)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)     ! last line is dp/dx
               end if
            else if (comp .eq. 2) then
               ! staggered y-component  ! for testing steady state problem
               if (prob_coeff .eq. 5) then 
                  val = -(2.d0*M_PI**2*(6.d0*cos(freq*x)**2*cos(freq*y)**2-6.d0*cos(freq*x)*sin(freq*y)-&
                       4.d0*cos(freq*x)**2-sin(freq*y)*cos(freq*y)+2.d0*cos(freq*x)**2*cos(freq*y)*sin(freq*y)&
                       +2.d0*cos(freq*x)*cos(freq*y)+1.d0-cos(freq*y)**2)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)     ! last line is dp/dy
               else   
                  ! by default prob_coeff=1
                  val = -(-12.d0*cos(freq*x)*M_PI**2*sin(freq*y)+&
                       4.d0*cos(freq*x)*M_PI**2*cos(freq*y)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)     ! last line is dp/dy
               end if
            end if

         case(3) ! gamma expression refer to init_mat

            if (comp .eq. 1) then
               ! staggered x-component
               if (prob_coeff .eq. 5) then 
                  val=-(2.d0/3.d0)*M_PI**2*(-11.d0*sin(freq*x)*cos(freq*x)+14.d0*cos(freq*x)*cos(freq*y)**2&
                       *sin(freq*x)-14.d0*sin(freq*x)*sin(freq*y)-2.d0*cos(freq*x)*cos(freq*y)*sin(freq*x)*&
                       sin(freq*y)-2.d0*sin(freq*x)*cos(freq*y)+12.d0*cos(freq*x)**2-12.d0*cos(freq*x)**2*&
                       cos(freq*y)**2+12.d0*cos(freq*x)**2*sin(freq*y)*cos(freq*y)-6.d0+6.d0*cos(freq*y)**2&
                       -6.d0*sin(freq*y)*cos(freq*y)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)     ! this line is dp/dx
               else 
                  ! by default prob_coeff=1, note that div not free
                  val=-(-4.d0*sin(freq*x)*M_PI**2*(3.d0*sin(freq*y)+cos(freq*y))) &
                       -4.d0*M_PI*sin(4.d0*M_PI*x)    &    ! the line is dp/dx , alpha=gamma=1, div not free      
                       +(4.d0/3.d0)*sin(freq*x)*M_PI**2*sin(freq*y)+(4.d0/3.d0)*sin(freq*x)*M_PI**2*cos(freq*y)
               end if
            else if (comp .eq. 2) then
               ! staggered y-component  ! for testing steady state problem
               if (prob_coeff .eq. 5) then 
                  val=-(2.d0/3.d0)*M_PI**2*(14.d0*cos(freq*x)**2*cos(freq*y)**2-14.d0*cos(freq*x)*&
                       sin(freq*y)-10.d0*cos(freq*x)**2-3.d0*sin(freq*y)*cos(freq*y)+2.d0*cos(freq*x)**2*&
                       sin(freq*y)*cos(freq*y)+2.d0*cos(freq*x)*cos(freq*y)+3.d0-3.d0*cos(freq*y)**2+&
                       12.d0*cos(freq*x)*cos(freq*y)*sin(freq*x)*sin(freq*y)+12.d0*cos(freq*x)*&
                       cos(freq*y)**2*sin(freq*x)-6.d0*sin(freq*x)*cos(freq*x))  &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)       ! this line is dp/dy   
               else
                  ! by default prob_coeff=1, note that div not free
                  val=-(-12.d0*cos(freq*x)*M_PI**2*sin(freq*y)+&
                       4.d0*cos(freq*x)*M_PI**2*cos(freq*y)) &
                       -4.d0*M_PI*sin(4.d0*M_PI*y)  &   ! the line is dp/dy, alpha=gamma=1, div not free
                       -(4.d0/3.d0)*cos(freq*x)*M_PI**2*cos(freq*y)+(4.d0/3.d0)*cos(freq*x)*M_PI**2*sin(freq*y)
               end if
            end if

         case default
            call bl_error('ERROR:init_rhs_2d, bad or unimplemented choice for prob_sol for testing accuracy')
         end select

         if (comp .eq. 3) then
            ! cell-centered pressure
            ! b_p=-Div, not divergence free
            val = -(2.d0*cos(freq*x)*M_PI*sin(freq*y)+2.d0*cos(freq*x)*cos(freq*y)*M_PI)
         end if

      elseif ( (abs(prob_sol)==32) .and. (prob_coeff ==8) ) then 
         ! smoothly variable density, for checking discretization accuracy of smooth density elliptic operator 

         freq  = 2.d0*M_PI
         select case(abs(visc_type))  
         case(1,2,3)
            if (comp .eq. 1) then
               ! staggered x-component              
               val = (sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)/fixed_dt)*(1.0d0+0.5d0*cos(freq*x)*sin(freq*y))&
                    -4.d0*M_PI*sin(4.d0*M_PI*x)

            else if (comp .eq. 2) then
               ! staggered y-component  ! for testing steady state problem
               val = (cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)/fixed_dt)*(1.0d0+0.5d0*cos(freq*x)*sin(freq*y))&
                    -4.d0*M_PI*sin(4.d0*M_PI*y)
            end if

         case default
            call bl_error('ERROR:init_rhs_2d, bad or unimplemented choice for prob_sol for testing accuracy')
         end select

         if (comp .eq. 3) then
            ! cell-centered pressure
            ! b_p=-Div, not divergence free
            val = -(2.d0*cos(freq*x)*M_PI*sin(freq*y)+2.d0*cos(freq*x)*cos(freq*y)*M_PI)
         end if

      elseif ( (abs(prob_sol)==20) .and. (abs(prob_coeff)==1) ) then 
         ! random rhs, Driven cavity BC

         call UniformRNG(val)

      elseif ( abs(prob_sol)==100 ) then 

         call UniformRNG(val)

      else 

         call bl_error('ERROR: bad or unimplemented choice for prob_sol/prob_coeff/theta_fac in rhs_val_2d')

      end if

    end function rhs_val_2d

    function rhs_val_3d(comp,x,y,z) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y,z
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: visc_coef, L, freq
      real(kind=dp_t) :: ufac,pfac,hx,hy,hz

      visc_coef=coeff_mag(2)/coeff_mag(1)

      if( (abs(prob_sol)==1) .and. (abs(prob_coeff)==1) .and. (theta_fac==0.0d0) ) then 
         ! Constant-coefficient Taylor vortex ! quasi-2D Taylor vortex

         L = prob_hi(1)-prob_lo(1)
         freq  = 2.d0*M_PI

         select case(abs(visc_type))      
         case(1,2,3)  ! as Taylor vortex is divergence free, visc_type=1, or 2, or 3, no difference
            if (comp .eq. 1) then
               ! staggered x-component
               ! for - div (T(u, p))
               val=-4.d0*M_PI*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*(-4.d0*M_PI*visc_coef*cos(freq*(-x+time)/L)*&
                    sin(freq*(-y+time)/L)+exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(4.d0*M_PI*(-x+time)/L)*L)/L**2   
               ! if only viscous term is concerned, uncomment the following code 
               ! val=-16.d0*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*cos(freq*(x-time)/L)*M_PI**2*&
               !               sin(freq*(y-time)/L)/L**2      
            else if (comp .eq. 2) then
               ! staggered y-component
               ! for - div(\nu D(u,p))
               val=-4.d0*M_PI*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*(4.d0*M_PI*visc_coef*sin(freq*(-x+time)/L)*&
                    cos(freq*(-y+time)/L)+exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(4.d0*M_PI*(-y+time)/L)*L)/L**2
               ! if only viscous term is concerned, uncomment the following code 
               ! val = 16.d0*exp(-8.d0*M_PI**2*visc_coef*time/L**2)*sin(freq*(x-time)/L)*M_PI**2*&
               !                cos(freq*(y-time)/L)/L**2
            else if (comp .eq. 3) then
               ! staggered z-component
               val = 0.d0
            else if (comp .eq. 4) then
               ! b_p
               val = 0.d0
            end if
         case default
            call bl_error('ERROR:rhs_val_3d, bad or unimplemented choice for visc_type')
         end select

      elseif( (abs(prob_sol)==2) .and. (abs(prob_coeff)==1) .and. (theta_fac==0.0d0) ) then

         freq  = 2.d0*M_PI
         ufac = dexp(-freq*freq*visc_coef*time)
         pfac = dexp(-2.0d0*freq*freq*visc_coef*time)
         hx = ABC_coefs(1)
         hy = ABC_coefs(2)
         hz = ABC_coefs(3)

         select case(abs(visc_type))      
         case(1,2,3)  ! as ABC flow is divergence free, visc_type=1, or 2, or 3, no difference
            if (comp .eq. 1) then
               ! staggered x-component
               val=4.d0*ufac*hz*cos(freq*(y-time))*M_PI**2+4.d0*ufac*hx*sin(freq*(z-time))*M_PI**2-ufac**2*&
                    (2.d0*hx*hy*cos(freq*(x-time))*M_PI*cos(freq*(z-time))-2.d0*hy*hz*sin(freq*(x-time))*M_PI*sin(freq*(y-time)))
            else if (comp .eq. 2) then
               ! staggered y-component
               val=4.d0*ufac*hy*sin(freq*(x-time))*M_PI**2+4.d0*ufac*hx*cos(freq*(z-time))*M_PI**2-ufac**2*&
                    (-2.d0*hx*hz*sin(freq*(y-time))*M_PI*sin(freq*(z-time))+2.d0*hy*hz*cos(freq*(x-time))*M_PI*cos(freq*(y-time)))
            else if (comp .eq. 3) then
               ! staggered z-component
               val = 4.d0*ufac*hy*cos(freq*(x-time))*M_PI**2+4.d0*ufac*hz*sin(freq*(y-time))*M_PI**2-ufac**2*&
                    (2.d0*hx*hz*cos(freq*(y-time))*M_PI*cos(freq*(z-time))-2.d0*hx*hy*sin(freq*(x-time))*M_PI*sin(freq*(z-time)))
            else if (comp .eq. 4) then
               ! b_p
               val = 0.d0
            end if
         case default
            call bl_error('ERROR: rhs_val_3d, bad or unimplemented choice for visc_type')
         end select

      elseif( (abs(prob_sol)==4) .or. (abs(prob_sol)==6) .or. (abs(prob_sol)==7) ) then   
         ! case=6 incompressible mode,7 stripe density

         val = 0.d0

      elseif( (abs(prob_sol)==3) .and. (abs(prob_coeff)==5) ) then 
         ! density=1, smooth viscosity sol     

         freq = 2.d0*M_PI

         if (comp .eq. 1) then
            ! staggered x-component
            val = -4.0d0*M_PI*exp(-8.0d0*M_PI**2*time)*(-2.0d0*sin(freq*x)*M_PI*cos(freq*x)+&
                 2.0d0*sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2+exp(-8.0d0*M_PI**2*time)*sin(4.0d0*M_PI*x))
         else if (comp .eq. 2) then
            ! staggered y-component
            val = -8.0d0*sin(freq*x)*M_PI**2*sin(freq*y)*exp(-8.0d0*M_PI**2*time)*&
                 (cos(freq*x)*sin(freq*y)+1.0d0)
         else
            ! staggered z-component OR pressure
            val = 0.d0
         end if

      elseif ( (abs(prob_sol)==31) .and. ((prob_coeff ==5) .or. (prob_coeff ==1)) )then 
         ! smoothly variable viscosity or vis=1, divergence not free solution, for checking discretization accuracy

         ! pre-assume that prob_dir=1
         freq  = 2.d0*M_PI

         if (comp .eq. 1) then
            ! staggered x-component

            if (abs(visc_type) .eq. 1) then 
               if (prob_coeff .eq. 5) then 
                  val = freq*(4.d0*sin(freq*x)*M_PI*cos(freq*x)-5.d0*sin(freq*x)*M_PI*&
                       cos(freq*x)*cos(freq*z)**2-5.d0*sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2+6.d0*&
                       sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2*cos(freq*z)**2+6.d0*sin(freq*x)* &
                       sin(freq*y)*sin(freq*z)*M_PI-2.d0*sin(4.d0*M_PI*x))
               else ! by default prob_coeff=1
                  val = 12.d0*sin(freq*x)*M_PI**2*sin(freq*y)*sin(freq*z) - 4.d0*sin(4.d0*M_PI*x)*M_PI
               end if
            else if (abs(visc_type) .eq. 2) then
               if (prob_coeff .eq. 5) then 
                  val = freq*(6.d0*sin(freq*x)*M_PI*cos(freq*x)-7.d0*sin(freq*x)*M_PI*&
                       cos(freq*x)*cos(freq*z)**2-7.d0*sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2+8.d0*&
                       sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2*cos(freq*z)**2+8.d0*sin(freq*x)*sin(freq*y)&
                       *sin(freq*z)*M_PI+2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)-2.d0*&
                       cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)*cos(freq*z)**2+2.d0*sin(freq*x)*&
                       cos(freq*y)*M_PI*sin(freq*z)-freq*cos(freq*x)**2*cos(freq*z)**2+freq*&
                       cos(freq*x)**2*(cos(freq*y)**2)*cos(freq*z)**2+2.d0*cos(freq*x)*M_PI*sin(freq*y)*sin(freq*z)&
                       +M_PI*cos(freq*x)**2-M_PI*(cos(freq*x)**2)*(cos(freq*y)**2)-2.d0*sin(4.d0*M_PI*x))
               else ! by default prob_coeff=1
                  val = 4.d0*M_PI*(4.d0*sin(freq*x)*sin(freq*y)*sin(freq*z)*M_PI+&
                       sin(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+cos(freq*x)*M_PI*sin(freq*y)&
                       *sin(freq*z)-sin(4.d0*M_PI*x))                              
               end if
            else if (abs(visc_type) .eq. 3) then
               if (prob_coeff .eq. 5) then 
                  val = freq*(6.d0*sin(freq*x)*M_PI*cos(freq*x)-7.d0*sin(freq*x)*M_PI*&
                       cos(freq*x)*cos(freq*z)**2-7.d0*sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2+8.d0*&
                       sin(freq*x)*M_PI*cos(freq*x)*cos(freq*y)**2*cos(freq*z)**2+8.d0*sin(freq*x)*sin(freq*y)&
                       *sin(freq*z)*M_PI+2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)-2.d0*&
                       cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)*cos(freq*z)**2+2.d0*sin(freq*x)*&
                       cos(freq*y)*M_PI*sin(freq*z)-freq*cos(freq*x)**2*cos(freq*z)**2+freq*&
                       cos(freq*x)**2*(cos(freq*y)**2)*cos(freq*z)**2+2.d0*cos(freq*x)*M_PI*sin(freq*y)*sin(freq*z)&
                       +M_PI*cos(freq*x)**2-M_PI*cos(freq*x)**2*(cos(freq*y)**2)-2.d0*sin(4.d0*M_PI*x))  &  ! the following line is diff(vis*div, x)
                       -(4.d0/3.d0)*sin(freq*z)*(8.d0*cos(freq*x)**2*sin(freq*z)-8.d0*cos(freq*x)**2*sin(freq*z)*&
                       cos(freq*y)**2+6.d0*cos(freq*x)**2*sin(freq*y)*sin(freq*z)*cos(freq*y)-4.d0*cos(freq*x)&
                       *sin(freq*z)*sin(freq*x)+4.d0*cos(freq*x)*sin(freq*z)*sin(freq*x)*cos(freq*y)**2+2.d0*&
                       sin(freq*x)*sin(freq*y)*sin(freq*z)*cos(freq*x)*cos(freq*y)-4.d0*sin(freq*z)+4.d0*&
                       sin(freq*z)*cos(freq*y)**2-3.d0*sin(freq*y)*sin(freq*z)*cos(freq*y)+2.d0*sin(freq*x)*&
                       sin(freq*y)+2.d0*sin(freq*x)*cos(freq*y)+2.d0*cos(freq*x)*sin(freq*y))*M_PI**2
               else ! by default prob_coeff=1
                  val = 4.d0*M_PI*(4.d0*sin(freq*x)*sin(freq*y)*sin(freq*z)*M_PI+&
                       sin(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+cos(freq*x)*M_PI*sin(freq*y) &
                       *sin(freq*z)-sin(4.d0*M_PI*x))    & ! the following is from div term
                       +(4.d0/3.d0)*M_PI**2*sin(freq*z)*(sin(freq*x)*sin(freq*y)+&
                       sin(freq*x)*cos(freq*y)+cos(freq*x)*sin(freq*y))
               end if
            else
               call bl_error('ERROR: init_rhs_3d, we only support visc_type <=3')
            end if

         else if (comp .eq. 2) then

            ! staggered y-component
            if (abs(visc_type) .eq. 1) then 
               if (prob_coeff .eq. 5) then 
                  val = -2.d0*sin(freq*x)**2*sin(freq*y)**2*sin(freq*z)**2*M_PI**2+(12.d0*(1.d0+0.5d0*&
                       cos(freq*x)*sin(freq*y)*sin(freq*z)))*cos(freq*x)*M_PI**2*sin(freq*y)*sin(freq*z)-2.d0*&
                       cos(freq*x)**2*cos(freq*y)**2*M_PI**2*sin(freq*z)**2-2.d0*cos(freq*x)**2*sin(freq*y)**2*&
                       cos(freq*z)**2*M_PI**2-4.d0*sin(4.d0*M_PI*y)*M_PI
               else ! by default prob_coeff=1
                  val = 12.d0*cos(freq*x)*M_PI**2*sin(freq*y)*sin(freq*z)-4.d0*sin(4.d0*M_PI*y)*M_PI
               end if
            else if (abs(visc_type) .eq. 2) then
               if (prob_coeff .eq. 5) then 
                  val = -freq*(2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)*&
                       cos(freq*z)**2-cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)+6.d0*M_PI*&
                       cos(freq*x)**2*cos(freq*z)**2-8.d0*M_PI*cos(freq*x)**2*cos(freq*y)**2*cos(freq*z)**2-5.d0*M_PI&
                       *cos(freq*x)**2+M_PI+2.d0*sin(4.d0*M_PI*y)-M_PI*cos(freq*z)**2-M_PI*cos(freq*y)**2-2.d0*sin(freq*x)*&
                       cos(freq*y)*M_PI*sin(freq*z)-8.d0*cos(freq*x)*M_PI*sin(freq*y)*sin(freq*z)-sin(freq*y)*&
                       M_PI*cos(freq*y)-2.d0*sin(freq*y)*M_PI*cos(freq*y)*cos(freq*x)**2*cos(freq*z)**2+M_PI*&
                       cos(freq*y)**2*cos(freq*z)**2+2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+sin(freq*y)*&
                       M_PI*cos(freq*y)*cos(freq*z)**2+2.d0*sin(freq*y)*M_PI*cos(freq*y)*cos(freq*x)**2&
                       +7.d0*M_PI*cos(freq*x)**2*cos(freq*y)**2)
               else ! by default prob_coeff=1
                  val = -4.d0*M_PI*(cos(freq*x)*cos(freq*y)*M_PI*sin(freq*z)-4.d0*cos(freq*x)*M_PI*&
                       sin(freq*y)*sin(freq*z)-sin(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+sin(4.d0*M_PI*y))
               end if
            else if (abs(visc_type) .eq. 3) then
               if (prob_coeff .eq. 5) then 
                  val = -freq*(2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)*&
                       cos(freq*z)**2-cos(freq*x)*cos(freq*y)*M_PI*sin(freq*x)*sin(freq*y)+6.d0*M_PI*&
                       cos(freq*x)**2*cos(freq*z)**2-8.d0*M_PI*cos(freq*x)**2*cos(freq*y)**2*cos(freq*z)**2-5.d0*M_PI&
                       *cos(freq*x)**2+M_PI+2.d0*sin(4.d0*M_PI*y)-M_PI*cos(freq*z)**2-M_PI*cos(freq*y)**2-2.d0*sin(freq*x)*&
                       cos(freq*y)*M_PI*sin(freq*z)-8.d0*cos(freq*x)*M_PI*sin(freq*y)*sin(freq*z)-sin(freq*y)*&
                       M_PI*cos(freq*y)-2.d0*sin(freq*y)*M_PI*cos(freq*y)*cos(freq*x)**2*cos(freq*z)**2+M_PI*&
                       cos(freq*y)**2*cos(freq*z)**2+2.d0*cos(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+sin(freq*y)*&
                       M_PI*cos(freq*y)*cos(freq*z)**2+2.d0*sin(freq*y)*M_PI*cos(freq*y)*cos(freq*x)**2&
                       +7.d0*M_PI*cos(freq*x)**2*cos(freq*y)**2)   &  ! the following line is diff(vis*div, y)
                       -(4.d0/3.d0)*sin(freq*z)*(8.d0*sin(freq*x)*sin(freq*y)*sin(freq*z)*cos(freq*x)*cos(freq*y)&
                       +6.d0*cos(freq*x)*sin(freq*z)*sin(freq*x)*cos(freq*y)**2-6.d0*sin(freq*y)*sin(freq*z)*&
                       cos(freq*y)+4.d0*cos(freq*x)**2*sin(freq*y)*sin(freq*z)*cos(freq*y)-2.d0*cos(freq*x)**2&
                       *sin(freq*z)*cos(freq*y)**2-3.d0*cos(freq*x)*sin(freq*z)*sin(freq*x)-2.d0*cos(freq*x)&
                       *cos(freq*y)+2.d0*cos(freq*x)*sin(freq*y)+2.d0*sin(freq*x)*cos(freq*y)&
                       +cos(freq*x)**2*sin(freq*z))*M_PI**2
               else ! by default prob_coeff=1
                  val = -4.d0*M_PI*(cos(freq*x)*cos(freq*y)*M_PI*sin(freq*z)-4.d0*cos(freq*x)*M_PI*&
                       sin(freq*y)*sin(freq*z)-sin(freq*x)*cos(freq*y)*M_PI*sin(freq*z)+sin(4.d0*M_PI*y)) & ! the following line is diff(vis*div, y)
                       +(4.d0/3.d0)*M_PI**2*sin(freq*z)*(-cos(freq*x)*cos(freq*y)+cos(freq*x)*&
                       sin(freq*y)+sin(freq*x)*cos(freq*y))                       
               end if
            else
               call bl_error('ERROR: init_rhs_3d, we only support visc_type <=3')
            end if

         else if (comp .eq. 3) then

            ! staggered z-component
            if (abs(visc_type) .eq. 1) then 
               if (prob_coeff .eq. 5) then 
                  val = freq*(5*sin(freq*x)*sin(freq*z)*M_PI*cos(freq*x)*cos(freq*z)-6.d0*&
                       sin(freq*x)*sin(freq*z)*M_PI*cos(freq*x)*cos(freq*z)*cos(freq*y)**2+6.d0*sin(freq*x)*&
                       sin(freq*y)*cos(freq*z)*M_PI-2.d0*sin(4.d0*M_PI*z))
               else ! by default prob_coeff=1
                  val = 12.d0*sin(freq*x)*M_PI**2*sin(freq*y)*cos(freq*z)-4.d0*sin(4.d0*M_PI*z)*M_PI
               end if
            else if (abs(visc_type) .eq. 2) then
               if (prob_coeff .eq. 5) then 
                  val = freq*(sin(freq*z)*M_PI*cos(freq*z)-sin(freq*z)*M_PI*cos(freq*z)*&
                       cos(freq*y)**2-2.d0*sin(freq*z)*M_PI*cos(freq*z)*cos(freq*x)**2+2.d0*sin(freq*z)*M_PI*&
                       cos(freq*z)*cos(freq*x)**2*cos(freq*y)**2+7.d0*sin(freq*x)*sin(freq*z)*M_PI*cos(freq*x)*&
                       cos(freq*z)-8.d0*cos(freq*x)*cos(freq*y)**2*M_PI*sin(freq*z)*sin(freq*x)*cos(freq*z)-2.d0*&
                       cos(freq*x)*M_PI*sin(freq*y)*cos(freq*z)+8.d0*sin(freq*x)*sin(freq*y)*cos(freq*z)*M_PI-& 
                       2.d0*cos(freq*x)**2*cos(freq*y)*M_PI*sin(freq*z)*sin(freq*y)*cos(freq*z)-2.d0*cos(freq*x)*&
                       M_PI*cos(freq*y)*cos(freq*z)-2.d0*sin(4.d0*M_PI*z))
               else ! by default prob_coeff=1
                  val = -4.d0*M_PI*(-4.d0*sin(freq*x)*sin(freq*y)*cos(freq*z)*M_PI+cos(freq*x)*&
                       M_PI*sin(freq*y)*cos(freq*z)+cos(freq*x)*M_PI*cos(freq*y)*cos(freq*z)+sin(4.d0*M_PI*z))
               end if
            else if (abs(visc_type) .eq. 3) then
               if (prob_coeff .eq. 5) then 
                  val = freq*(sin(freq*z)*M_PI*cos(freq*z)-sin(freq*z)*M_PI*cos(freq*z)*&
                       cos(freq*y)**2-2.d0*sin(freq*z)*M_PI*cos(freq*z)*cos(freq*x)**2+2.d0*sin(freq*z)*M_PI*&
                       cos(freq*z)*cos(freq*x)**2*cos(freq*y)**2+7.d0*sin(freq*x)*sin(freq*z)*M_PI*cos(freq*x)*&
                       cos(freq*z)-8.d0*cos(freq*x)*cos(freq*y)**2*M_PI*sin(freq*z)*sin(freq*x)*cos(freq*z)-2.d0*&
                       cos(freq*x)*M_PI*sin(freq*y)*cos(freq*z)+8.d0*sin(freq*x)*sin(freq*y)*cos(freq*z)*M_PI-& 
                       2.d0*cos(freq*x)**2*cos(freq*y)*M_PI*sin(freq*z)*sin(freq*y)*cos(freq*z)-2.d0*cos(freq*x)*&
                       M_PI*cos(freq*y)*cos(freq*z)-2.d0*sin(4.d0*M_PI*z)) &  ! the following line is diff(vis*div, z)
                       +(8.d0/3.d0)*cos(freq*z)*(-4.d0*cos(freq*x)*sin(freq*z)*sin(freq*x)+4.d0*cos(freq*x)*&
                       sin(freq*z)*sin(freq*x)*cos(freq*y)**2-3.d0*sin(freq*x)*sin(freq*y)*sin(freq*z)*&
                       cos(freq*x)*cos(freq*y)+3.d0*sin(freq*z)-3.d0*sin(freq*z)*cos(freq*y)**2-2.d0*&
                       cos(freq*x)**2*sin(freq*z)+2.d0*cos(freq*x)**2*sin(freq*z)*cos(freq*y)**2+cos(freq*x)**2&
                       *sin(freq*y)*sin(freq*z)*cos(freq*y)+cos(freq*x)*sin(freq*y)+cos(freq*x)*cos(freq*y)&
                       -sin(freq*x)*sin(freq*y))*M_PI**2
               else ! by default prob_coeff=1
                  val = -4.d0*M_PI*(-4.d0*sin(freq*x)*sin(freq*y)*cos(freq*z)*M_PI+cos(freq*x)*&
                       M_PI*sin(freq*y)*cos(freq*z)+cos(freq*x)*M_PI*cos(freq*y)*cos(freq*z)+sin(4.d0*M_PI*z)) &! the following line is diff(vis*div, z)
                       +(4.d0/3.d0)*M_PI**2*cos(freq*z)*(-cos(freq*x)*sin(freq*y)-cos(freq*x)*cos(freq*y)&
                       +sin(freq*x)*sin(freq*y))
               end if
            else
               call bl_error('ERROR: init_rhs_3d, we only support visc_type <=3')
            end if

         else if (comp .eq. 4) then

            ! -Div term, not divergence free
            if ((abs(visc_type) .lt. 1) .or. (abs(visc_type) .gt. 3)) then
               call bl_error('ERROR: init_rhs_3d, un-supportted visc_type ')
            else 
               val = -freq*sin(freq*z)*(cos(freq*x)*sin(freq*y)+&
                    cos(freq*x)*cos(freq*y)-sin(freq*x)*sin(freq*y))
            end if

         end if

      elseif ( (abs(prob_sol)==32) .and. (prob_coeff ==8) )then 
         ! smoothly variable density, for checking discretization accuracy of smooth density elliptic operator 

         freq  = 2.d0*M_PI
         select case(abs(visc_type))  
         case(1,2,3)
            if (comp .eq. 1) then
               ! staggered x-component
               val = (sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)/fixed_dt)*&
                    (1.0d0+0.5d0*cos(freq*x)*sin(freq*y)*sin(freq*z)) - 4.d0*M_PI*sin(4.d0*M_PI*x)
            elseif (comp .eq. 2) then
               ! staggered y-component  
               val = (cos(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*sin(2.d0*M_PI*z)/fixed_dt)*&
                    (1.0d0+0.5d0*cos(freq*x)*sin(freq*y)*sin(freq*z)) - 4.d0*M_PI*sin(4.d0*M_PI*y)
            elseif (comp .eq. 3) then 
               ! staggered z-component  
               val = (sin(2.d0*M_PI*x)*sin(2.d0*M_PI*y)*cos(2.d0*M_PI*z)/fixed_dt)*&
                    (1.0d0+0.5d0*cos(freq*x)*sin(freq*y)*sin(freq*z)) - 4.d0*M_PI*sin(4.d0*M_PI*z)    
            end if

         case default
            call bl_error('ERROR:init_rhs_2d, bad or unimplemented choice for prob_sol for testing accuracy')
         end select

         if (comp .eq. 4) then
            ! -Div term, not divergence free
            if ((abs(visc_type) .lt. 1) .or. (abs(visc_type) .gt. 3)) then
               call bl_error('ERROR: init_rhs_3d, un-supportted visc_type ')
            else 
               val = -freq*sin(freq*z)*(cos(freq*x)*sin(freq*y)+&
                    cos(freq*x)*cos(freq*y)-sin(freq*x)*sin(freq*y))
            end if
         end if

      elseif ( (abs(prob_sol)==20) .and. (abs(prob_coeff)==1) ) then 
         ! random rhs, Driven cavity BC

         call UniformRNG(val)

      elseif ( abs(prob_sol)==100 ) then 

         call UniformRNG(val)

      else 

         call bl_error('ERROR: bad or unimplemented choice for prob_sol/prob_coeff/theta_fac in rhs_val_3d')

      end if

    end function rhs_val_3d

  end subroutine init_rhs

  subroutine init_mat(mla,alpha,beta,gamma,dx,time,the_bc_level)

    use probin_module       , only: prob_coeff, prob_sol, smoothing_width, coeff_ratio, &
                                    var_coeff_mag, coeff_mag
    use probin_common_module, only: n_cells, max_grid_size, visc_type, prob_lo, prob_hi

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: alpha(:)
    type(multifab) , intent(inout) :: beta(:)
    type(multifab) , intent(inout) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: ap(:,:,:,:), bp(:,:,:,:), cp(:,:,:,:), gp(:,:,:,:)

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: n,nlevs,i,dm,ng_a,ng_b,ng_g

    type(box)      :: bx
    type(boxarray) :: ba_onebox,ba_manybox
    type(layout)   :: la_onebox,la_manybox
    type(multifab) :: coeff_onebox(mla%nlevel),coeff_manybox(mla%nlevel)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_a = alpha(1)%ng
    ng_b = beta(1)%ng
    ng_g = gamma(1)%ng

    ! check the consistency of visc_type: only when prob_coeff=1 or 3, we use constant viscosity
    if ((abs(var_coeff_mag(1)) .gt. 0.d0) .or. (abs(var_coeff_mag(2)) .gt. 0.d0) .or. (abs(var_coeff_mag(3)) .gt. 0.d0)) then 
       if (visc_type >0) then
          call bl_error('ERROR: make sure the sign of visc_type in input files is correct')
       end if
    end if

    ! Deterministic part 
    do n=1,nlevs
       do i = 1, nfabs(alpha(n))
          ap  => dataptr(alpha(n),i)
          bp  => dataptr(beta (n),i)
          gp  => dataptr(gamma(n),i)
          lo =  lwb(get_box(alpha(n),i))
          hi =  upb(get_box(alpha(n),i))
          select case (dm)
          case (2)
             call init_mat_2d(ap(:,:,1,1), ng_a, &
                              bp(:,:,1,1), ng_b, &
                              gp(:,:,1,1), ng_g, &
                              lo, hi, dx(n,:), time)
          case (3)
             call init_mat_3d(ap(:,:,:,1), ng_a, &
                              bp(:,:,:,1), ng_b, &
                              gp(:,:,:,1), ng_g, &
                              lo, hi, dx(n,:), time)
          end select
       end do
    end do

    ! Random part 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! special code to fill random coefficients
    !
    ! make a new multifab with the same grids, but with grown-by-one boxes and 0 ghost cells
    ! then make another new multifab with only 1 box, covers the entire domain grown-by-one,
    !   and 0 ghost cells.  fill this with random.
    ! then copy random from one big box to special many-box multifab
    ! then copy into the original multifab
    ! we have to do this since parallel copies only work from valid-to-valid region

    if (nlevs .gt. 1) then
       call bl_error("init_mat for random coefficients does not work for multilevel")
    end if

    ! create a box from (0,0) to (n_cells-1,n_cells-1)
    lo(1:dm) = 0
    hi(1:dm) = n_cells(1:dm)-1
    bx = make_box(lo,hi)

    ! build a boxarray with a single box
    call boxarray_build_bx(ba_manybox,bx)

    ! chop up the box to respect max_grid_size
    call boxarray_maxsize(ba_manybox,max_grid_size)

    ! grow each box by 1
    call boxarray_grow_n(ba_manybox,1)

    ! build a layout with the same proc mapping as the standard multifabs build with mla
    call layout_build_ba(la_manybox,ba_manybox,bx,mla%pmask,explicit_mapping=get_proc(mla%la(1)))

    ! don't need this anymore - free up memory
    call destroy(ba_manybox)

    ! build the a many-boxed multifab, but with grown boxes and no ghost cells
    call multifab_build(coeff_manybox(1),la_manybox,3,0)

    ! build a boxarray with a single box
    call boxarray_build_bx(ba_onebox,bx)

    ! grow the box by 1
    call boxarray_grow_n(ba_onebox,1)

    ! build a layout for a multifab with 1 box with valid region grown by 1
    call layout_build_ba(la_onebox,ba_onebox,bx,mla%pmask)

    ! don't need this anymore - free up memory
    call destroy(ba_onebox)

    ! build the special multifab with 1 box with valid region grown by 1
    call multifab_build(coeff_onebox(1),la_onebox,3,0)

    ! fill coeff_onebox with random numbers
    do n=1,nlevs
       do i=1,nfabs(coeff_onebox(n))
          cp => dataptr(coeff_onebox(n),i)
          select case (dm)
          case (2)
             call UniformRNGs(cp(:,:,1,:), size(cp(:,:,1,:)))
             cp(:,:,1,1) = cp(:,:,1,1)*var_coeff_mag(1)
             cp(:,:,1,2) = cp(:,:,1,2)*var_coeff_mag(2)
             cp(:,:,1,3) = cp(:,:,1,3)*var_coeff_mag(3)
          case (3)
             call UniformRNGs(cp(:,:,:,:), size(cp(:,:,:,:)))
             cp(:,:,:,1) = cp(:,:,:,1)*var_coeff_mag(1)
             cp(:,:,:,2) = cp(:,:,:,2)*var_coeff_mag(2)
             cp(:,:,:,3) = cp(:,:,:,3)*var_coeff_mag(3)
          end select
       end do
    end do

    ! copy data from one box into many boxes
    do n=1,nlevs
       call multifab_copy_c(coeff_manybox(n),1,coeff_onebox(n),1,3,0)
    end do

    ! copy data from many boxes into original multifab
    do n=1,nlevs
       do i=1,nfabs(coeff_manybox(n))
          cp => dataptr(coeff_manybox(n),i)
          ap => dataptr(alpha(n),i)
          bp => dataptr(beta (n),i)
          gp => dataptr(gamma(n),i)
          select case (dm)
          case (2)
             ap(:,:,1,1) = ap(:,:,1,1)*(1.d0+cp(:,:,1,1))*coeff_mag(1)
             bp(:,:,1,1) = bp(:,:,1,1)*(1.d0+cp(:,:,1,2))*coeff_mag(2)
             gp(:,:,1,1) = gp(:,:,1,1)*(1.d0+cp(:,:,1,3))*coeff_mag(3)
          case (3)
             ap(:,:,:,1) = ap(:,:,:,1)*(1.d0+cp(:,:,:,1))*coeff_mag(1)
             bp(:,:,:,1) = bp(:,:,:,1)*(1.d0+cp(:,:,:,2))*coeff_mag(2)
             gp(:,:,:,1) = gp(:,:,:,1)*(1.d0+cp(:,:,:,3))*coeff_mag(3)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(coeff_onebox(n))
       call multifab_destroy(coeff_manybox(n))
    end do
    call destroy(la_onebox)
    call destroy(la_manybox)

    ! end special code for random coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       call multifab_fill_boundary(alpha(n))
       call multifab_fill_boundary(beta(n))
       call multifab_fill_boundary(gamma(n))
       call multifab_coefbc(alpha(n),the_bc_level(n))
       call multifab_coefbc( beta(n),the_bc_level(n))
       call multifab_coefbc(gamma(n),the_bc_level(n))
    enddo

  contains

    subroutine init_mat_2d (alpha,ng_a,beta,ng_b,gamma,ng_g,lo,hi,dx,time)

      integer, intent(in) :: lo(:), hi(:), ng_a, ng_b, ng_g
      real(kind=dp_t), intent(inout) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(inout) ::  beta(lo(1)-ng_b:,lo(2)-ng_b:)
      real(kind=dp_t), intent(inout) :: gamma(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      ! local 
      real(kind=dp_t) :: x, y
      integer :: i, j

      do j=lo(2)-ng_a,hi(2)+ng_a
         y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
         do i=lo(1)-ng_a,hi(1)+ng_a
            x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
            alpha(i,j) = mat_val_2d(1,x,y,dx)
         end do
      end do

      do j=lo(2)-ng_b,hi(2)+ng_b
         y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
         do i=lo(1)-ng_b,hi(1)+ng_b
            x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
            beta(i,j) = mat_val_2d(2,x,y,dx)
         end do
      end do

      do j=lo(2)-ng_g,hi(2)+ng_g
         y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
         do i=lo(1)-ng_g,hi(1)+ng_g
            x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
            gamma(i,j) = mat_val_2d(3,x,y,dx)
         end do
      end do

    end subroutine init_mat_2d

    subroutine init_mat_3d (alpha,ng_a,beta,ng_b,gamma,ng_g,lo,hi,dx,time)

      integer, intent(in) :: lo(:), hi(:), ng_a, ng_b, ng_g
      real(kind=dp_t), intent(inout) :: alpha(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout) ::  beta(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
      real(kind=dp_t), intent(inout) :: gamma(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      ! local 
      real(kind=dp_t) :: x, y, z
      integer :: i, j, k

      do k=lo(3)-ng_a,hi(3)+ng_a
         z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
         do j=lo(2)-ng_a,hi(2)+ng_a
            y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
            do i=lo(1)-ng_a,hi(1)+ng_a
               x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
               alpha(i,j,k) = mat_val_3d(1,x,y,z,dx)
            end do
         end do
      end do

      do k=lo(3)-ng_b,hi(3)+ng_b
         z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
         do j=lo(2)-ng_b,hi(2)+ng_b
            y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
            do i=lo(1)-ng_b,hi(1)+ng_b
               x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
               beta(i,j,k) = mat_val_3d(2,x,y,z,dx)
            end do
         end do
      end do

      do k=lo(3)-ng_g,hi(3)+ng_g
         z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
         do j=lo(2)-ng_g,hi(2)+ng_g
            y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
            do i=lo(1)-ng_g,hi(1)+ng_g
               x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
               gamma(i,j,k) = mat_val_3d(3,x,y,z,dx)
            end do
         end do
      end do

    end subroutine init_mat_3d

    function mat_val_2d(comp,x,y,dx) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y,dx(:)
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: freq,r,one_third_domain1,one_third_domain2,y1,y2,interface_dist

      select case (abs(prob_coeff))
      case (1,2,3,4)

         ! set alpha, beta, and gamma to constant=coeff_mag
         val = 1.d0

      case (5)

         ! density = 1, smooth viscosity and gamma
         freq = 2.d0*M_PI

         if (comp .eq. 1) then
            ! alpha
            val = 1.d0
         else if (comp .eq. 2) then
            ! beta
            val = 1.0d0+0.5d0*cos(freq*x)*sin(freq*y)
         else if (comp .eq. 3) then
            ! gamma
            if (prob_sol .eq. 31) then
               val = sin(freq*x)*sin(freq*y)
            else
               val = 1.d0
            end if
         end if

      case (6)

         ! density and viscosity jump at a bubble interface, multiphase flow
         r = sqrt( (x-0.5*(prob_lo(1)+prob_hi(1)))**2&
              +(y-0.5*(prob_lo(2)+prob_hi(2)))**2)
         interface_dist=min(0.25d0*(prob_hi(1)-prob_lo(1)), &
              0.25d0*(prob_hi(2)-prob_lo(2)))

         select case (comp) 
         case (1)  
            ! alpha 
            val = 0.5d0*(coeff_ratio(1)+1.0d0)+0.5d0*(coeff_ratio(1)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case (2)
            ! beta
            val = 0.5d0*(coeff_ratio(2)+1.0d0)+0.5d0*(coeff_ratio(2)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case (3)
            ! gamma
            val = 0.5d0*(coeff_ratio(3)+1.0d0)+0.5d0*(coeff_ratio(3)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case default
            call bl_error('ERROR in mat_val_2d: comp should not be greater than 3')
         end select

      case (7)

         ! mcai: bilayer interface (stripe), for both density function and viscosity
         one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
         one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)
         y1 = y - one_third_domain1
         y2 = y - one_third_domain2

         select case (comp) 
         case (1)  
            ! alpha  
            val = 0.5d0*(coeff_ratio(1)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case (2) 
            ! beta  
            val = 0.5d0*(coeff_ratio(2)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case (3) 
            ! gamma
            val = 0.5d0*(coeff_ratio(3)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case default
            call bl_error('ERROR in mat_val_3d: comp should not be greater than 3')
         end select

      case (8) ! for checking variable density discretization accuracy
         !  smooth density 
         freq = 2.d0*M_PI
         if (comp .eq. 1) then
            ! alpha
            val = 1.0d0+0.5d0*cos(freq*x)*sin(freq*y)
         else if (comp .eq. 2) then
            ! beta
            val = 0.0d0
         else if (comp .eq. 3) then
            ! gamma
            val = 0.d0            
         end if

      case default
         call bl_error('ERROR in mat_val_3d: bad or unimplemented choice for density/vis in 2D')
      end select

    end function mat_val_2d

    function mat_val_3d(comp,x,y,z,dx) result(val)

      use bl_constants_module

      integer        , intent(in   ) :: comp
      real(kind=dp_t), intent(in   ) :: x,y,z,dx(:)
      real(kind=dp_t)                :: val

      ! local
      real(kind=dp_t) :: freq,r,one_third_domain1,one_third_domain2,y1,y2,interface_dist

      select case (abs(prob_coeff))
      case (1,2,3,4)

         ! set alpha, beta, and gamma to 1.0
         val = 1.d0
         !val = 0.d0   ! testing random coefficient case, 1/rho jump coefficient

         ! Note: for case=2, we will later overwrite beta with random
         !       for case=3, we will later overwrite alpha with random
         !       for case=4, we will later overwrite, alpha/beta/gamma with random

      case (5)

         ! density = 1, smooth viscosity and gamma
         freq = 2.d0*M_PI

         if (comp .eq. 1) then
            ! alpha
            val = 1.d0
         else if (comp .eq. 2) then
            ! beta
            val = 1.d0+0.5d0*cos(freq*x)*sin(freq*y)*sin(freq*z)
         else if (comp .eq. 3) then
            ! gamma
            if (prob_sol .eq. 31) then
               val = sin(freq*x)*sin(freq*y)*sin(freq*z)
            else
               val = 1.d0
            end if
         end if

      case (6)

         ! density and viscosity jump at a bubble interface, multiphase flow
         r = sqrt( (x-0.5*(prob_lo(1)+prob_hi(1)))**2 &
              +(y-0.5*(prob_lo(2)+prob_hi(2)))**2 &
              +(z-0.5*(prob_lo(3)+prob_hi(3)))**2)
         interface_dist=min(0.25d0*(prob_hi(1)-prob_lo(1)), &
              0.25d0*(prob_hi(2)-prob_lo(2)), &
              0.25d0*(prob_hi(2)-prob_lo(2)))

         select case (comp) 
         case (1)  
            ! alpha   
            val = 0.5d0*(coeff_ratio(1)+1.0d0)+0.5d0*(coeff_ratio(1)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case (2) 
            ! beta  
            val = 0.5d0*(coeff_ratio(2)+1.0d0)+0.5d0*(coeff_ratio(2)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case (3) 
            ! gamma
            val = 0.5d0*(coeff_ratio(3)+1.0d0)+0.5d0*(coeff_ratio(3)-1.0d0)*tanh((r-interface_dist)/(smoothing_width*dx(1)))
         case default
            call bl_error('ERROR in mat_val_3d: comp should not be greater than 3')
         end select

      case (7)

         ! mcai: bilayer interface (stripe), for both density function and viscosity
         one_third_domain1=2.0d0/3.0d0*prob_lo(2)+1.0d0/3.0d0*prob_hi(2)
         one_third_domain2=1.0d0/3.0d0*prob_lo(2)+2.0d0/3.0d0*prob_hi(2)
         y1 = y - one_third_domain1
         y2 = y - one_third_domain2

         select case (comp) 
         case (1)
            ! alpha          
            val = 0.5d0*(coeff_ratio(1)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case (2)
            ! beta  
            val = 0.5d0*(coeff_ratio(2)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case (3) 
            ! gamma
            val = 0.5d0*(coeff_ratio(3)-1.d0)*(tanh(y1/(smoothing_width*dx(2)))-tanh(y2/(smoothing_width*dx(2))))+1.d0
         case default
            call bl_error('ERROR in mat_val_2d: comp should not be greater than 3')
         end select

      case (8) ! for checking variable density discretization accuracy
         !  smooth density 
         freq = 2.d0*M_PI
         if (comp .eq. 1) then
            ! alpha
            val = 1.0d0+0.5d0*cos(freq*x)*sin(freq*y)*sin(freq*z)
         else if (comp .eq. 2) then
            ! beta
            val = 0.0d0
         else if (comp .eq. 3) then
            ! gamma
            val = 0.d0            
         end if

      case default
         call bl_error('ERROR in mat_val_3d: bad or unimplemented choice for density/vis in 2D')
      end select

    end function mat_val_3d

  end subroutine init_mat

  subroutine multifab_coefbc(s,the_bc_level)

    ! this fills ghost cells for transport coefficients at walls by
    ! averaging ghost and interior into ghost

    type(multifab) , intent(inout) :: s
    type(bc_level) , intent(in   ) :: the_bc_level
   
    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          call coefbc_2d(sp(:,:,1,1), lo, hi, ng, &
                         the_bc_level%phys_bc_level_array(i,:,:))
       case (3)
          call coefbc_3d(sp(:,:,:,1), lo, hi, ng, &
                         the_bc_level%phys_bc_level_array(i,:,:))
       end select
    end do
 
  end subroutine multifab_coefbc

  subroutine coefbc_2d(s,lo,hi,ng,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = 0.5d0*(s(lo(1),j)+s(lo(1)-1,j))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = 0.5d0*(s(hi(1),j)+s(hi(1)+1,j))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = 0.5d0*(s(i,lo(2))+s(i,lo(2)-1))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = 0.5d0*(s(i,hi(2))+s(i,hi(2)+1))
       end do
    end if

  end subroutine coefbc_2d

  subroutine coefbc_3d(s,lo,hi,ng,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j,k

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = 0.5d0*(s(lo(1),j,k)+s(lo(1)-1,j,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = 0.5d0*(s(hi(1),j,k)+s(hi(1)+1,j,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = 0.5d0*(s(i,lo(2),k)+s(i,lo(2)-1,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = 0.5d0*(s(i,hi(2),k)+s(i,hi(2)+1,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. SLIP_WALL .or. bc(3,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = 0.5d0*(s(i,j,lo(3))+s(i,j,lo(3)-1))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. SLIP_WALL .or. bc(3,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = 0.5d0*(s(i,j,hi(3))+s(i,j,hi(3)+1))
          end do
       end do
    end if

  end subroutine coefbc_3d

end module init_module
