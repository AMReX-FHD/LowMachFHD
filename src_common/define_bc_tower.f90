module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module

  implicit none

  type bc_level

     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: max_level_built = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()
     integer       , pointer :: domain_bc(:,:) => Null()

  end type bc_tower

  private

  public :: bc_level, bc_tower, &
            initialize_bc, bc_tower_init, bc_tower_level_build, bc_tower_destroy

contains

  subroutine initialize_bc(the_bc_tower,num_levs,dm,pmask)

     use bc_module

     use probin_common_module, only : bc_lo, bc_hi

     type(bc_tower), intent(  out) :: the_bc_tower
     integer       , intent(in   ) :: num_levs,dm
     logical       , intent(in   ) :: pmask(:)

     integer :: domain_phys_bc(dm,2)

     ! Define the physical boundary conditions on the domain
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bc_lo(1)
     domain_phys_bc(1,2) = bc_hi(1)
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bc_lo(1) .ne. -1 .or. bc_hi(1) .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bc_lo(2)
        domain_phys_bc(2,2) = bc_hi(2)
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bc_lo(2) .ne. -1 .or. bc_hi(2) .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bc_lo(3)
        domain_phys_bc(3,2) = bc_hi(3)
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bc_lo(3) .ne. -1 .or. bc_hi(3) .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.
     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

  end subroutine initialize_bc

  subroutine bc_tower_init(bct,num_levs,dm,phys_bc_in)

    type(bc_tower ), intent(  out) :: bct
    integer        , intent(in   ) :: num_levs
    integer        , intent(in   ) :: dm
    integer        , intent(in   ) :: phys_bc_in(:,:)

    allocate(bct%bc_tower_array(num_levs))
    allocate(bct%domain_bc(dm,2))

    bct%domain_bc(:,:) = phys_bc_in(:,:)

  end subroutine bc_tower_init

  subroutine bc_tower_level_build(bct,n,la)

    type(bc_tower ), intent(inout) :: bct
    integer        , intent(in   ) :: n
    type(layout)   , intent(in   ) :: la

    integer :: ngrids,dm
    integer :: default_value

    if (associated(bct%bc_tower_array(n)%phys_bc_level_array)) then
      deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
      deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
      deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
      bct%bc_tower_array(n)%phys_bc_level_array => NULL()
      bct%bc_tower_array(n)%adv_bc_level_array => NULL()
      bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end if

    ngrids = layout_nlocal(la)
    dm = layout_dim(la)

    allocate(bct%bc_tower_array(n)%phys_bc_level_array(0:ngrids,dm,2))
    default_value = INTERIOR
    call phys_bc_level_build(bct%bc_tower_array(n)%phys_bc_level_array,la, &
                             bct%domain_bc,default_value)

    ! Here we allocate dm components for x_u,
    !                  1 component for x_p,
    !                  1 component for coefficients (alpha, beta, gamma, rho, rho1)
    allocate(bct%bc_tower_array(n)%adv_bc_level_array(0:ngrids,dm,2,dm+2))
    default_value = INTERIOR
    call adv_bc_level_build(bct%bc_tower_array(n)%adv_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

    ! Here we allocate dm components for x_u,
    !                  1 component for x_p,
    !                  1 component for coefficients (alpha, beta, gamma, rho, rho1)
    allocate(bct%bc_tower_array(n)%ell_bc_level_array(0:ngrids,dm,2,dm+2))
    default_value = BC_INT
    call ell_bc_level_build(bct%bc_tower_array(n)%ell_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

     bct%max_level_built = n

  end subroutine bc_tower_level_build

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: n

    do n = 1,bct%max_level_built
       deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
       bct%bc_tower_array(n)%phys_bc_level_array => NULL()
       bct%bc_tower_array(n)%adv_bc_level_array => NULL()
       bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end do
    deallocate(bct%bc_tower_array)

    deallocate(bct%domain_bc)

  end subroutine bc_tower_destroy

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: d,i

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do d = 1,layout_dim(la_level)
       phys_bc_level(i,d,1) = domain_bc(d,1)
       phys_bc_level(i,d,2) = domain_bc(d,2)
    end do

    do i = 1,layout_nlocal(la_level)
       bx = layout_get_box(la_level,global_index(la_level,i))
       do d = 1,layout_dim(la_level)
          if (lwb(bx,d) == lwb(pd,d)) phys_bc_level(i,d,1) = domain_bc(d,1)
          if (upb(bx,d) == upb(pd,d)) phys_bc_level(i,d,2) = domain_bc(d,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm,igrid,d,lohi

    ! these boundary conditions are written for a stokes system where we 
    ! solve for the INCREMENT to velocity and pressure

    adv_bc_level = default_value

    dm = size(adv_bc_level,dim=2)

    do igrid = 0, size(adv_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2

       if (phys_bc_level(igrid,d,lohi) == SLIP) then

          ! for normal velocity we impose a Dirichlet velocity condition
          ! for transverse velocity we imposie a Dirichlet stress condition
          adv_bc_level(igrid,d,lohi,1:dm) = DIR_TRACT ! transverse velocity
          adv_bc_level(igrid,d,lohi,d   ) = DIR_VEL   ! normal velocity
          adv_bc_level(igrid,d,lohi,dm+1) = FOEXTRAP  ! pressure
          adv_bc_level(igrid,d,lohi,dm+2) = FOEXTRAP  ! cell-centered prims and coeffs

       else if (phys_bc_level(igrid,d,lohi) == NO_SLIP) then

          ! for normal velocity we impose a Dirichlet velocity condition
          ! for transverse velocity we imposie a Dirichlet velocity condition
          adv_bc_level(igrid,d,lohi,1:dm) = DIR_VEL    ! transverse velocity
          adv_bc_level(igrid,d,lohi,1:dm) = DIR_VEL    ! normal velocity
          adv_bc_level(igrid,d,lohi,dm+1) = FOEXTRAP   ! pressure
          adv_bc_level(igrid,d,lohi,dm+2) = FOEXTRAP   ! cell-centered prims and coeffs

       else if (phys_bc_level(igrid,d,lohi) == INLET) then

          ! for normal velocity we impose a Dirichlet velocity condition
          ! for transverse velocity we imposie a Dirichlet velocity condition
          adv_bc_level(igrid,d,lohi,1:dm) = DIR_VEL    ! transverse velocity
          adv_bc_level(igrid,d,lohi,1:dm) = DIR_VEL    ! normal velocity
          adv_bc_level(igrid,d,lohi,dm+1) = FOEXTRAP   ! pressure
          adv_bc_level(igrid,d,lohi,dm+2) = EXT_DIR    ! cell-centered prims and coeffs

       else if (phys_bc_level(igrid,d,lohi) == PERIODIC .or. &
                phys_bc_level(igrid,d,lohi) == INTERIOR ) then

          ! retain the default value of INTERIOR

       else

          print*,'adv_bc_level_build',igrid,d,lohi,phys_bc_level(igrid,d,lohi)
          call bl_error('BC TYPE NOT SUPPORTED')

       end if

    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm
    integer :: igrid,d,lohi
    integer :: press_comp

    ell_bc_level = default_value
 
    dm = size(ell_bc_level,dim=2)

    press_comp = dm+1

    do igrid = 0, size(ell_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2

       if (phys_bc_level(igrid,d,lohi) == SLIP) then

          ell_bc_level(igrid,d,lohi,press_comp) = BC_NEU ! pressure

       else if (phys_bc_level(igrid,d,lohi) == NO_SLIP) then

          ell_bc_level(igrid,d,lohi,press_comp) = BC_NEU ! pressure

       else if (phys_bc_level(igrid,d,lohi) == INLET) then

          ell_bc_level(igrid,d,lohi,press_comp) = BC_NEU ! pressure

       else if (phys_bc_level(igrid,d,lohi) == PERIODIC) then

          ell_bc_level(igrid,d,lohi,press_comp) = BC_PER ! pressure

       else if (phys_bc_level(igrid,d,lohi) == INTERIOR) then

          ! retain the default value of INTERIOR

       else

          print*,'ell_bc_level_build',igrid,d,lohi,phys_bc_level(igrid,d,lohi)
          call bl_error('BC TYPE NOT SUPPORTED')

       end if

    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
