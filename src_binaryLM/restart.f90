module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use probin_common_module, only: dim_in

  implicit none

  private

  public :: fill_restart_data

contains

  subroutine fill_restart_data(restart_int,mba,chkdata, &
                               chkdata_edgex,chkdata_edgey,chkdata_edgez, &
                               chkdata_vel_bc_t1,chkdata_vel_bc_t2,chkdata_vel_bc_t3, &
                               chkdata_vel_bc_t4,chkdata_vel_bc_t5,chkdata_vel_bc_t6, &
                               time,dt)


    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chkdata_edgex(:)
    type(multifab)   , pointer        :: chkdata_edgey(:)
    type(multifab)   , pointer        :: chkdata_edgez(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t1(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t2(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t3(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t4(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t5(:)
    type(multifab)   , pointer        :: chkdata_vel_bc_t6(:)
    character(len=9)                  :: sd_name
    integer                           :: n,nlevs,dm
    integer                           :: rrs(10)

    dm = dim_in

    write(unit=sd_name,fmt='("chk",i6.6)') restart_int

    if ( parallel_IOProcessor() ) then
       print *,'Reading ',sd_name,' to get state data for restart'
    end if

    if (dm .eq. 2) then
       call checkpoint_read_2d(chkdata, chkdata_edgex, chkdata_edgey, &
                               chkdata_vel_bc_t1, chkdata_vel_bc_t2, &
                               sd_name, rrs, time, dt, nlevs)
    else if (dm .eq. 3) then
       call checkpoint_read_3d(chkdata, chkdata_edgex, chkdata_edgey, chkdata_edgez, &
                               chkdata_vel_bc_t1, chkdata_vel_bc_t2, chkdata_vel_bc_t3, &
                               chkdata_vel_bc_t4, chkdata_vel_bc_t5, chkdata_vel_bc_t6, &
                               sd_name, rrs, time, dt, nlevs)

    end if

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = rrs(n-1)
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

  end subroutine fill_restart_data

end module restart_module
