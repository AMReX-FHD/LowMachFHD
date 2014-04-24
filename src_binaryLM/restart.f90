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
                               time,dt)


    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chkdata_edgex(:)
    type(multifab)   , pointer        :: chkdata_edgey(:)
    type(multifab)   , pointer        :: chkdata_edgez(:)
    character(len=9)                  :: sd_name
    integer                           :: n,nlevs,dm
    integer                           :: rrs(10)

    dm = dim_in

    write(unit=sd_name,fmt='("chk",i6.6)') restart_int

    if ( parallel_IOProcessor() ) then
       print *,'Reading ',sd_name,' to get state data for restart'
    end if

    call checkpoint_read(chkdata, chkdata_edgex, chkdata_edgey, chkdata_edgez, &
                         sd_name, rrs, time, dt, nlevs)

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
