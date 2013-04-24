module checkpoint_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: write_checkfile, checkpoint_write, checkpoint_read

contains

  subroutine write_checkfile(mla,s_in,m_in,time,dt,istep_to_write)

    use ml_layout_module
    use probin_module, only: nscal
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s_in(:)
    type(multifab) , intent(in   ) :: m_in(:,:)
    integer        , intent(in   ) :: istep_to_write
    real(kind=dp_t), intent(in   ) :: time,dt

    type(multifab), pointer ::  chkdata(:)
    type(multifab), pointer ::  chkdata_nodal(:,:)

    integer :: n,nlevs,dm,i

    character(len=9) :: sd_name

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(chkdata(nlevs))
    allocate(chkdata_nodal(nlevs,dm))
    do n = 1,nlevs
       call multifab_build(chkdata(n), mla%la(n), nscal, 0)
       call multifab_copy_c(chkdata(n), 1, s_in(n), 1, nscal)
       do i=1,dm
          call multifab_build_edge(chkdata_nodal(n,i), mla%la(n), 1, 0, i)
          call multifab_copy_c(chkdata_nodal(n,i), 1, m_in(n,i), 1, 1)
       end do
    end do
    write(unit=sd_name,fmt='("chk",i6.6)') istep_to_write

    call checkpoint_write(nlevs, sd_name, chkdata, chkdata_nodal, mla%mba%rr, time, dt)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
       do i=1,dm
          call multifab_destroy(chkdata_nodal(n,i))
       end do
    end do
    deallocate(chkdata)

  end subroutine write_checkfile

  subroutine checkpoint_write(nlevs_in, dirname, mfs, mfs_nodal, rrs, time_in, dt_in)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
    use parallel
    use probin_module, only: verbose

    integer         , intent(in) :: nlevs_in
    type(multifab)  , intent(in) :: mfs(:)
    type(multifab)  , intent(in) :: mfs_nodal(:,:)
    integer         , intent(in) :: rrs(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t) , intent(in) :: time_in, dt_in

    integer :: n
    character(len=128) :: header, sd_name
    integer :: un, dm

    integer         :: nlevs
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ time
    namelist /chkpoint/ dt
    namelist /chkpoint/ nlevs

    dm = get_dim(mfs(1))
    if ( parallel_IOProcessor() ) call fabio_mkdir(dirname)

    call parallel_barrier() ! All CPUs have to wait till the directory is built.

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      print *,'Writing cc state to checkpoint file ',trim(sd_name)
      print *,' '
    end if

    write(unit=sd_name, fmt='(a,"/State_nodalx")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal(:,1), rrs(:,1), sd_name)

    write(unit=sd_name, fmt='(a,"/State_nodaly")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal(:,2), rrs(:,1), sd_name)

    if (dm .eq. 3) then
       write(unit=sd_name, fmt='(a,"/State_nodalz")') trim(dirname)
       call fabio_ml_multifab_write_d(mfs_nodal(:,3), rrs(:,1), sd_name)
    end if

    time  = time_in
    dt = dt_in
    nlevs = nlevs_in

    if (parallel_IOProcessor()) then
       header = "Header"
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un, nml = chkpoint)
       do n = 1,nlevs-1
          write(unit=un,fmt=*) rrs(n,1)
       end do
       close(un)
    end if

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodalx, mfs_nodaly, mfs_nodalz, dirname, &
                             rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    type(multifab)  ,                pointer :: mfs_nodalx(:)
    type(multifab)  ,                pointer :: mfs_nodaly(:)
    type(multifab)  ,                pointer :: mfs_nodalz(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, dm, nlevs
    integer, pointer :: rrs(:)

    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)
    allocate(rrs(nlevs-1))
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
    close(un)

     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

    dm = get_dim(mfs(1))

    write(unit=sd_name, fmt='(a,"/State_nodalx")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_nodalx, sd_name)

    write(unit=sd_name, fmt='(a,"/State_nodaly")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_nodaly, sd_name)

    if (dm .eq. 3) then
       write(unit=sd_name, fmt='(a,"/State_nodalz")') trim(dirname)
       call fabio_ml_multifab_read_d(mfs_nodalz, sd_name)
    end if

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

  end subroutine checkpoint_read

end module checkpoint_module
