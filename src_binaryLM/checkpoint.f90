module checkpoint_module

  use parallel
  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_IO_module
  use fab_module
  use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
  use probin_common_module, only: dim_in
  use probin_binarylm_module, only: algorithm_type

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(mla,sold,mold,pres,rhoc_fluxdiv,time,dt,istep_to_write)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)           ! cell-centered densities
    type(multifab) , intent(in   ) :: mold(:,:)         ! edge-based velocities
    type(multifab) , intent(in   ) :: pres(:)           ! cell-centered pressure
    type(multifab) , intent(in   ) :: rhoc_fluxdiv(:)   ! cell-centered diffusive mass fluxes
    integer        , intent(in   ) :: istep_to_write
    real(kind=dp_t), intent(in   ) :: time,dt

    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chkdata_edge(:,:)

    integer :: n,nlevs,dm,i

    character(len=9) :: sd_name

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(chkdata(nlevs))
    allocate(chkdata_edge(nlevs,dm))
    do n = 1,nlevs
       if (algorithm_type .eq. 0) then
          ! 2 densities + 1 pressure + rhoc mass fluxes
          call multifab_build(chkdata(n), mla%la(n), 4, 0)
       else
          ! 2 densities + 1 pressure
          call multifab_build(chkdata(n), mla%la(n), 3, 0)
       end if
       call multifab_copy_c(chkdata(n), 1, sold(n), 1, 2)  ! copy densities
       call multifab_copy_c(chkdata(n), 3, pres(n), 1, 1)  ! copy pressure
       if (algorithm_type .eq. 0) then
          ! copy mass fluxes
          call multifab_copy_c(chkdata(n), 4, rhoc_fluxdiv(n), 1, 1)
       end if
       do i=1,dm
          ! 1 velocity component and 1 normal bc component for each face
          call multifab_build_edge(chkdata_edge(n,i), mla%la(n), 1, 0, i)
          ! copy velocities
          call multifab_copy_c(chkdata_edge(n,i), 1, mold(n,i), 1, 1)
       end do
    end do
    write(unit=sd_name,fmt='("chk",i6.6)') istep_to_write

    call checkpoint_write_doit(nlevs, sd_name, chkdata, chkdata_edge, mla%mba%rr, time, dt)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
       do i=1,dm
          call multifab_destroy(chkdata_edge(n,i))
       end do
    end do
    deallocate(chkdata)
    deallocate(chkdata_edge)

  contains

    subroutine checkpoint_write_doit(nlevs_in, dirname, mfs, mfs_edge, rrs, time_in, dt_in)
      
      integer         , intent(in) :: nlevs_in
      type(multifab)  , intent(in) :: mfs(:)
      type(multifab)  , intent(in) :: mfs_edge(:,:)
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

      dm = dim_in
      if ( parallel_IOProcessor() ) call fabio_mkdir(dirname)

      call parallel_barrier() ! All CPUs have to wait till the directory is built.

      write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

      if (parallel_IOProcessor()) then
         print *,'Writing cc state to checkpoint file ',trim(sd_name)
         print *,' '
      end if

      write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_edge(:,1), rrs(:,1), sd_name)

      write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_edge(:,2), rrs(:,1), sd_name)

      if (dm .eq. 3) then
         write(unit=sd_name, fmt='(a,"/State_edgez")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_edge(:,3), rrs(:,1), sd_name)
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

    end subroutine checkpoint_write_doit

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_edgex, mfs_edgey, mfs_edgez, &
                             dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    type(multifab)  ,                pointer :: mfs_edgex(:)
    type(multifab)  ,                pointer :: mfs_edgey(:)
    type(multifab)  ,                pointer :: mfs_edgez(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, nlevs, dm
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

    write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgex, sd_name)

    write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgey, sd_name)

    if (dm .eq. 3) then
       write(unit=sd_name, fmt='(a,"/State_edgez")') trim(dirname)
       call fabio_ml_multifab_read_d(mfs_edgez, sd_name)
    end if

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

  end subroutine checkpoint_read

end module checkpoint_module
