module checkpoint_module

  use parallel
  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_IO_module
  use fab_module
  use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
  use probin_common_module, only: dim_in

  implicit none

  private

  public :: write_checkfile, checkpoint_read_2d, checkpoint_read_3d

contains

  subroutine write_checkfile(mla,s,p,umac,vel_bc_n,vel_bc_t,time,dt,istep_to_write)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)           ! cell-centered densities
    type(multifab) , intent(in   ) :: p(:)           ! cell-centered pressure
    type(multifab) , intent(in   ) :: umac(:,:)      ! edge-based velocities
    type(multifab) , intent(in   ) :: vel_bc_n(:,:)  ! edge-based normal bc's
    type(multifab) , intent(in   ) :: vel_bc_t(:,:)  ! in 2D this is nodal
                                                     ! in 3D this is nodal in 2 directions
    integer        , intent(in   ) :: istep_to_write
    real(kind=dp_t), intent(in   ) :: time,dt

    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chkdata_edge(:,:)
    type(multifab), pointer :: chkdata_vel_bc_t(:,:)

    integer :: n,nlevs,dm,i

    character(len=9) :: sd_name

    logical :: nodal_temp(3)

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(chkdata(nlevs))
    allocate(chkdata_edge(nlevs,dm))
    do n = 1,nlevs
       call multifab_build(chkdata(n), mla%la(n), 3, 0) ! 2 densities + 1 pressure
       call multifab_copy_c(chkdata(n), 1, s(n), 1, 2)  ! copy densities
       call multifab_copy_c(chkdata(n), 3, p(n), 1, 1)  ! copy pressure
       do i=1,dm
          ! 1 velocity component and 1 normal bc component for each face
          call multifab_build_edge(chkdata_edge(n,i), mla%la(n), 2, 0, i)
          ! copy velocities
          call multifab_copy_c(chkdata_edge(n,i), 1, umac(n,i), 1, 1)
          ! copy normal bc component
          call multifab_copy_c(chkdata_edge(n,i), 2, vel_bc_n(n,i), 1, 1)
       end do

       if (dm .eq. 2) then
          allocate(chkdata_vel_bc_t(nlevs,2))
          ! y-velocity bc on x-faces (nodal)
          call multifab_build_nodal(chkdata_vel_bc_t(n,1),mla%la(n),1,0)
          call multifab_copy_c(chkdata_vel_bc_t(n,1), 1, vel_bc_t(n,1), 1, 1)
          ! x-velocity bc on y-faces (nodal)
          call multifab_build_nodal(chkdata_vel_bc_t(n,2),mla%la(n),1,0)
          call multifab_copy_c(chkdata_vel_bc_t(n,2), 1, vel_bc_t(n,2), 1, 1)
       else if (dm .eq. 3) then
          allocate(chkdata_vel_bc_t(nlevs,6))
          ! y-velocity bc on x-faces (nodal in y and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(chkdata_vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,1), 1, vel_bc_t(n,1), 1, 1)
          ! z-velocity bc on x-faces (nodal in z and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(chkdata_vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,2), 1, vel_bc_t(n,2), 1, 1)
          ! x-velocity bc on y-faces (nodal in x and y)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(chkdata_vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,3), 1, vel_bc_t(n,3), 1, 1)
          ! z-velocity bc on y-faces (nodal in z and y)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(chkdata_vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,4), 1, vel_bc_t(n,4), 1, 1)
          ! x-velocity bc on z-faces (nodal in x and z)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(chkdata_vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,5), 1, vel_bc_t(n,5), 1, 1)
          ! y-velocity bc on z-faces (nodal in y and z)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(chkdata_vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
          call multifab_copy_c(chkdata_vel_bc_t(n,6), 1, vel_bc_t(n,6), 1, 1)
       end if
    end do
    write(unit=sd_name,fmt='("chk",i6.6)') istep_to_write

    call checkpoint_write(nlevs, sd_name, chkdata, chkdata_edge, chkdata_vel_bc_t, mla%mba%rr, time, dt)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
       do i=1,dm
          call multifab_destroy(chkdata_edge(n,i))
       end do
       do i=1,size(chkdata_vel_bc_t,dim=2)
          call multifab_destroy(chkdata_vel_bc_t(n,i))
       end do
    end do
    deallocate(chkdata)
    deallocate(chkdata_edge)

  contains

    subroutine checkpoint_write(nlevs_in, dirname, mfs, mfs_edge, mfs_vel_bc_t, rrs, time_in, dt_in)
      
      integer         , intent(in) :: nlevs_in
      type(multifab)  , intent(in) :: mfs(:)
      type(multifab)  , intent(in) :: mfs_edge(:,:)
      type(multifab)  , intent(in) :: mfs_vel_bc_t(:,:)
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

      write(unit=sd_name, fmt='(a,"/State_vel_bc_t1")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,1), rrs(:,1), sd_name)

      write(unit=sd_name, fmt='(a,"/State_vel_bc_t2")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,2), rrs(:,1), sd_name)

      if (dm .eq. 3) then

         write(unit=sd_name, fmt='(a,"/State_vel_bc_t3")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,3), rrs(:,1), sd_name)

         write(unit=sd_name, fmt='(a,"/State_vel_bc_t4")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,4), rrs(:,1), sd_name)

         write(unit=sd_name, fmt='(a,"/State_vel_bc_t5")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,5), rrs(:,1), sd_name)

         write(unit=sd_name, fmt='(a,"/State_vel_bc_t6")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_vel_bc_t(:,6), rrs(:,1), sd_name)

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

  end subroutine write_checkfile

  subroutine checkpoint_read_2d(mfs, mfs_edgex, mfs_edgey, mfs_vel_bc_t1, mfs_vel_bc_t2, &
                                dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    type(multifab)  ,                pointer :: mfs_edgex(:)
    type(multifab)  ,                pointer :: mfs_edgey(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t1(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t2(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, nlevs
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

    write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgex, sd_name)

    write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgey, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t1")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t1, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t2")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t2, sd_name)

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

  end subroutine checkpoint_read_2d

  subroutine checkpoint_read_3d(mfs, mfs_edgex, mfs_edgey, mfs_edgez, &
                                mfs_vel_bc_t1, mfs_vel_bc_t2, mfs_vel_bc_t3, &
                                mfs_vel_bc_t4, mfs_vel_bc_t5, mfs_vel_bc_t6, &
                                dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    type(multifab)  ,                pointer :: mfs_edgex(:)
    type(multifab)  ,                pointer :: mfs_edgey(:)
    type(multifab)  ,                pointer :: mfs_edgez(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t1(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t2(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t3(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t4(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t5(:)
    type(multifab)  ,                pointer :: mfs_vel_bc_t6(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, nlevs
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

    write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgex, sd_name)

    write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgey, sd_name)

    write(unit=sd_name, fmt='(a,"/State_edgez")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgez, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t1")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t1, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t2")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t2, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t3")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t3, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t4")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t4, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t5")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t5, sd_name)

    write(unit=sd_name, fmt='(a,"/State_vel_bc_t6")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_vel_bc_t6, sd_name)

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

  end subroutine checkpoint_read_3d

end module checkpoint_module
