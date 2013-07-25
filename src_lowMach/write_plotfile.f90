module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use probin_lowmach_module, only : nscal
  use probin_common_module , only : prob_lo, prob_hi
  use fabio_module
  use convert_stag_module

  implicit none

  private

  public :: write_plotfile

contains

  subroutine write_plotfile(mla,m_in,umac_in,s_in,p_in,dx,time,istep_to_write)

    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) ::    m_in(:,:)
    type(multifab)   , intent(in   ) :: umac_in(:,:)
    type(multifab)   , intent(in   ) ::    s_in(:)
    type(multifab)   , intent(in   ) ::    p_in(:)
    integer          , intent(in   ) :: istep_to_write
    real(kind=dp_t)  , intent(in   ) :: dx(:,:), time

    integer                 :: i,n,n_plot_comps

    character(len=9 ) :: sd_name
    character(len=11) :: sd_namex
    character(len=11) :: sd_namey
    character(len=11) :: sd_namez

    character(len=20), allocatable :: plot_names(:)
    character(len=20)              :: plot_names_stagx(2)
    character(len=20)              :: plot_names_stagy(2)
    character(len=20)              :: plot_names_stagz(2)

    ! These are only used if you want to coarsen your plotdata before writing
    ! Start crse
    type(multifab), allocatable :: plotdata(:)          ! cell-centered
    type(multifab), allocatable :: plotdata_stag(:,:)   ! staggered

    integer :: nlevs,dm

    nlevs = mla%nlevel
    dm = mla%dim

    ! Set up plot_names for writing plot files.
    ! velocity_averaged (dm)
    ! velocity_shifted (dm)
    ! momentum_averaged (dm) 
    ! momentum_shifted (dm)
    ! densities (nscal) 
    ! concentrations (nscal-1)
    ! pressure (1)
    n_plot_comps = 4*dm+2*nscal

    allocate(plot_names(n_plot_comps))

    plot_names(1) = "velx_averaged"
    plot_names(2) = "vely_averaged"
    if (dm > 2) plot_names(3) = "velz_averaged"

    plot_names(dm+1) = "velx_shifted"
    plot_names(dm+2) = "vely_shifted"
    if (dm > 2) plot_names(dm+3) = "velz_shifted"

    plot_names(2*dm+1) = "mx_averaged"
    plot_names(2*dm+2) = "my_averaged"
    if (dm > 2) plot_names(2*dm+3) = "mz_averaged"

    plot_names(3*dm+1) = "mx_shifted"
    plot_names(3*dm+2) = "my_shifted"
    if (dm > 2) plot_names(3*dm+3) = "mz_shifted"

    plot_names(4*dm+1) = "rho"
    plot_names(4*dm+2:4*dm+nscal) = "rho*c"
    plot_names(4*dm+nscal+1:4*dm+2*nscal-1) = "c"

    plot_names(4*dm+2*nscal) = "pressure"

    plot_names_stagx(1) = "velx"
    plot_names_stagx(2) = "mx"
    plot_names_stagy(1) = "vely"
    plot_names_stagy(2) = "my"
    plot_names_stagz(1) = "velz"
    plot_names_stagz(2) = "mz"

    allocate(plotdata(nlevs))
    allocate(plotdata_stag(nlevs,dm))

    do n=1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       do i=1,dm
          call multifab_build_edge(plotdata_stag(n,i), mla%la(n), 2, 0, i)
       end do
    end do

    ! vel averaged
    do i=1,dm
       call average_face_to_cc(mla,umac_in(:,i),1,plotdata,i,1)
    end do

    ! vel shifted
    do i=1,dm
       call shift_face_to_cc(mla,umac_in(:,i),1,plotdata,dm+i,1)
    end do

    ! m averaged
    do i=1,dm
       call average_face_to_cc(mla,m_in(:,i),1,plotdata,2*dm+i,1)
    end do

    ! m shifted
    do i=1,dm
       call shift_face_to_cc(mla,m_in(:,i),1,plotdata,3*dm+i,1)
    end do

    do n = 1,nlevs
       ! densities
       call multifab_copy_c(plotdata(n),4*dm+1,s_in(n),1,nscal)
       ! concentrations
       call multifab_copy_c(plotdata(n),4*dm+1+nscal,s_in(n),2,nscal-1)
       do i=1,nscal-1
          call multifab_div_div_c(plotdata(n),4*dm+nscal+i,s_in(n),1,1,0)
       end do
       ! pressure
       call multifab_copy_c(plotdata(n),4*dm+2*nscal,p_in(n),1,1,0)
    end do

    ! copy staggered velocity and momentum into plotdata_stag
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(plotdata_stag(n,i),1,umac_in(n,i),1,1)
          call multifab_copy_c(plotdata_stag(n,i),2,   m_in(n,i),1,1)
       end do
    end do
    
    write(unit=sd_name,fmt='("plt",i6.6)') istep_to_write
    if ( parallel_IOProcessor() ) then
      write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(sd_name)
      write(*,*)
    end if
    write(unit=sd_namex,fmt='("stagx",i6.6)') istep_to_write
    write(unit=sd_namey,fmt='("stagy",i6.6)') istep_to_write
    if (dm > 2) then
       write(unit=sd_namez,fmt='("stagz",i6.6)') istep_to_write
    end if

    ! cell-centered plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))

    ! staggered plotfiles
    ! the data appears "shifted" in amrvis, but you can post process these and still keep
    ! all the boundary data
    call fabio_ml_multifab_write_d(plotdata_stag(:,1), mla%mba%rr(:,1), sd_namex, &
                                   plot_names_stagx, mla%mba%pd(1), prob_lo, prob_hi, &
                                   time, dx(1,:))

    call fabio_ml_multifab_write_d(plotdata_stag(:,2), mla%mba%rr(:,1), sd_namey, &
                                   plot_names_stagy, mla%mba%pd(1), prob_lo, prob_hi, &
                                   time, dx(1,:))

    if (dm > 2) then
       call fabio_ml_multifab_write_d(plotdata_stag(:,3), mla%mba%rr(:,1), sd_namez, &
                                      plot_names_stagz, mla%mba%pd(1), prob_lo, prob_hi, &
                                      time, dx(1,:))
    end if

    do n=1,nlevs
      call multifab_destroy(plotdata(n))
      do i=1,dm
         call multifab_destroy(plotdata_stag(n,i))
      end do
    end do
    deallocate(plotdata)
    deallocate(plotdata_stag)
    deallocate(plot_names)

  end subroutine write_plotfile

end module write_plotfile_module
