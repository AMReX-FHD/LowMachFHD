module write_plotfile_module

  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_module
  use probin_module, only : prob_lo, prob_hi, nscal, center_snapshots, plot_stag
  use fabio_module
  use convert_module

  implicit none

  private

  public :: write_plotfile

contains

  subroutine write_plotfile(mla,m_in,umac_in,s_in,dx,time,istep_to_write)

    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) ::    m_in(:,:)
    type(multifab)   , intent(in   ) :: umac_in(:,:)
    type(multifab)   , intent(in   ) ::    s_in(:)
    integer          , intent(in   ) :: istep_to_write
    real(kind=dp_t)  , intent(in   ) :: dx(:,:), time

    integer                 :: i,n,n_plot_comps
    logical                 :: coarsen_plot_data

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
    type(box)                   :: pd_crse
    type(boxarray)              :: ba_crse
    type(layout)                :: la_crse
    type(multifab), allocatable :: mf_crse(:)
    real(dp_t),     allocatable :: dx_crse(:)
    integer   ,     allocatable :: ref_ratio(:)
    integer                     :: rr_for_write(mla%nlevel-1)
    integer                     :: coarsening_factor
    type(multifab), allocatable :: plotdata(:)          ! cell-centered
    type(multifab), allocatable :: plotdata_stag(:,:)   ! staggered

    integer :: nlevs,dm

    nlevs = mla%nlevel
    dm = mla%dim

    ! Set up plot_names for writing plot files.
    ! velocity (dm) + momentum (dm) + densities (nscal) + concentrations (nscal-1)
    allocate(plot_names(2*dm+2*nscal-1))

    plot_names(1) = "velx"
    plot_names(2) = "vely"
    if (dm > 2) plot_names(3) = "velz"
    plot_names(dm+1) = "mx"
    plot_names(dm+2) = "my"
    if (dm > 2) plot_names(dm+3) = "mz"
    plot_names(2*dm+1) = "rho"
    if (nscal > 1) plot_names(2*dm+2:2*dm+nscal) = "rho*c"
    plot_names(2*dm+nscal+1:2*dm+2*nscal-1) = "c"

    plot_names_stagx(1) = "velx"
    plot_names_stagx(2) = "mx"
    plot_names_stagy(1) = "vely"
    plot_names_stagy(2) = "my"
    plot_names_stagz(1) = "velz"
    plot_names_stagz(2) = "mz"

    allocate(ref_ratio(dm))
    allocate(dx_crse(dm))
    allocate(mf_crse(nlevs))
    !   End crse
  
    coarsen_plot_data = .false.
    coarsening_factor = 2

    allocate(plotdata(nlevs))
    allocate(plotdata_stag(nlevs,dm))

    ! velocity (dm) + momentum (dm) + densities (nscal) + concentrations (nscal-1)
    n_plot_comps = 2*dm+2*nscal-1

    do n=1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       do i=1,dm
          call multifab_build_edge(plotdata_stag(n,i), mla%la(n), 2, 0, i)
       end do
    end do

    ! copy staggered velocity and momentum into plotdata_stag
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(plotdata_stag(n,i),1,umac_in(n,i),1,1)
          call multifab_copy_c(plotdata_stag(n,i),2,   m_in(n,i),1,1)
       end do
    end do

    if(center_snapshots) then
       ! Average from faces to cell centers and render the centered (but also smoothed) vector field
       ! We do this for momentum only, and then divide by density, rather than averaging umac
       ! This ensures the momentum is strictly conserved by the interpolation operation

       do i=1,dm          
          call average_face_to_cc(mla,m_in(:,i),1,plotdata,dm+i)
       end do   

       do n=1,nlevs ! Calculate u=m/rho
          call multifab_copy_c(plotdata(n),1,plotdata(n),dm+1,dm)
          do i=1,dm
             call multifab_div_div_c(plotdata(n),i,s_in(n),1,1)
          end do   
       end do      

    else ! just shift face data to cell centers (OK for periodic)

       do i=1,dm
          call shift_face_to_cc(mla,umac_in(:,i),1,plotdata,i)
          call shift_face_to_cc(mla,m_in(:,i),1,plotdata,dm+i)
       end do

    end if   

    do n = 1,nlevs
       ! rho and rho*c
       call multifab_copy_c(plotdata(n),2*dm+1,s_in(n),1,nscal)
       ! c
       call multifab_copy_c(plotdata(n),2*dm+1+nscal,s_in(n),2,nscal-1)
       do i=1,nscal-1
          call multifab_div_div_c(plotdata(n),2*dm+nscal+i,s_in(n),1,1,0)
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

    if (coarsen_plot_data) then

       ! We have only implemented this for nlevs = 1 right now
       ref_ratio(1:dm) = coarsening_factor
       rr_for_write(:) = coarsening_factor
       dx_crse(:) = dx(1,:) / ref_ratio

       pd_crse = coarsen(mla%mba%pd(1),ref_ratio)

       do n = 1, nlevs
          call boxarray_build_copy(ba_crse,get_boxarray(mla%la(n)))
          call boxarray_coarsen(ba_crse,ref_ratio)
          call layout_build_ba(la_crse,ba_crse,coarsen(mla%mba%pd(n),ref_ratio))
          call print(mla%la(n),'LA FINE')
          call print(la_crse,'LA CRSE')
          call multifab_build(mf_crse(n), la_crse, n_plot_comps, 0)
          call destroy(ba_crse)

          call ml_cc_restriction(mf_crse(n),plotdata(n),ref_ratio)
       end do

       call fabio_ml_multifab_write_d(mf_crse, rr_for_write, sd_name, plot_names, &
                                      pd_crse, prob_lo, prob_hi, time, dx_crse)
    else

       call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                      mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))

       if (plot_stag) then
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
       end if

    end if

    do n=1,nlevs
      call multifab_destroy(plotdata(n))
      do i=1,dm
         call multifab_destroy(plotdata_stag(n,i))
      end do
    end do
    deallocate(plotdata)
    deallocate(plotdata_stag)

    if (coarsen_plot_data) then
       do n = 1, nlevs
          call destroy(mf_crse(n))
       end do
       deallocate(mf_crse)
       deallocate(ref_ratio)
       deallocate(dx_crse)
    end if

    deallocate(plot_names)

  end subroutine write_plotfile

end module write_plotfile_module
