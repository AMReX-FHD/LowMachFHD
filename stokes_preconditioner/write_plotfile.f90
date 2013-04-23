module write_plotfile_module

  use ml_layout_module
  use ml_restriction_module
  use multifab_module
  use probin_module, only : prob_lo, prob_hi
  use fabio_module
  use convert_module

  implicit none

  private

  public :: write_plotfile

contains

  subroutine write_plotfile(mla,umac,pres,alpha,beta,gamma,dx,time,istep)

    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: umac(:,:)
    type(multifab)   , intent(in   ) :: pres(:)
    type(multifab)   , intent(in   ) :: alpha(:)
    type(multifab)   , intent(in   ) :: beta(:)
    type(multifab)   , intent(in   ) :: gamma(:)
    real(kind=dp_t)  , intent(in   ) :: dx(:,:), time
    integer          , intent(in   ) :: istep

    integer                 :: i,n,nlevs,dm

    character(len=9 ) :: sd_name
    character(len=20) :: plot_names(2*mla%dim+4) ! Set up plot_names for writing plot files.

    type(multifab)    :: plotdata(mla%nlevel)    ! cell-centered

    nlevs = mla%nlevel
    dm = mla%dim

    plot_names(1) = "velx"
    plot_names(2) = "vely"
    if (dm > 2) plot_names(3) = "velz"
    plot_names(dm+1) = "pressure"
    plot_names(dm+2) = "alpha"
    plot_names(dm+3) = "beta"
    plot_names(dm+4) = "gamma"

    do n=1,nlevs
       call multifab_build(plotdata(n), mla%la(n), dm+4, 0)
    end do

    ! average mac velocities to cell centers
    do i=1,dm          
       call average_face_to_cc(mla,umac(:,i),1,plotdata,i)
    end do

    ! pressure
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),dm+1,pres(n),1,1)
    end do

    ! alpha, beta, and gamma
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),dm+2,alpha(n),1,1)
       call multifab_copy_c(plotdata(n),dm+3,beta(n),1,1)
       call multifab_copy_c(plotdata(n),dm+4,gamma(n),1,1)
    end do

    write(unit=sd_name,fmt='("plt",i6.6)') istep
    if ( parallel_IOProcessor() ) then
      write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(sd_name)
      write(*,*)
    end if

    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))

    do n=1,nlevs
      call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile

end module write_plotfile_module
