module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_rhoh_to_h_module
  use eos_model_wrapper_module
  use probin_multispecies_module, only: nspecies, plot_stag
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

contains
  
  subroutine write_plotfile(mla,rho,rhotot,rhoh,Temp,umac,pi,p0,istep,dx,time)

    type(ml_layout),    intent(in)    :: mla
    type(multifab),     intent(inout) :: rho(:)
    type(multifab),     intent(in)    :: rhotot(:)
    type(multifab),     intent(inout) :: rhoh(:)
    type(multifab),     intent(in)    :: Temp(:)
    type(multifab),     intent(in)    :: umac(:,:)
    type(multifab),     intent(in)    :: pi(:)
    integer,            intent(in)    :: istep
    real(kind=dp_t),    intent(in)    :: p0,dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plot_names_stagx(1)
    character(len=20)               :: plot_names_stagy(1)
    character(len=20)               :: plot_names_stagz(1)
    character(len=20)               :: plotfile_name
    character(len=20)               :: plotfile_namex
    character(len=20)               :: plotfile_namey
    character(len=20)               :: plotfile_namez
    integer                         :: i,dm,n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)
    type(multifab), allocatable     :: plotdata_stag(:,:)

    type(multifab) :: conc(mla%nlevel)
    type(multifab) :: h(mla%nlevel)

    ! pressure from the EOS
    type(multifab) :: Peos(mla%nlevel)

    nlevs = mla%nlevel
    dm = mla%dim
  
    ! cell-centered quantities

    ! rho
    ! rho_i
    ! c_i
    ! rhoh
    ! h
    ! Temp
    ! umac averaged
    ! umac shifted  
    ! pi
    ! p0
    ! pi
    ! Peos
    ! (Peos-p0)
    allocate(plot_names(2*nspecies+2*dm+8))
    allocate(plotdata(nlevs))
    allocate(plotdata_stag(nlevs,dm))
 
    plot_names(1) = "rho"
    do n=1,nspecies
       write(plot_names(n+1),'(a,i0)') "rho", n
    enddo
    do n=1,nspecies
       write(plot_names(nspecies+n+1),'(a,i0)') "c", n
    enddo
    plot_names(2*nspecies+2) = "rhoh"
    plot_names(2*nspecies+3) = "h"
    plot_names(2*nspecies+4) = "Temp"
    plot_names(2*nspecies+5) = "averaged_velx"
    plot_names(2*nspecies+6) = "averaged_vely"
    if (dm > 2) plot_names(2*nspecies+7) = "averaged_velz"
    plot_names(2*nspecies+dm+5) = "shifted_velx"
    plot_names(2*nspecies+dm+6) = "shifted_vely"
    if (dm > 2) plot_names(2*nspecies+dm+7) = "shifted_velz"
    plot_names(2*nspecies+2*dm+5) = "pi"
    plot_names(2*nspecies+2*dm+6) = "p0"
    plot_names(2*nspecies+2*dm+7) = "Peos"
    plot_names(2*nspecies+2*dm+8) = "Peos_minus_p0"

    plot_names_stagx(1) = "velx"
    plot_names_stagy(1) = "vely"
    plot_names_stagz(1) = "velz"

    ! for Peos
    do n=1,nlevs
       call multifab_build(Peos(n),mla%la(n),1,0)
    end do

    ! compute concentrations
    do n=1,nlevs
       call multifab_build(conc(n),mla%la(n),nspecies,0)
    end do
    call convert_rhoc_to_c(mla,rho,rhotot,conc,.true.)

    ! compute h
    do n=1,nlevs
       call multifab_build(h(n),mla%la(n),1,0)
    end do
    call convert_rhoh_to_h(mla,rhoh,rhotot,h,.true.)

    ! build plotdata for 2*nspecies+2*dm+8 and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),2*nspecies+2*dm+8,0)
       do i=1,dm
          call multifab_build_edge(plotdata_stag(n,i), mla%la(n), 1, 0, i)
       end do
    enddo
    
    ! copy rhotot, rho, conc, rhoh, h, and Temp into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1           ,rhotot(n),1,       1,0)
       call multifab_copy_c(plotdata(n),2           ,rho(n)   ,1,nspecies,0)
       call multifab_copy_c(plotdata(n),nspecies+2  ,conc(n)  ,1,nspecies,0)
       call multifab_copy_c(plotdata(n),2*nspecies+2,rhoh(n)  ,1,1       ,0)
       call multifab_copy_c(plotdata(n),2*nspecies+3,h(n)     ,1,1       ,0)
       call multifab_copy_c(plotdata(n),2*nspecies+4,Temp(n)  ,1,1       ,0)
    enddo

    ! vel averaged
    do i=1,dm
       call average_face_to_cc(mla,umac(:,i),1,plotdata,2*nspecies+4+i,1)
    end do

    ! vel shifted
    do i=1,dm
       call shift_face_to_cc(mla,umac(:,i),1,plotdata,2*nspecies+dm+4+i,1)
    end do

    ! pi
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),2*nspecies+2*dm+5,pi(n),1,1,0)
    enddo

    ! p0
    do n=1,nlevs
       call multifab_setval_c(plotdata(n),p0,2*nspecies+2*dm+6,1)
    end do

    ! compute P_eos
    call compute_p(mla,rhotot,Temp,conc,Peos)
    do n = 1,nlevs
       ! Peos
       call multifab_copy_c(plotdata(n),2*nspecies+2*dm+7,Peos(n),1,1,0)
       call multifab_sub_sub_s_c(Peos(n),1,p0,1,0)
       ! Peos - p0
       call multifab_copy_c(plotdata(n),2*nspecies+2*dm+8,Peos(n),1,1,0)
    enddo

    ! copy staggered velocity and momentum into plotdata_stag
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(plotdata_stag(n,i),1,umac(n,i),1,1)
       end do
    end do
    
    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='(a,i8.8)') trim(plot_base_name), istep
    if ( parallel_IOProcessor() ) then
      write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(plotfile_name)
      write(*,*)
    end if
    write(unit=plotfile_namex,fmt='("stagx",i8.8)') istep
    write(unit=plotfile_namey,fmt='("stagy",i8.8)') istep
    if (dm > 2) then
       write(unit=plotfile_namez,fmt='("stagz",i8.8)') istep
    end if
    
    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), trim(plotfile_name), & 
                                   plot_names, mla%mba%pd(1), prob_lo, prob_hi, & 
                                   time, dx(1,1:dm))

    if (plot_stag) then
       ! staggered plotfiles
       ! the data appears "shifted" in amrvis, but you can post process these and still keep
       ! all the boundary data
       call fabio_ml_multifab_write_d(plotdata_stag(:,1), mla%mba%rr(:,1), trim(plotfile_namex), &
                                      plot_names_stagx, mla%mba%pd(1), prob_lo, prob_hi, &
                                      time, dx(1,1:dm))

       call fabio_ml_multifab_write_d(plotdata_stag(:,2), mla%mba%rr(:,1), trim(plotfile_namey), &
                                      plot_names_stagy, mla%mba%pd(1), prob_lo, prob_hi, &
                                      time, dx(1,1:dm))

       if (dm > 2) then
          call fabio_ml_multifab_write_d(plotdata_stag(:,3), mla%mba%rr(:,1), trim(plotfile_namez), &
                                         plot_names_stagz, mla%mba%pd(1), prob_lo, prob_hi, &
                                         time, dx(1,1:dm))
       end if
    end if

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(Peos(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(h(n))
       call multifab_destroy(plotdata(n))
      do i=1,dm
         call multifab_destroy(plotdata_stag(n,i))
      end do
    enddo
    
    deallocate(plotdata)
    deallocate(plotdata_stag)
    deallocate(plot_names)

  end subroutine write_plotfile

end module write_plotfile_module
