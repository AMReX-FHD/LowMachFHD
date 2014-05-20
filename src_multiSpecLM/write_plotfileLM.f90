module write_plotfileLM_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use convert_stag_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

contains
  
  subroutine write_plotfileLM(mla,name,rho,rhotot,Temp,umac,istep,dx,time)

    type(ml_layout),    intent(in)  :: mla
    character(len=*),   intent(in)  :: name
    type(multifab),     intent(in)  :: rho(:)
    type(multifab),     intent(in)  :: rhotot(:)
    type(multifab),     intent(in)  :: umac(:,:)
    type(multifab),     intent(in)  :: Temp(:)
    integer,            intent(in)  :: istep
    real(kind=dp_t),    intent(in)  :: dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plotfile_name
    integer                         :: i,dm,n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)

    nlevs = mla%nlevel
    dm = mla%dim
  
    allocate(plot_names(nspecies+2*dm+2)) ! rho + species + Temp + dm (averaged umac) + dm (shifted umac)
    allocate(plotdata(nlevs))
 
    plot_names(1) = "rho"
    do n=1,nspecies
       write(plot_names(n+1),'(a,i0)') "rho", n
    enddo
    plot_names(nspecies+2) = "Temp"
    plot_names(nspecies+3) = "averaged_velx"
    plot_names(nspecies+4) = "averaged_vely"
    if (dm > 2) plot_names(nspecies+5) = "averaged_velz"
    plot_names(nspecies+dm+3) = "shifted_velx"
    plot_names(nspecies+dm+4) = "shifted_vely"
    if (dm > 2) plot_names(nspecies+dm+5) = "shifted_vely"

    ! build plotdata for nspecies+2*dm+2 and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),nspecies+2*dm+2,0)
    enddo
    
    ! copy rhotot, rho, and Temp into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1         ,rhotot(n),1,       1,0)
       call multifab_copy_c(plotdata(n),2         ,rho(n)   ,1,nspecies,0)
       call multifab_copy_c(plotdata(n),nspecies+2,Temp(n)  ,1,1       ,0)
    enddo

    ! vel averaged
    do i=1,dm
       call average_face_to_cc(mla,umac(:,i),1,plotdata,nspecies+2+i,1)
    end do

    ! vel shifted
    do i=1,dm
       call shift_face_to_cc(mla,umac(:,i),1,plotdata,nspecies+dm+2+i,1)
    end do
    
    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='(a,i6.6)') name, istep
    !if ( parallel_IOProcessor() ) then
    !  write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(plotfile_name)
    !  write(*,*)
    !end if
    
    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), trim(plotfile_name), & 
                                   plot_names, mla%mba%pd(1), prob_lo, prob_hi, & 
                                   time, dx(1,:))

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    enddo
    
    deallocate(plotdata)
    deallocate(plot_names)

  end subroutine write_plotfileLM

end module write_plotfileLM_module
