module write_plotfileLM_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

contains
  
  subroutine write_plotfileLM(mla,name,rho,rhotot,Temp,istep,dx,time)

    type(ml_layout),    intent(in)  :: mla
    character(len=*),   intent(in)  :: name
    type(multifab),     intent(in)  :: rho(:)
    type(multifab),     intent(in)  :: rhotot(:)
    type(multifab),     intent(in)  :: Temp(:)
    integer,            intent(in)  :: istep
    real(kind=dp_t),    intent(in)  :: dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plotfile_name
    integer                         :: n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)

    nlevs = mla%nlevel
  
    allocate(plot_names(nspecies+2)) ! rho + species + Temp
    allocate(plotdata(nlevs))
 
    plot_names(1) = "rho"
    do n=1,nspecies
       write(plot_names(n+1),'(a,i0)') "rho", n
    enddo
    plot_names(nspecies+2) = "Temp"

    ! build plotdata for nspecies+2 and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),nspecies+2,0)
    enddo
    
    ! copy the state into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1         ,rhotot(n),1,       1,0)
       call multifab_copy_c(plotdata(n),2         ,rho(n)   ,1,nspecies,0)
       call multifab_copy_c(plotdata(n),nspecies+2,Temp(n)  ,1,1       ,0)
    enddo
    
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
