module write_plotfile1_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use probin_multispecies_module

  implicit none

contains
  
  subroutine write_plotfile1(mla,name,Temp,istep,dx,time,prob_lo,prob_hi)

    type(ml_layout),    intent(in)  :: mla
    character(len=*),   intent(in)  :: name
    type(multifab),     intent(in)  :: Temp(:)
    integer,            intent(in)  :: istep
    real(kind=dp_t),    intent(in)  :: dx(:,:),time
    real(kind=dp_t),    intent(in)  :: prob_lo(Temp(1)%dim), prob_hi(Temp(1)%dim)

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plotfile_name
    integer                         :: n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)

    nlevs = mla%nlevel
  
    allocate(plot_names(1))
    allocate(plotdata(nlevs))
 
    ! write density with species index
    write(plot_names(1),'(a,i0)') "Scalar"

    ! build plotdata for nspecies and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),1,0)
    enddo
    
    ! copy the state into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1,Temp(n),1,1,0)
    enddo
    
    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='(a,i5.5)') name, istep
        
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

  end subroutine write_plotfile1

end module write_plotfile1_module
