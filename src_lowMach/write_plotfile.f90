module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use fluid_charge_module
  use probin_multispecies_module, only: plot_stag
  use probin_common_module, only: prob_lo, prob_hi, nspecies, plot_base_name
  use probin_charged_module, only: use_charged_fluid

  implicit none

contains
  
  subroutine write_plotfile(mla,rho,rhotot,Temp,umac,pres,Epot,grad_Epot, &
                            gradPhiApprox,istep,dx,time)

    type(ml_layout),    intent(in)    :: mla
    type(multifab),     intent(inout) :: rho(:)
    type(multifab),     intent(in)    :: rhotot(:)
    type(multifab),     intent(in)    :: Temp(:)
    type(multifab),     intent(in)    :: umac(:,:)
    type(multifab),     intent(in)    :: pres(:)
    type(multifab),     intent(in)    :: Epot(:)
    type(multifab),     intent(in)    :: grad_Epot(:,:)
    type(multifab),     intent(in)    :: gradPhiApprox(:,:)
    integer,            intent(in)    :: istep
    real(kind=dp_t),    intent(in)    :: dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20), allocatable  :: plot_names_stagx(:)
    character(len=20), allocatable  :: plot_names_stagy(:)
    character(len=20), allocatable  :: plot_names_stagz(:)
    character(len=128)              :: plotfile_name
    character(len=128)              :: plotfile_namex
    character(len=128)              :: plotfile_namey
    character(len=128)              :: plotfile_namez
    integer                         :: i,dm,n,nlevs

    ! multifab of size nlevs  
    type(multifab) :: plotdata(mla%nlevel)
    type(multifab) :: plotdata_stag(mla%nlevel,mla%dim)

    type(multifab) :: cc_temp(mla%nlevel)

    integer :: nvarsCC, nvarsStag, counter

    nlevs = mla%nlevel
    dm = mla%dim

    ! cell-centered quantities

    ! rho
    ! rho_i
    ! c_i
    ! Temp
    ! umac averaged
    ! umac shifted
    ! pressure
    nvarsCC = 2*nspecies + 2*dm + 3

    if (use_charged_fluid) then
       ! charge
       ! Epot
       ! grad_Epot averaged
       ! gradPhiApprox averaged
       nvarsCC = nvarsCC + 2 + 2*dm
    end if

    allocate(plot_names(nvarsCC))
 
    counter = 1

    plot_names(counter) = "rho"
    counter = counter + 1
    do n=1,nspecies
       write(plot_names(counter),'(a,i0)') "rho", n
       counter = counter + 1
    enddo
    do n=1,nspecies
       write(plot_names(counter),'(a,i0)') "c", n
       counter = counter + 1
    enddo

    plot_names(counter) = "Temp"
    counter = counter + 1

    plot_names(counter) = "averaged_velx"
    counter = counter + 1
    plot_names(counter) = "averaged_vely"
    counter = counter + 1
    if (dm > 2) then
       plot_names(counter) = "averaged_velz"
       counter = counter + 1
    end if

    plot_names(counter) = "shifted_velx"
    counter = counter + 1
    plot_names(counter) = "shifted_vely"
    counter = counter + 1
    if (dm > 2) then
       plot_names(counter) = "shifted_velz"
       counter = counter + 1
    end if

    plot_names(counter) = "pres"
    counter = counter + 1

    if (use_charged_fluid) then
       plot_names(counter) = "charge_density"
       counter = counter + 1

       plot_names(counter) = "Epot"
       counter = counter + 1

       plot_names(counter) = "av_gradEpotx"
       counter = counter + 1
       plot_names(counter) = "av_gradEpoty"
       counter = counter + 1
       if (dm > 2) then
          plot_names(counter) = "av_gradEpotz"
          counter = counter + 1
       end if

       plot_names(counter) = "av_gradPhiApproxx"
       counter = counter + 1
       plot_names(counter) = "av_gradPhiApproxy"
       counter = counter + 1
       if (dm > 2) then
          plot_names(counter) = "av_gradPhiApproxz"
          counter = counter + 1
       end if
    end if

    ! staggered quantities

    ! vel
    nvarsStag = 1
    if (use_charged_fluid) then
       ! grad_Epot
       ! gradPhiApprox
       nvarsStag = nvarsStag + 2
    end if

    allocate(plot_names_stagx(nvarsStag))
    allocate(plot_names_stagy(nvarsStag))
    allocate(plot_names_stagz(nvarsStag))

    counter = 1

    plot_names_stagx(counter) = "velx"
    plot_names_stagy(counter) = "vely"
    plot_names_stagz(counter) = "velz"
    counter = counter + 1

    if (use_charged_fluid) then
       plot_names_stagx(counter) = "grad_Epotx"
       plot_names_stagy(counter) = "grad_Epoty"
       plot_names_stagz(counter) = "grad_Epotz"
       counter = counter + 1

       plot_names_stagx(counter) = "gradPhiApproxx"
       plot_names_stagy(counter) = "gradPhiApproxy"
       plot_names_stagz(counter) = "gradPhiApproxz"
       counter = counter + 1
    end if

    do n=1,nlevs
       ! temporary to help with variable conversions
       call multifab_build(cc_temp(n),mla%la(n),nspecies,0)
       ! build plotdata for nvarsCC components
       call multifab_build(plotdata(n),mla%la(n),nvarsCC,0)
       do i=1,dm
          ! plotdata_stag for nvarsStag components
          call multifab_build_edge(plotdata_stag(n,i), mla%la(n), nvarsStag, 0, i)
       end do
    enddo

    counter = 1
    
    ! rhotot
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),counter,rhotot(n),1,1,0)
    end do
    counter = counter + 1

    ! rho
    do n=1,nlevs
       call multifab_copy_c(plotdata(n),counter,rho(n),1,nspecies,0)
    end do
    counter = counter + nspecies

    ! compute concentrations
    call convert_rhoc_to_c(mla,rho,rhotot,cc_temp,.true.)
    do n=1,nlevs
       call multifab_copy_c(plotdata(n),counter,cc_temp(n),1,nspecies,0)
    end do
    counter = counter + nspecies

    ! Temp
    do n=1,nlevs
       call multifab_copy_c(plotdata(n),counter,Temp(n),1,1,0)
    enddo
    counter = counter + 1

    ! vel averaged
    do i=1,dm
       call average_face_to_cc(mla,umac(:,i),1,plotdata,counter,1)
       counter = counter + 1
    end do

    ! vel shifted
    do i=1,dm
       call shift_face_to_cc(mla,umac(:,i),1,plotdata,counter,1)
       counter = counter + 1
    end do

    ! pressure
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),counter,pres(n),1,1,0)
    enddo
    counter = counter + 1

    if (use_charged_fluid) then

       ! compute total charge, then copy into the correct component
       call dot_with_z(mla,rho,cc_temp)
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),counter,cc_temp(n),1,1,0)
       end do
       counter = counter + 1

       ! Epot
       do n = 1,nlevs
          call multifab_copy_c(plotdata(n),counter,Epot(n),1,1,0)
       enddo
       counter = counter + 1

       ! grad_Epot_averaged
       do i=1,dm
          call average_face_to_cc(mla,grad_Epot(:,i),1,plotdata,counter,1)
          counter = counter + 1
       end do

       ! gradPhiApprox_averaged
       do i=1,dm
          call average_face_to_cc(mla,gradPhiApprox(:,i),1,plotdata,counter,1)
          counter = counter + 1
       end do
    end if

    counter = 1

    ! copy staggered velocity
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(plotdata_stag(n,i),counter,umac(n,i),1,1)
       end do
    end do
    counter = counter + 1

    if (use_charged_fluid) then

       ! grad_Epot
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(plotdata_stag(n,i),counter,grad_Epot(n,i),1,1)
          end do
       end do
       counter = counter + 1

       ! gradPhiApprox
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(plotdata_stag(n,i),counter,gradPhiApprox(n,i),1,1)
          end do
       end do
       counter = counter + 1

    end if
    
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
       call multifab_destroy(cc_temp(n))
       call multifab_destroy(plotdata(n))
      do i=1,dm
         call multifab_destroy(plotdata_stag(n,i))
      end do
    enddo
    
    deallocate(plot_names)
    deallocate(plot_names_stagx)
    deallocate(plot_names_stagy)
    deallocate(plot_names_stagz)

  end subroutine write_plotfile

end module write_plotfile_module
