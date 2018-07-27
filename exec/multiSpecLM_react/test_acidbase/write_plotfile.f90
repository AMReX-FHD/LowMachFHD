! this is a custom version of write plotfile for problems with the 4 species, in order:
! HCl NaOH NaCl H2O
! and writes out the the plotfile the ionic concentraions (for rho_i and c_i), in order:
! Na+ Cl- H+ OH- H2O

module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use fluid_charge_module
  use eos_check_module
  use probin_multispecies_module, only: plot_stag, is_nonisothermal
  use probin_common_module, only: prob_lo, prob_hi, nspecies, plot_base_name, &
                                  algorithm_type, rho0, plot_umac_tavg, plot_Epot_tavg, & 
                                  plot_rho_tavg, plot_avg_gradPhiApprox, plot_shifted_vel, &
                                  plot_gradEpot
  use probin_charged_module, only: use_charged_fluid

  implicit none

contains
  
  subroutine write_plotfile(mla,rho,rho_avg,rhotot,Temp,umac,umac_avg,pres,Epot,Epot_avg, &  
                              grad_Epot, gradPhiApprox,istep,dx,time)

    type(ml_layout),    intent(in)    :: mla
    type(multifab),     intent(inout) :: rho(:)
    type(multifab),     intent(inout) :: rho_avg(:) 
    type(multifab),     intent(in)    :: rhotot(:)
    type(multifab),     intent(in)    :: Temp(:)
    type(multifab),     intent(in)    :: umac(:,:)
    type(multifab),     intent(in)    :: umac_avg(:,:)
    type(multifab),     intent(in)    :: pres(:)
    type(multifab),     intent(in)    :: Epot(:)
    type(multifab),     intent(in)    :: Epot_avg(:) 
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

    integer :: nvarsCC, nvarsStag, counter, nspecies_eff
    logical :: boussinesq

    real(kind=dp_t) :: massfrac
    real(kind=dp_t), parameter :: m_Na=3.817540700000000d-023, &	
                                  m_Cl=5.887108600000000d-023, &
                                  m_H =1.673723600000000d-024, &
                                  m_OH=2.824068560000000d-023

    nlevs = mla%nlevel
    dm = mla%dim

    boussinesq = (algorithm_type==6).or.(algorithm_type==4)

    ! split HCl NaOH NaCl H2O -> Na Cl H OH H2O
    if (nspecies .eq. 4) then
       nspecies_eff = nspecies + 1
    else
       nspecies_eff = nspecies
    end if

    ! cell-centered quantities

    ! rho                  :1
    ! rho_i                :nspecies (partial densities)
    ! rho_avg_i            :nspecies
    ! c_i                  :nspecies (mass fractions, only if not Boussinesq)
    ! Temp                 :1 (temperature, only if non isothermal)
    ! umac averaged        :dm
    ! umac shifted         :dm (optional)
    ! umac_avg averaged    :dm (optional)
    ! umac_avg shifted     :dm (optional)
    ! pressure             :1 (pressure)
    if(boussinesq) then
       ! For Boussinesq we don't have to write mass fractions since density is constant
       nvarsCC = nspecies_eff + dm + 2    
    else 
       nvarsCC = 2*nspecies_eff + dm + 2
    end if  
    if(is_nonisothermal) then ! Add temperature
       nvarsCC = nvarsCC + 1
    end if 
    if(plot_shifted_vel) then
       nvarsCC = nvarsCC + dm
    end if
    if (plot_rho_tavg) then  ! time-averaged rho
       nvarsCC = nvarsCC + nspecies
    end if 
    if (plot_umac_tavg) then ! time-averaged umac (cc and shifted)
       nvarsCC = nvarsCC + dm
       if(plot_shifted_vel) then
          nvarsCC = nvarsCC + dm
       end if
    end if

    if (use_charged_fluid) then
       ! charge                   :1
       ! Epot                     :1
       ! Epot_avg                 :1  (optional)
       ! grad_Epot averaged       :dm (optional)
       ! gradPhiApprox averaged   :dm (optional)

       nvarsCC = nvarsCC + 2
       if (plot_gradEpot) then
          nvarsCC = nvarsCC + dm
       end if
       if (plot_Epot_tavg) then 
          nvarsCC = nvarsCC + 1
       end if 
       if (plot_avg_gradPhiApprox) then
          nvarsCC = nvarsCC + dm
       end if
    end if

    if(boussinesq) then ! Boussinesq
       ! add rho_eos-rho0
       nvarsCC = nvarsCC+1
    end if

    allocate(plot_names(nvarsCC))

    counter = 1

    plot_names(counter) = "rho"
    counter = counter + 1
    do n=1,nspecies_eff
       write(plot_names(counter),'(a,i0)') "rho", n
       counter = counter + 1
    enddo
    if (plot_rho_tavg) then
       call bl_error("this local write_plotfile.f90 does not support plot_rho_tavg")
    endif
    if (.not.boussinesq) then
       do n=1,nspecies_eff
          write(plot_names(counter),'(a,i0)') "c", n
          counter = counter + 1
       enddo
    end if

    if(is_nonisothermal) then  
       plot_names(counter) = "Temp"
       counter = counter + 1
    end if   

    plot_names(counter) = "averaged_velx"
    counter = counter + 1
    plot_names(counter) = "averaged_vely"
    counter = counter + 1
    if (dm > 2) then
       plot_names(counter) = "averaged_velz"
       counter = counter + 1
    end if

    if(plot_shifted_vel) then
       plot_names(counter) = "shifted_velx"
       counter = counter + 1
       plot_names(counter) = "shifted_vely"
       counter = counter + 1
       if (dm > 2) then
          plot_names(counter) = "shifted_velz"
          counter = counter + 1
       end if
    end if

    if (plot_umac_tavg) then 
       plot_names(counter) = "tavg_averaged_velx"
       counter = counter + 1
       plot_names(counter) = "tavg_averaged_vely"
       counter = counter + 1
       if (dm > 2) then
          plot_names(counter) = "tavg_averaged_velz"
          counter = counter + 1
       end if

       if (plot_shifted_vel) then
          plot_names(counter) = "tavg_shifted_velx"
          counter = counter + 1
          plot_names(counter) = "tavg_shifted_vely"
          counter = counter + 1
          if (dm > 2) then
             plot_names(counter) = "tavg_shifted_velz"
             counter = counter + 1
          end if
       end if
    end if

    plot_names(counter) = "pres"
    counter = counter + 1

    if (use_charged_fluid) then
       plot_names(counter) = "charge_density"
       counter = counter + 1

       plot_names(counter) = "Epot"
       counter = counter + 1

       if (plot_Epot_tavg) then
          plot_names(counter) = "tavg_Epot" 
          counter = counter + 1
       end if

       if (plot_gradEpot) then
          plot_names(counter) = "averaged_Ex"
          counter = counter + 1
          plot_names(counter) = "averaged_Ey"
          counter = counter + 1
          if (dm > 2) then
             plot_names(counter) = "averaged_Ez"
             counter = counter + 1
          end if
       end if

       if (plot_avg_gradPhiApprox) then
          plot_names(counter) = "av_gradPhiApproxx"
          counter = counter + 1
          plot_names(counter) = "av_gradPhiApproxy"
          counter = counter + 1
          if (dm > 2) then
             plot_names(counter) = "av_gradPhiApproxz"
             counter = counter + 1
          end if
       end if

    end if

    if(boussinesq) then
       plot_names(counter) = "rho_eos"
       counter = counter+1
    end if

    ! staggered quantities

    ! vel
    nvarsStag = 1
    if (use_charged_fluid) then
       ! electric field (Ex, Ey, Ez)
       if (plot_gradEpot) then
          nvarsStag = nvarsStag + 1
       end if
       ! gradPhiApprox
       if (plot_avg_gradPhiApprox) then 
          nvarsStag = nvarsStag + 1
       end if
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
       if (plot_gradEpot) then
          plot_names_stagx(counter) = "Ex"
          plot_names_stagy(counter) = "Ey"
          plot_names_stagz(counter) = "Ez"
          counter = counter + 1
       end if

       if (plot_avg_gradPhiApprox) then 
          plot_names_stagx(counter) = "gradPhiApproxx"
          plot_names_stagy(counter) = "gradPhiApproxy"
          plot_names_stagz(counter) = "gradPhiApproxz"
          counter = counter + 1
       end if
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

    if (nspecies .eq. 4) then
       ! concentrations
       call convert_rhoc_to_c(mla,rho,rhotot,cc_temp,.true.)
       
       ! w_Na = m_Na / (m_Na+m_OH) * w_NaOH + m_Na / (m_Na+m_Cl) * w_NaCl
       ! w_Cl = m_Cl / (m_H+m_Cl) * w_HCl + m_Cl / (m_Na+m_Cl) * w_NaCl
       ! w_H = m_H / (m_H+m_Cl) * w_HCl
       ! w_OH = m_OH / (m_Na+m_OH) * w_NaOH

       ! Na
       do n=1,nlevs
          call multifab_setval_c(plotdata(n),0.d0,counter,1,all=.true.)
          ! from NaOH
          massfrac = m_Na / (m_Na+m_OH)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),2,1)
          ! from NaCl
          massfrac = m_Na / (m_Na+m_Cl)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),3,1)
       end do
       counter = counter+1

       ! Cl
       do n=1,nlevs
          call multifab_setval_c(plotdata(n),0.d0,counter,1,all=.true.)
          ! from HCl
          massfrac = m_Cl / (m_H+m_Cl)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),1,1)
          ! from NaCl
          massfrac = m_Cl / (m_Na+m_Cl)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),3,1)
       end do
       counter = counter+1

       ! H
       do n=1,nlevs
          call multifab_setval_c(plotdata(n),0.d0,counter,1,all=.true.)
          ! from HCl
          massfrac = m_H / (m_H+m_Cl)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),1,1)
       end do
       counter = counter+1

       ! OH
       do n=1,nlevs
          call multifab_setval_c(plotdata(n),0.d0,counter,1,all=.true.)
          ! from NaOH
          massfrac = m_OH / (m_Na+m_OH)
          call multifab_saxpy_3_cc(plotdata(n),counter,massfrac,cc_temp(n),2,1)
       end do
       counter = counter+1

       ! H2O
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),counter,cc_temp(n),4,1,0)
       end do
       counter = counter+1

       if (boussinesq) then
          ! we only want rho's.  multiply c's by rhotot
          do n=1,nlevs
             do i=1,5
                call multifab_mult_mult_c(plotdata(n),counter-i,rhotot(n),1,1,0)
             end do
          end do
       else
          ! we want rho's and c's.
          do n=1,nlevs
             ! copy c's into next indices
             call multifab_copy_c(plotdata(n),counter,plotdata(n),counter-5,5,0)
             ! multiply first set of c's by rhotot
             do i=1,5
                call multifab_mult_mult_c(plotdata(n),counter-i,rhotot(n),1,1,0)
             end do
             counter = counter+5
          end do
       end if

    else

       ! rho
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),counter,rho(n),1,nspecies_eff,0)
       end do
       counter = counter + nspecies_eff

       if(.not.boussinesq) then  
          ! compute concentrations
          call convert_rhoc_to_c(mla,rho,rhotot,cc_temp,.true.)
          do n=1,nlevs
             call multifab_copy_c(plotdata(n),counter,cc_temp(n),1,nspecies_eff,0)
          end do
          counter = counter + nspecies_eff
       end if

    end if

    ! Temp
    if(is_nonisothermal) then
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),counter,Temp(n),1,1,0)
       enddo
       counter = counter + 1
    end if   

    ! vel averaged
    do i=1,dm
       call average_face_to_cc(mla,umac(:,i),1,plotdata,counter,1)
       counter = counter + 1
    end do

    if(plot_shifted_vel) then
       ! vel shifted
       do i=1,dm
          call shift_face_to_cc(mla,umac(:,i),1,plotdata,counter,1)
          counter = counter + 1
       end do
    end if   

    ! time-averaged vel averaged and time-averaged vel shifted
    if (plot_umac_tavg) then
       do i=1,dm
          call average_face_to_cc(mla,umac_avg(:,i),1,plotdata,counter,1)
          counter = counter + 1
       end do
       if (plot_shifted_vel) then
          do i=1,dm
             call shift_face_to_cc(mla,umac_avg(:,i),1,plotdata,counter,1)
             counter = counter + 1
          end do
       end if
    end if

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

       ! time-averaged Epot 
       if (plot_Epot_tavg) then 
          do n = 1,nlevs
             call multifab_copy_c(plotdata(n),counter,Epot_avg(n),1,1,0)
          enddo
          counter = counter + 1
       end if

       ! averaged electric field
       if (plot_gradEpot) then
          do i=1,dm
             call average_face_to_cc(mla,grad_Epot(:,i),1,plotdata,counter,1)
             do n = 1,nlevs
                call multifab_mult_mult_s_c(plotdata(n),counter,-1.d0,1,0)
             end do
             counter = counter + 1
          end do
       end if

       ! averaged gradPhiApprox
       if (plot_avg_gradPhiApprox) then 
          do i=1,dm
             call average_face_to_cc(mla,gradPhiApprox(:,i),1,plotdata,counter,1)
             counter = counter + 1
          end do
       end if 
    end if

    if(boussinesq) then
       ! rho_eos - rho0
       call compute_rhotot_eos(mla,rho,plotdata,counter)
       do n=1,nlevs
          call multifab_sub_sub_s_c(plotdata(n),counter,rho0,1)
       end do
       counter = counter+1
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

       ! electric field
       if (plot_gradEpot) then
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(plotdata_stag(n,i),counter,grad_Epot(n,i),1,1)
                call multifab_mult_mult_s_c(plotdata_stag(n,i),counter,-1.d0,1,0)
             end do
          end do
          counter = counter + 1
       end if

       ! gradPhiApprox
       if (plot_avg_gradPhiApprox) then 
          do n=1,nlevs
             do i=1,dm
                call multifab_copy_c(plotdata_stag(n,i),counter,gradPhiApprox(n,i),1,1)
             end do
          end do
          counter = counter + 1
       end if

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
