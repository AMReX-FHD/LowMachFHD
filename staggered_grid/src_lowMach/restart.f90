module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use define_bc_module
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: dim_in, advection_type, restart, advection_type, &
                                  use_bl_rng, nspecies, &
                                  seed_momentum, seed_diffusion, seed_reaction, &
                                  check_base_name
  use probin_chemistry_module, only: nreactions
  use probin_charged_module, only: use_charged_fluid

  implicit none

  private

  public :: initialize_from_restart

contains

  subroutine initialize_from_restart(mla,time,dt,rho,rho_sum,rhotot,pi,umac,umac_sum,pmask,Epot,Epot_sum,grad_Epot)
 
     type(ml_layout),intent(out)   :: mla
     real(dp_t)    , intent(  out) :: time,dt
     type(multifab), intent(inout) :: rho(:)
     type(multifab), intent(inout) :: rho_sum(:)  
     type(multifab), intent(inout) :: rhotot(:)
     type(multifab), intent(inout) :: pi(:)
     type(multifab), intent(inout) :: Epot(:)
     type(multifab), intent(inout) :: Epot_sum(:) 
     type(multifab), intent(inout) :: umac(:,:)
     type(multifab), intent(inout) :: umac_sum(:,:)
     type(multifab), intent(inout) :: grad_Epot(:,:)
     logical       , intent(in   ) :: pmask(:)

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(multifab), pointer   :: chkdata_edgex(:)
     type(multifab), pointer   :: chkdata_edgey(:)
     type(multifab), pointer   :: chkdata_edgez(:)
     type(layout)              :: la

     integer :: n,nlevs,i,dm,ng_s

     type(bl_prof_timer), save :: bpt

     call build(bpt,"initialize_from_restart")

     dm = dim_in

     if (advection_type .eq. 0) then
        ng_s = 2 ! centered advection
     else if (advection_type .le. 3) then
        ng_s = 3 ! bilinear bds or unlimited quadratic bds
     else if (advection_type .eq. 4) then
        ng_s = 4 ! limited quadratic bds
     end if

     call fill_restart_data(mba,chkdata,chkdata_edgex,chkdata_edgey,chkdata_edgez,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     do n = 1,nlevs
        call multifab_build(rho(n)   , mla%la(n), nspecies, ng_s)
        call multifab_build(rho_sum(n), mla%la(n), nspecies, ng_s) 
        call multifab_build(rhotot(n), mla%la(n),        1, ng_s)
        call multifab_build(pi(n)    , mla%la(n),        1, 1)
        if (use_charged_fluid) then
           call multifab_build(Epot(n)    , mla%la(n),        1, 1)  
           call multifab_build(Epot_sum(n)    , mla%la(n),        1, 1)  
        endif
        do i=1,dm
           call multifab_build_edge(umac(n,i), mla%la(n), 1, 1, i)
           call multifab_build_edge(umac_sum(n,i), mla%la(n), 1, 1, i)

           ! with mixed boundary conditions some of the corner umac ghost cells that
           ! never affect the solution aren't filled, causing segfaults on some compilers
           ! this prevents segfaults
           call setval(umac(n,i),0.d0,all=.true.)
           call setval(umac_sum(n,i),0.d0,all=.true.)

           if (use_charged_fluid) then
              call multifab_build_edge(grad_Epot(n,i), mla%la(n), 1, 1, i)
              call setval(grad_Epot(n,i),0.d0,all=.true.)                     ! for ghost cells as well
           endif

        end do
     end do

     ! cell-centered data
     do n = 1,nlevs
        call multifab_copy_c(rho(n)   , 1,chkdata(n) ,1           ,nspecies) 
        call multifab_copy_c(rho_sum(n), 1,chkdata(n),nspecies+1 ,nspecies)
        call multifab_copy_c(rhotot(n), 1,chkdata(n) ,2*nspecies+1  ,1)
        call multifab_copy_c(pi(n)    , 1,chkdata(n) ,2*nspecies+2  ,1)
        if (use_charged_fluid) then
           call multifab_copy_c(Epot(n)    , 1,chkdata(n) ,2*nspecies+3  ,1)
           call multifab_copy_c(Epot_sum(n)    , 1,chkdata(n) ,2*nspecies+4  ,1)
        endif
     end do

     ! edge data
     do n=1,nlevs
        call multifab_copy_c(umac(n,1),1,chkdata_edgex(n),1,1)
        call multifab_copy_c(umac(n,2),1,chkdata_edgey(n),1,1)
        call multifab_copy_c(umac_sum(n,1),1,chkdata_edgex(n),2,1)      
        call multifab_copy_c(umac_sum(n,2),1,chkdata_edgey(n),2,1)
        if (use_charged_fluid) then
           call multifab_copy_c(grad_Epot(n,1),1,chkdata_edgex(n),3,1)  
           call multifab_copy_c(grad_Epot(n,2),1,chkdata_edgey(n),3,1)
        endif
        if (dm .eq. 3) then
           call multifab_copy_c(umac(n,3),1,chkdata_edgez(n),1,1)
           call multifab_copy_c(umac_sum(n,3),1,chkdata_edgez(n),2,1)   
           if (use_charged_fluid) then
              call multifab_copy_c(grad_Epot(n,3),1,chkdata_edgez(n),3,1)  
           endif
        end if
     end do

     do n=1,nlevs
        !
        ! The layout for chkdata is built standalone, level
        ! by level, and need to be destroy()d as such as well.
        !
        la = get_layout(chkdata(n))
        call multifab_destroy(chkdata(n))
        call destroy(la)
        la = get_layout(chkdata_edgex(n))
        call multifab_destroy(chkdata_edgex(n))
        call destroy(la)
        la = get_layout(chkdata_edgey(n))
        call multifab_destroy(chkdata_edgey(n))
        call destroy(la)
        if (dm .eq. 3) then
           la = get_layout(chkdata_edgez(n))
           call multifab_destroy(chkdata_edgez(n))
           call destroy(la)
        end if
     end do

     deallocate(chkdata)
     deallocate(chkdata_edgex)
     deallocate(chkdata_edgey)
     if (dm .eq. 3) then
        deallocate(chkdata_edgez)
     end if
     call destroy(mba)

     call destroy(bpt)

  end subroutine initialize_from_restart

  subroutine fill_restart_data(mba,chkdata, &
                               chkdata_edgex,chkdata_edgey,chkdata_edgez, &
                               time,dt)


    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chkdata_edgex(:)
    type(multifab)   , pointer        :: chkdata_edgey(:)
    type(multifab)   , pointer        :: chkdata_edgez(:)
    character(len=8)                  :: check_index
    character(len=128)                :: sd_name
    character(len=128)                :: rand_name
    integer                           :: n,nlevs,dm
    integer                           :: rrs(10)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_restart_data")

    dm = dim_in

    write(unit=check_index,fmt='(i8.8)') restart
    sd_name = trim(check_base_name) // check_index

    if ( parallel_IOProcessor() ) then
       print *,'Reading ',sd_name,' to get state data for restart'
    end if

    call checkpoint_read(chkdata, chkdata_edgex, chkdata_edgey, chkdata_edgez, &
                         sd_name, rrs, time, dt, nlevs)

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = rrs(n-1)
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

     ! random state
    if (use_bl_rng) then

       if (seed_momentum .eq. -1) then
          rand_name = trim(sd_name) // '/rng_eng_mom'
          call bl_rng_restore_engine(rng_eng_momentum, rand_name)
       end if

       if (seed_diffusion .eq. -1) then
          rand_name = trim(sd_name) // '/rng_eng_diff'
          call bl_rng_restore_engine(rng_eng_diffusion, rand_name)
       end if

       if (seed_reaction .eq. -1) then
          rand_name = trim(sd_name) // '/rng_eng_react'
          call bl_rng_restore_engine(rng_eng_reaction, rand_name)
       end if

    end if

    call destroy(bpt)

  end subroutine fill_restart_data

end module restart_module
