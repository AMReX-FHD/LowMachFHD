module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use define_bc_module
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: dim_in, restart
  use probin_reactdiff_module, only: nspecies, use_bl_rng

  implicit none

  private

  public :: initialize_from_restart

contains

  subroutine initialize_from_restart(mla,time,dt,n_old,pmask,ng_s)
 
     type(ml_layout),intent(  out) :: mla
     real(dp_t)    , intent(  out) :: time,dt
     type(multifab), intent(inout) :: n_old(:)
     logical       , intent(in   ) :: pmask(:)
     integer       , intent(in   ) :: ng_s

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(layout)              :: la

     integer :: n,nlevs

     type(bl_prof_timer), save :: bpt

     call build(bpt,"initialize_from_restart")

     call fill_restart_data(mba,chkdata,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     do n = 1,nlevs
        call multifab_build(n_old(n), mla%la(n), nspecies, ng_s)
     end do
     do n = 1,nlevs
        call multifab_copy_c(n_old(n), 1, chkdata(n), 1, nspecies)
        !
        ! The layout for chkdata is built standalone, level
        ! by level, and need to be destroy()d as such as well.
        !
        la = get_layout(chkdata(n))
        call multifab_destroy(chkdata(n))
        call destroy(la)
     end do
     deallocate(chkdata)
     call destroy(mba)

     call destroy(bpt)

  end subroutine initialize_from_restart

  subroutine fill_restart_data(mba,chkdata,time,dt)


    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    character(len=11)                 :: sd_name
    character(len=40)                 :: rand_name
    integer                           :: n,nlevs
    integer                           :: rrs(10)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_restart_data")

    write(unit=sd_name,fmt='("chk",i8.8)') restart

    if ( parallel_IOProcessor() ) then
       print *,'Reading ',sd_name,' to get state data for restart'
    end if

    call checkpoint_read(chkdata, sd_name, rrs, time, dt, nlevs)

    call build(mba,nlevs,dim_in)
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
       rand_name = sd_name//'/rng_binomial_diffusion'
       call bl_rng_restore(rng_binomial_diffusion, rand_name)

       rand_name = sd_name//'/rng_normal_diffusion'
       call bl_rng_restore(rng_normal_diffusion, rand_name)

       rand_name = sd_name//'/rng_poisson_reaction'
       call bl_rng_restore(rng_poisson_reaction, rand_name)

       rand_name = sd_name//'/rng_normal_reaction'
       call bl_rng_restore(rng_normal_reaction, rand_name)

       rand_name = sd_name//'/rng_uniform_real_reaction'
       call bl_rng_restore(rng_uniform_real_reaction, rand_name)
    end if

    call destroy(bpt)

  end subroutine fill_restart_data

end module restart_module
