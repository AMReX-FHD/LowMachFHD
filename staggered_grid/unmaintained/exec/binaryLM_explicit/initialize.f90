module initialize_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use init_module
  use box_util_module
  use make_new_grids_module

  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, &
            initialize_with_adaptive_grids, make_new_state, initialize_bc

contains

  subroutine initialize_from_restart(mla,restart,time,dt,dx,sold,mold,pold,the_bc_tower)
 
     use probin_module, only : dim_in, nlevs, nscal, ng_scal, ng_mom, pmask

     type(ml_layout),intent(out)   :: mla
     integer       , intent(in   ) :: restart
     real(dp_t)    , intent(  out) :: time,dt
     real(dp_t)    , pointer       :: dx(:,:)
     type(multifab), pointer       :: sold(:)
     type(multifab), pointer       :: mold(:,:)
     type(multifab), pointer       :: pold(:)
     type(bc_tower), intent(  out) :: the_bc_tower

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(multifab), pointer   :: chkdata_edgex(:)
     type(multifab), pointer   :: chkdata_edgey(:)
     type(multifab), pointer   :: chkdata_edgez(:)
     type(layout)              :: la

     integer :: n,dm,i

     dm = dim_in

     call fill_restart_data(restart,mba,chkdata,chkdata_edgex, &
                            chkdata_edgey,chkdata_edgez,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     allocate(mold(nlevs,dm),sold(nlevs))

     do n = 1,nlevs
        do i=1,dm
           call multifab_build_edge(mold(n,i), mla%la(n), 1, ng_mom, i)
        end do
        call multifab_build(sold(n), mla%la(n), nscal, ng_scal)
     end do
     do n = 1,nlevs

        call multifab_copy_c(sold(n),1,chkdata(n),1,nscal)
        call multifab_copy_c(mold(n,1),1,chkdata_edgex(n),1,1)
        call multifab_copy_c(mold(n,2),1,chkdata_edgey(n),1,1)
        if (dm .eq. 3) then
           call multifab_copy_c(mold(n,3),1,chkdata_edgez(n),1,1)
        end if
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

     call initialize_dx(dx,mba,nlevs)

     call initialize_bc(the_bc_tower,nlevs,dm)
     do n = 1,nlevs
        call bc_tower_level_build(the_bc_tower,n,mla%la(n))
     end do

     call destroy(mba)

  end subroutine initialize_from_restart

  subroutine initialize_with_fixed_grids(mla,dx,mold,sold,the_bc_tower)

     use probin_module, only : dim_in, nlevs, nscal, ng_scal, ng_mom, fixed_grids, pmask

     type(ml_layout),intent(out   ):: mla
     real(dp_t)    , pointer       :: dx(:,:)
     type(multifab), pointer       :: mold(:,:)
     type(multifab), pointer       :: sold(:)
     type(bc_tower), intent(  out) :: the_bc_tower

     type(ml_boxarray)         :: mba

     integer :: n,dm,i

     dm = dim_in

     call read_a_hgproj_grid(mba, fixed_grids)
     call ml_layout_build(mla,mba,pmask)

     ! check for proper nesting
     if (.not. ml_boxarray_properly_nested(mla%mba, ng_scal, pmask)) &
         call bl_error('fixed_grids not properly nested')

     nlevs = mla%nlevel
     allocate(mold(nlevs,dm),sold(nlevs))

     do n = 1,nlevs
        do i=1,dm
           call multifab_build_edge(mold(n,i), mla%la(n), 1, ng_mom, i)
        end do
        call multifab_build(sold(n), mla%la(n), nscal, ng_scal)
     end do

     call initialize_dx(dx,mba,nlevs)

     call initialize_bc(the_bc_tower,nlevs,dm)

     do n = 1,nlevs
        call bc_tower_level_build( the_bc_tower,n,mla%la(n))
     end do

     call initdata(mold,sold,dx,mla,0.d0)

     call destroy(mba)

  end subroutine initialize_with_fixed_grids

  subroutine initialize_with_adaptive_grids(mla,dx,mold,sold,the_bc_tower)

     use probin_module, only : dim_in, nlevs, nscal, n_cells, &
                               regrid_int, amr_buf_width, max_grid_size, &
                               ref_ratio, max_levs, verbose, pmask

     type(ml_layout),intent(inout)  :: mla
     real(dp_t)    , pointer        :: dx(:,:)
     type(multifab), pointer        :: mold(:,:)
     type(multifab), pointer        :: sold(:)
     type(bc_tower), intent(  out)  :: the_bc_tower

     type(layout)                   :: la_array(max_levs)
     type(box)                      :: bxs
     type(ml_boxarray)              :: mba

     logical :: new_grid
     integer :: lo(dim_in), hi(dim_in)
     integer :: n, nl, dm, i

     dm = dim_in

     ! set up hi & lo to carry indexing info
     lo = 0
     hi(1) = n_cells(1)-1
     if (dm > 1) then   
        hi(2) = n_cells(2) - 1        
        if (dm > 2)  then
           hi(3) = n_cells(3) -1
        endif
     endif

     ! mba is big enough to hold max_levs levels
     call ml_boxarray_build_n(mba,max_levs,dm)
     do n = 1, max_levs-1
        mba%rr(n,:) = ref_ratio
     enddo

     allocate(mold(max_levs,dm),sold(max_levs))

       ! Build the level 1 boxarray
     call box_build_2(bxs,lo,hi)
     call boxarray_build_bx(mba%bas(1),bxs)
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     if ( parallel_IOprocessor() ) then
        call print(mba%bas(1), 'Level 1 BoxArray')
        print*,''
     end if

     ! build pd(:)
     mba%pd(1) = bxs
     do n = 2, max_levs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     enddo

     ! Need to build pd before making dx
     call initialize_dx(dx,mba,max_levs)

     ! Initialize bc's.
     call initialize_bc(the_bc_tower,max_levs,dm)

     ! Build the level 1 layout.
     call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

     ! Build the level 1 data only.
     call make_new_state(la_array(1),mold(1,:),sold(1))

     ! Define bc_tower at level 1.
     call bc_tower_level_build(the_bc_tower,1,la_array(1))

     nlevs = 1

     if (max_levs > 1) then

        ! Initialize the level 1 data only.
        call initdata_on_level(mold(1,:),sold(1),dx(1,:),0.d0)

        new_grid = .true.
        nl = 1

        do while ( (nl .lt. max_levs) .and. (new_grid) )

           ! Do we need finer grids?
           ! AJN - fixme for max_grid_size(dm)
!           call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
!                               amr_buf_width,ref_ratio,nl,max_grid_size)
        
           if (new_grid) then

              call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

              ! Build the level nl+1 data only.
              call make_new_state(la_array(nl+1),mold(nl+1,:),sold(nl+1))

              ! Define bc_tower at level nl+1.
              call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
            
             ! fills the physical region of each level with problem data (blob now)
              call initdata_on_level(mold(nl+1,:),sold(nl+1),dx(nl+1,:),0.d0)

              nlevs = nl+1
              nl = nl + 1

           endif ! if (new_grid) 

      enddo          

      do n = 1,nlevs
         do i=1,dm
            call destroy(mold(n,i))
         end do
         call destroy(sold(n))
      end do

      nlevs = nl

      if (nlevs .ge. 3) then

         ! check for proper nesting
         ! AJN - fixme for max_grid_size(dm)
!         call enforce_proper_nesting(mba,la_array,max_grid_size)

         ! enforce_proper_nesting can create new grids at coarser levels
         ! this makes sure the boundary conditions are properly defined everywhere
         do n = 2,nlevs
            call bc_tower_level_build(the_bc_tower,n,la_array(n))
         end do

      end if

   else
      
      do i=1,dm
         call destroy(mold(1,i))
      end do
      call destroy(sold(1))

   end if ! end if (maxlev > 1)

   call ml_layout_restricted_build(mla,mba,nlevs,pmask)

   nlevs = mla%nlevel

   do n = 1, nlevs
      call destroy(la_array(n))
   end do

   do n = 1,nlevs
      call make_new_state(mla%la(n),mold(n,:),sold(n))
   end do

   call initdata(mold,sold,dx,mla,0.d0)

   call destroy(mba)

  end subroutine initialize_with_adaptive_grids

  subroutine make_new_state(la_loc,mold_loc,sold_loc)

    use probin_module, only : dim_in, nscal, ng_scal, ng_mom
    use bl_constants_module

    type(layout)  , intent(in   ) :: la_loc
    type(multifab), intent(inout) :: mold_loc(:)
    type(multifab), intent(inout) :: sold_loc
 
    integer :: dm,i

    dm = dim_in

    do i=1,dm
       call multifab_build_edge(mold_loc(i), la_loc, 1, ng_mom, i)
       call setval(mold_loc(i), ZERO, all=.true.)
    end do
    call multifab_build( sold_loc, la_loc, nscal, ng_scal)
    call setval(sold_loc,ZERO, all=.true.)

  end subroutine make_new_state

  subroutine initialize_bc(the_bc_tower,num_levs,dm)

     use bc_module

     use probin_module, only : bc_lo, bc_hi, pmask

     type(bc_tower), intent(  out) :: the_bc_tower
     integer       , intent(in   ) :: num_levs,dm

     integer :: domain_phys_bc(dm,2)

     ! Define the physical boundary conditions on the domain
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bc_lo(1)
     domain_phys_bc(1,2) = bc_hi(1)
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bc_lo(1) .ne. -1 .or. bc_hi(1) .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bc_lo(2)
        domain_phys_bc(2,2) = bc_hi(2)
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bc_lo(2) .ne. -1 .or. bc_hi(2) .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bc_lo(3)
        domain_phys_bc(3,2) = bc_hi(3)
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bc_lo(3) .ne. -1 .or. bc_hi(3) .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.

     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

  end subroutine initialize_bc

  subroutine initialize_dx(dx,mba,num_levs)

     use probin_module, only : dim_in, prob_lo, prob_hi
  
     real(dp_t)       , pointer     :: dx(:,:)
     type(ml_boxarray), intent(in ) :: mba
     integer          , intent(in ) :: num_levs

     integer :: i,n,dm

     dm = dim_in

     allocate(dx(num_levs,dm))

     do i = 1,dm
        dx(1,i) = (prob_hi(i)-prob_lo(i)) / float(mba%pd(1)%hi(i)-mba%pd(1)%lo(i)+1)
     end do
     do n = 2,num_levs
        dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
     end do
     
     do n=1,num_levs
        do i=2,dm
           if (dx(n,1) .ne. dx(n,i)) then
              call bl_error('initialize.f90: dx must be the same in each direction')
           end if
        enddo
     enddo

  end subroutine initialize_dx

end module initialize_module
