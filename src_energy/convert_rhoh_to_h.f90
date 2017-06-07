module convert_rhoh_to_h_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module

  implicit none

  private

  public :: convert_rhoh_to_h, fill_rhoh_ghost_cells, &
            fill_h_ghost_cells, fill_Temp_ghost_cells
  
contains

  subroutine convert_rhoh_to_h(mla,rhoh,rhotot,h,rhoh_to_h)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rhoh(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: h(:)
    logical        , intent(in   ) :: rhoh_to_h

    ! local
    integer :: n,nlevs

    nlevs = mla%nlevel

    if (rhoh_to_h) then

       ! rhoh to h - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(h(n),1,rhoh(n),1,1,0)
          call multifab_div_div_c(h(n),1,rhotot(n),1,1,0)
       end do

    else

       ! h to rhoh- VALID + GHOST (CAN CHANGE TO DO ONLY GHOST TO SAVE COMPUTATION)
       do n=1,nlevs
          call multifab_copy_c(rhoh(n),1,h(n),1,1,rhoh(n)%ng)
          call multifab_mult_mult_c(rhoh(n),1,rhotot(n),1,1,rhoh(n)%ng)
       end do

    end if

  end subroutine convert_rhoh_to_h

  ! ghost cells for rhotot must be filled
  subroutine fill_rhoh_ghost_cells(mla,rhoh,rhotot,dx,the_bc_tower)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rhoh(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    type(multifab) :: enth(mla%nlevel)

    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(enth(n),mla%la(n),1,rhoh(n)%ng)
    end do

    ! compute h from rhoh in VALID REGION
    call convert_rhoh_to_h(mla,rhoh,rhotot,enth,.true.)

    ! fill h ghost cells
    call fill_h_ghost_cells(mla,enth,dx,the_bc_tower)

    ! compute h from rhoh - INCLUDING GHOST CELLS
    call convert_rhoh_to_h(mla,rhoh,rhotot,enth,.false.)

    do n=1,nlevs
       call multifab_destroy(enth(n))
    end do

  end subroutine fill_rhoh_ghost_cells

  subroutine fill_h_ghost_cells(mla,enth,dx,the_bc_tower)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: enth(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    nlevs = mla%nlevel

    do n=1,nlevs
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(enth(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(enth(n),1,h_bc_comp,1, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

  end subroutine fill_h_ghost_cells

  subroutine fill_Temp_ghost_cells(mla,Temp,dx,the_bc_tower)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: Temp(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    nlevs = mla%nlevel

    do n=1,nlevs
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(Temp(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(Temp(n),1,temp_bc_comp,1, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

  end subroutine fill_Temp_ghost_cells

end module convert_rhoh_to_h_module
