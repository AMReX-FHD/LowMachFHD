module convert_variables_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: convert_m_to_umac, convert_cons_to_prim
  
contains

  subroutine convert_m_to_umac(mla,s_fc,m,umac,m_to_umac)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) ::    s_fc(:,:)
    type(multifab) , intent(inout) ::    m(:,:)  
    type(multifab) , intent(inout) :: umac(:,:)
    logical        , intent(in   ) :: m_to_umac

    ! local
    integer :: n,i,dm,nlevs

    dm = mla%dim
    nlevs = mla%nlevel

    if (m_to_umac) then

       ! compute umac = m / rho - NO GHOST CELLS
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac(n,i), 1, m(n,i), 1, 1, 0)
             call multifab_div_div_c(umac(n,i), 1, s_fc(n,i), 1, 1, 0)
          end do
       end do

    else

       ! compute m = rho * umac - INCLUDING GHOST CELLS
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(m(n,i), 1, umac(n,i), 1, 1, m(n,i)%ng)
             call multifab_mult_mult_c(m(n,i), 1, s_fc(n,i), 1, 1, m(n,i)%ng)
          end do
       end do

    end if

  end subroutine convert_m_to_umac

  subroutine convert_cons_to_prim(mla,s,prim,cons_to_prim)
    
    use probin_lowmach_module, only: nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) ::    s(:)
    type(multifab) , intent(inout) :: prim(:)
    logical        , intent(in   ) :: cons_to_prim

    ! local
    integer :: nlevs,n,i

    nlevs = mla%nlevel

    if (cons_to_prim) then

       ! cons to prim - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(prim(n),1,s(n),1,nscal,0)
          do i=2,nscal
             call multifab_div_div_c(prim(n),i,s(n),1,1,0)
          end do
       end do

    else

       ! prim to cons - INCLUDING GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(s(n),1,prim(n),1,nscal,s(n)%ng)
          do i=2,nscal
             call multifab_mult_mult_c(s(n),i,s(n),1,1,s(n)%ng)
          end do
       end do

    end if
    
  end subroutine convert_cons_to_prim

end module convert_variables_module
