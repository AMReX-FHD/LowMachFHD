module multifab_fill_random_module

  use multifab_module
  use BoxLibRNGs
  use bl_rng_module
  use bl_random_module
  use bl_error_module
  use probin_reactdiff_module, only: use_bl_rng

  implicit none

  private

  public :: multifab_fill_random

contains

  ! fill a multifab with random numbers
  subroutine multifab_fill_random(mfab, comp, variance, variance_mfab, variance_mfab2, &
                                  init_in)
    type(multifab), intent(inout)           :: mfab(:)
    integer       , intent(in   ), optional :: comp ! Only one component
    real(dp_t)    , intent(in   ), optional :: variance
    type(multifab), intent(in   ), optional :: variance_mfab(:)
    type(multifab), intent(in   ), optional :: variance_mfab2(:)
    logical       , intent(in   ), optional :: init_in

    integer :: n,box

    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    logical :: init

    init = .false.
    if (present(init_in)) init = init_in

    !--------------------------------------
    do n=1,size(mfab)
       do box = 1, nfabs(mfab(n))
          if(present(comp)) then
             fp => dataptr(mfab(n),box,comp,1)
          else
             fp => dataptr(mfab(n),box)
          end if

          ! Fill the whole grid with random numbers
          if (use_bl_rng) then
             if (init) then
                call NormalRNGs(fp, size(fp), rng_eng_init%p)
             else
                call NormalRNGs(fp, size(fp), rng_eng_diffusion%p)
             end if
          else
             call NormalRNGs(fp, size(fp))
          end if

          if(present(variance_mfab)) then 
             ! Must have same distribution
             if (mfab(n)%ng .ne. variance_mfab(n)%ng) then
                call bl_error('ERROR: multifab_fill_random; mfab and variance_mfab require same # of ghost cells')
             end if
             fpvar => dataptr(variance_mfab(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if

          if(present(variance_mfab2)) then 
             ! Must have same distribution
             if (mfab(n)%ng .ne. variance_mfab2(n)%ng) then
                call bl_error('ERROR: multifab_fill_random; mfab and variance_mfab2 require same # of ghost cells')
             end if
             fpvar => dataptr(variance_mfab2(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if

          if(present(variance)) then
             fp=sqrt(variance)*fp
          end if

       end do
    end do

  end subroutine multifab_fill_random

end module multifab_fill_random_module
