module compute_rhoh_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use eos_model_wrapper_module
  use convert_stag_module
  use div_and_grad_module
  use bc_module
  use probin_common_module, only: nspecies

  implicit none

  private

  public :: compute_rhoh_fluxdiv

contains

  subroutine compute_rhoh_fluxdiv(mla,lambda,Temp,diff_mass_flux,rhotot,diff_rhoh_fluxdiv, &
                                  dx,time,the_bc_tower)
       
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: lambda(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: diff_mass_flux(:,:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: diff_rhoh_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n,nlevs,i,dm,comp

    type(multifab) :: hk(mla%nlevel)
    type(multifab) :: misc_fc(mla%nlevel,mla%dim) ! used for various purposes
    type(multifab) :: gradT(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(hk(n),mla%la(n),nspecies,1)
       do i=1,dm
          call multifab_build_edge(misc_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(  gradT(n,i),mla%la(n),       1,0,i)
       end do
    end do

    ! average lambda to faces and store in misc_fc
    call average_cc_to_face(nlevs,lambda,misc_fc,1,tran_bc_comp,1, &
                            the_bc_tower%bc_tower_array)

    ! compute grad(T)
    call compute_grad(mla,Temp,gradT,dx,1,temp_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! set misc_fc = lambda*grad(T)
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(misc_fc(n,i),1,gradT(n,i),1,1,0)
       end do
    end do

    ! set diff_rhoh_fluxdiv = div(lambda grad T)
    call compute_div(mla,misc_fc,diff_rhoh_fluxdiv,dx,1,1,1)

    ! compute a multifab holding h_k
    call compute_hk(mla,Temp,hk)

    ! average all components h_k to faces and store in misc_fc
    call average_cc_to_face(nlevs,hk,misc_fc,1,c_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array)

    ! set (misc_fc)_k = -h_k F_k = -h_k (-rho*W*chi*Gamma*grad(x) - ... )
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(misc_fc(n,i),1,diff_mass_flux(n,i),1,nspecies,0)
       end do
    end do

    ! add sum (div (-h_k F_k)) to diff_rhoh_fluxdiv
    do comp=1,nspecies

       ! add divergence
       call compute_div(mla,misc_fc,diff_rhoh_fluxdiv,dx,comp,1,1,increment_in=.true.)

    end do

    ! add external heating
    call add_external_heating(mla,rhotot,diff_rhoh_fluxdiv,dx,time)

    do n=1,nlevs
       call multifab_destroy(hk(n))
       do i=1,dm
          call multifab_destroy(misc_fc(n,i))
          call multifab_destroy(gradT(n,i))
       end do
    end do

  end subroutine compute_rhoh_fluxdiv
  
end module compute_rhoh_fluxdiv_module
