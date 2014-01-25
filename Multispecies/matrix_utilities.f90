module matrix_utilities

  use bl_types
  use probin_multispecies_module

  implicit none

  private

  public :: Dbar2chi_iterative
  
  ! To avoid zero mass or mole fractions we round up using a tolerance:
  real(kind=dp_t), public :: fraction_tolerance = 2*epsilon(1.0_dp_t)

contains

    ! nspec is number of species
    ! num_iterations is the number of terms in the sum to use: 3-5 are reasonable values
    ! Dbar is matrix of Maxwell-Stefan binary diffusion coefficient
    ! chi is the multispecies diffusion matrix
    ! Wk is vector of molecular weights (molecular masses will work as well)
    ! Xk is mole fractions
    subroutine Dbar2chi_iterative(nspec,num_iterations,Dbar,Wk,Xk,chi)
      integer, intent(in) :: nspec, num_iterations
      real(kind=dp_t), intent(in) :: Xk(1:nspec), Wk(1:nspec), Dbar(1:nspec,1:nspec)
      real(kind=dp_t), intent(out) :: chi(1:nspec,1:nspec)
      
      ! Local variables
      real(kind=dp_t) :: term1, term2, MWmix
      real(kind=dp_t) :: Di(1:nspec)
      real(kind=dp_t) :: Deltamat(1:nspec,1:nspec), Zmat(1:nspec,1:nspec)
      real(kind=dp_t), dimension(1:nspec,1:nspec) :: Pmat, Jmat
      real(kind=dp_t), dimension(1:nspec) :: Minv, Mmat
      real(kind=dp_t), dimension(1:nspec,1:nspec) :: PJ, matrix1, matrix2
      real(kind=dp_t) :: scr
      real(kind=dp_t) :: Ykp(1:nspec), Xkp(1:nspec)

      integer :: i, j, k, ii, jj

      ! mole fractions correction - EGLIB
      do ii = 1, nspec
       Xkp(ii) = Xk(ii) + fraction_tolerance*(sum(Xk(:))/dble(nspec)-Xk(ii))
      enddo

      ! molecular weight of mixture - EGLIB
      Mwmix = 0.0d0
      do ii = 1, nspec
       MWmix = MWmix + Xkp(ii)*Wk(ii)
      enddo

      ! mass fractions correction - EGLIB
      do ii = 1, nspec
       Ykp(ii) = Wk(ii)/MWmix*Xkp(ii)
      enddo


      ! Find Di matrix 
      do i = 1, nspec
       term2 = 0.0d0
       do j = 1, nspec
        if(j.ne.i) then
          term2 = term2 + Xkp(j)/Dbar(i,j)
        endif
       enddo   
       Di(i) = (1.d0-Ykp(i))/term2 
      enddo   

      ! Compute Mmat and Minv
      do i = 1, nspec
       Mmat(i) = Xkp(i)/Di(i)
       Minv(i) = Di(i)/Xkp(i)
      enddo

      ! Compute P matrix
      Pmat = 0.0d0
      do i = 1, nspec
       do j = 1, nspec
         Pmat(i,j) = - Ykp(j) 
         if(i.eq.j) then
          Pmat(i,j) =  Pmat(i,j) + 1.0d0  
         endif
       enddo
      enddo

      ! Compute Deltamat
      Deltamat = 0.0d0 
      do i = 1, nspec
       do j = 1, nspec
         if(i.eq.j) then
          term1 = 0.0d0
          do k = 1, nspec
           if(k.ne.i) then
            term1 = term1 + Xkp(i)*Xkp(k)/Dbar(i,k)
           endif
          enddo  
          Deltamat(i,i) = term1
         else
          Deltamat(i,j) = -Xkp(i)*Xkp(j)/Dbar(i,j) 
         endif  
          Zmat(i,j) = -Deltamat(i,j)
       enddo
      enddo  

      ! Compute Zmat
      do i = 1, nspec
        Zmat(i,i) = Zmat(i,i) + Mmat(i)
      enddo  

      ! Compute Jmat
      do i = 1, nspec
       do j = 1, nspec
         Jmat(i,j) = Minv(i)*Zmat(i,j)
        enddo
       enddo

      ! Compute PJ
      PJ = 0.0d0
      do i = 1, nspec
       do j = 1, nspec
        do k = 1, nspec
         PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
        enddo
       enddo
      enddo

      ! Compute P M^-1 Pt; store it in matrix2
      do i = 1, nspec
       do j = 1, nspec
        scr = 0.d0
        do k = 1, nspec
         scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
            ! notice the change in indices for Pmat to represent Pmat^t
        enddo
         matrix2(i,j) = scr
         chi(i,j) = scr
       enddo
      enddo

      do jj = 1,num_iterations
       do i = 1, nspec
        do j = 1, nspec
         scr = 0.d0
         do k = 1, nspec
            scr = scr + PJ(i,k)*chi(k,j)
         enddo
          matrix1(i,j) = scr+matrix2(i,j)
        enddo
       enddo 
       chi=matrix1
      enddo

  end subroutine

end module
