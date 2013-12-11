
!----------------------------------------------------------------------

          subroutine D_GIOVANGIGLI(nspec,Dbin,molecular_weight,Xk,Diff_ij)

! nspec is number of species, Dbin is matrix of binary diffusion coefficient
!  molecular_weight is vector of molecular weights (molecular masses will work as well)
!  Xk is mole fractions

           implicit none 
 
          integer :: nspec, i, j, k, jj

          real(kind=8) :: Xk(1:nspec)
          real(kind=8) :: molecular_weight(1:nspec)
          real(kind=8) :: Dbin(1:nspec,1:nspec), term1, term2
          real(kind=8) :: Di(1:nspec), Diff_ij(1:nspec,1:nspec)
          real(kind=8) :: Deltamat(1:nspec,1:nspec), Zmat(1:nspec,1:nspec)
          real(kind=8), dimension(1:nspec,1:nspec) :: Pmat, Jmat
          real(kind=8), dimension(1:nspec) :: Minv, Mmat
          real(kind=8), dimension(1:nspec,1:nspec) :: PJ, matrix1, matrix2
          real(kind=8) :: scr
          real(kind=8), dimension(1:nspec) :: xxtr, yytr
          real(kind=8) :: Ykp(1:nspec), Xkp(1:nspec)

          integer :: jmax


          jmax = 3 

          ! mole fractions correction - EGLIB
          do ii = 1, nspec
           Xkp(ii) = Xk(ii) + (1.0d-16)*(sum(Xk(:))/dble(nspec)-Xk(ii))
          enddo

          ! molecular weight of mixture - EGLIB
          Mwmix = 0.0d0
         do ii = 1, nspec
          MWmix = MWmix + Xkp(ii)*molecular_weight(ii)
         enddo
          wbar = MWmix

          ! mass fractions correction - EGLIB
         do ii = 1, nspec
          Ykp(ii) = molecular_weight(ii)/MWmix*Xkp(ii)
         enddo


          ! Find Di matrix 
          do i = 1, nspec
           term2 = 0.0d0
           do j = 1, nspec
            if(j.ne.i) then
              term2 = term2 + Xkp(j)/Dbin(i,j)
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
                term1 = term1 + Xkp(i)*Xkp(k)/Dbin(i,k)
               endif
              enddo  
              Deltamat(i,i) = term1
             else
              Deltamat(i,j) = -Xkp(i)*Xkp(j)/Dbin(i,j) 
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
             Diff_ij(i,j) = scr
           enddo
          enddo

          if(jmax.gt.0)then

          do jj = 1,jmax
          

!         matrix1=0
          do i = 1, nspec
           do j = 1, nspec
            scr = 0.d0
            do k = 1, nspec
               scr = scr + PJ(i,k)*Diff_ij(k,j)
            enddo
             matrix1(i,j) = scr+matrix2(i,j)
           enddo
          enddo 

          Diff_ij=matrix1

          enddo

          endif

        return
        end subroutine


