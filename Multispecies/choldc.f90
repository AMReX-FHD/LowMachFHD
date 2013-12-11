! this fragment zeros rows where mass fractions are really small
!  make cholesky clean

             do ns=1,NCOMP
                if(abs(yy(ns)) + abs(yyp(ns)) .le. 1.d-12)then
                   DijY_edge(ns,1:NCOMP)=0.d0
                   DijY_edge(1:NCOMP,ns)=0.d0
                endif
             enddo


!  a is input matrix.  sqda is cholesky factor
!  p is used for diagonal 
!  other things in list are diagnostics and testing

           subroutine choldc(a,np,sqda)

       integer :: np

       real(kind=8) :: a(np,np), p(np), sqda(np,np), sqda2(np,np), dij(np,np)
       real(kind=8) :: dd(np,np)
       real(kind=8) :: yy(np), mwmix 

       integer :: i, j, k, ii, jj

       real(kind=8) :: sum1

       real(kind=8) :: small_number



       integer :: idiag,ising



        small_number = 0.d0



        idiag = 0


       do i = 1, np

           ising = 0

        do j = i, np

           sum1 = a(i,j)

           do k = i-1, 1, -1

              sum1 = sum1 - a(i,k)*a(j,k)

           enddo

           if(i.eq.j) then

             if(sum1.le.small_number) then

             p(i) = 0.d0

             ising = 1

             else

             p(i) = sqrt(sum1)

             endif

           else

             if(ising.eq.0)then

                a(j,i) = sum1/p(i)

             else

                a(j,i) = 0.d0

             endif

           endif

        enddo

       enddo

        
         sqda = 0.0d0

         do i = 1, np

          do j = i-1, 1, -1

           sqdA(i,j) = a(i,j)

          enddo

          sqdA(i,i) = p(i)

         enddo

           

         if(idiag.eq.1)then



              print*, 'a: '

              do ii = 1, np

               write(6,1000) (a(ii,jj),jj=1,np)
1000           format(9e24.16)

              enddo
            
              call flush()


              print*, 'sqdA: '

              do ii = 1, np

               write(6,1000) (sqdA(ii,jj),jj=1,np)

              enddo

              call flush()


              do j=1,np

              do i=1,np

                  sum = 0.d0

                  do k=1,np

                     sum = sum + sqdA(i,k)*sqdA(j,k)

                  enddo

                  a(i,j) = sum

              enddo

              enddo

              print*, 'sqdA **2'

              do ii = 1, np

               write(6,1000) (a(ii,jj),jj=1,np)

              enddo

              call flush()


              print*, 'dtilde '

              do ii = 1, np

               write(6,1000) (dij(ii,jj),jj=1,np)

              enddo

              call flush()

              print*, 'yy, mwmix: ', yy(1:np), mwmix

              call flush()


              ! compute d
              do ii = 1, np
               do jj = 1, np 
                 if(yy(ii).gt.0.0d0) then
                  dd(ii,jj) = dij(ii,jj)/yy(ii)
                 else
                  dd(ii,jj) = 0.0d0 
                 endif   
               enddo            
              enddo  

              print*, 'd_ij: '
              do ii = 1, np

               write(6,1000) (dd(ii,jj),jj=1,np)

              enddo

              call flush() 



         endif


        return

        end subroutine

