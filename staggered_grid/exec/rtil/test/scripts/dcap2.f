           double precision cap(100),phi(100), dcap(100)
           double precision dlogx, sigma, h1,h2,a,b,c


           integer n


           do n=1,8
           read(5,*)phi(n),cap(n)
           enddo


           do n=1,7

           dcap(n) = (cap(n+1)-cap(n))/(phi(n+1)-phi(n))

           enddo

           do n=1,7
           write(6,*)0.5d0*(phi(n)+phi(n+1)),dcap(n)
           enddo

           end
