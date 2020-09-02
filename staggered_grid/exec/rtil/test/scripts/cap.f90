program capacitance
  
  double precision phi(100), cap(100), dcap(100), errorbar(100)
  double precision dlogx, sigma, h1,h2,a,b,c

  integer nFiles

  integer n

  read(5,*) nFiles
  
  do n=1,nFiles
     read(5,*)phi(n),cap(n),errorbar(n)
  enddo

  do n=1,nFiles-1
     dcap(n)     = (cap(n+1)-cap(n))/(phi(n+1)-phi(n))
     errorbar(n) = (errorbar(n+1)+errorbar(n))/(phi(n+1)-phi(n))
  enddo

  do n=1,nFiles-1
     write(6,*)0.5d0*(phi(n)+phi(n+1)),dcap(n),errorbar(n)
  enddo

end program
