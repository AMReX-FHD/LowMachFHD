program dphidn
  
  double precision x(3),phi(3),errorbar(3),eb_tot
  double precision dx, sigma

  integer n

  eb_tot = 0.

  do n=1,3
     read(5,*) x(n),phi(n),errorbar(n)
  enddo

  dx = x(3)-x(2)

  sigma = -9.2e-20*(9.d0*phi(2)-phi(3)-8.d0*phi(1))/(3.d0*dx)

  eb_tot = 9.2e-20*(9.d0*errorbar(2)+errorbar(3))/(3.d0*dx)
  
  write(6,*) phi(1),sigma,eb_tot

  ! throwaway
  do n=1,188     
     read(5,*) x(1),phi(1),errorbar(1)
  end do

  do n=1,3
     read(5,*) x(n),phi(n),errorbar(n)
  enddo

  sigma = 9.2e-20*(9.d0*phi(2)-phi(1)-8.d0*phi(3))/(3.d0*dx)

  eb_tot = 9.2e-20*(9.d0*errorbar(2)+errorbar(1))/(3.d0*dx)
  
  write(6,*) phi(3),sigma,eb_tot

end program
