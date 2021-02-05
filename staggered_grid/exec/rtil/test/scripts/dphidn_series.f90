program dphidn
  
  double precision x(200),phi(200),errorbar(200),eb_tot
  double precision dx, sigma, sigma2

  integer istart,iend,iinc
  double precision dt, toffset, time

  integer n,nf
  
  character*16 filename
  character*7 step

  write(6,*)"istart,iend,iinc,dt,toffset"
  read(5,*)istart,iend,iinc,dt, toffset
  open(2,file="dphidn_v_t",form='formatted')

  do nf = istart,iend,iinc

  write(step,'(i7.7)') nf
  filename = "Epot_" // step // ".dat"
  write(6,*)filename
  open(1,file=filename,form='formatted')
  rewind 1

  do n=1,98
     read(1,*) x(n),phi(n),errorbar(n)
  enddo
  close(1)

  dx = x(3)-x(2)

  sigma = -9.2e-20*(9.d0*phi(2)-phi(3)-8.d0*phi(1))/(3.d0*dx)
  sigma2 = 9.2e-20*(9.d0*phi(97)-phi(96)-8.d0*phi(98))/(3.d0*dx)

  time = dfloat(nf)*dt -toffset
  write(2,*)time,sigma,sigma2

  enddo

end program
