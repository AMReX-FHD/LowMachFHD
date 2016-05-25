program coarsen
   implicit none
   
   integer, parameter :: nx=256, ny=256
   integer, parameter :: ratio(2)=4
   integer, parameter :: nxc=nx/ratio(1), nyc=ny/ratio(2)
   
   integer :: i,j
   real :: original(nx,ny), coarsened(nxc,nyc)
   
   read(*,*) original
   
   do j=1, nyc
   do i=1, nxc
      coarsened(i,j)=sum(original( (i-1)*ratio(1)+1:(i-1)*ratio(1)+ratio(1), &
                                   (j-1)*ratio(2)+1:(j-1)*ratio(2)+ratio(2) )) / product(ratio)
      write(*,*) coarsened(i,j)                                   
   end do
   end do

end program
