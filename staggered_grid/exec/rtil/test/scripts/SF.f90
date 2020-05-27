! read in a vtk file for 192^2 dataset
! modify the vtk file to only include the coordinates and data for a given pair
! compute the expected value with a filter eliminating small signals
program SF

  double precision x_face(193)
  double precision x(192)
  double precision data(192,192)

  double precision row(9)
  double precision sum1,sum2,val,peak

  integer counter,a,b

  counter = 1
  do b=1,21
     read(*,*) row(1),row(2),row(3),row(4),row(5),row(6),row(7),row(8),row(9)
     do a=1,9
        x_face(counter) = row(a)
        counter = counter+1
     end do
  end do

  read(*,*) row(1),row(2),row(3),row(4)
  do a=1,4
     x_face(counter) = row(a)
     counter = counter+1
  end do

  do a=1,192
     x(a) = 0.5*(x_face(a)+x_face(a+1))
  end do

!  do a=1,193
!     print*,'hack',a,x_face(a)
!  end do

  peak = 0.
  counter = 0
  do b=1,4096
  
     read(*,*) row(1),row(2),row(3),row(4),row(5),row(6),row(7),row(8),row(9)

     do a=1,9
        i = mod(counter,192)+1
        j = counter/192+1
        data(i,j) = row(a)
        if (data(i,j) .gt. peak) then
           peak = data(i,j)
        end if
!        if (data(i,j) .lt. 1.d-24) then
!           data(i,j) = data(i,j)*0
!        end if
        counter = counter+1
     end do

  end do
  
  sum1 = 0.
  sum2 = 0.
  counter = 0
  do a=1,36864
     i = mod(counter,192)+1
     j = counter/192+1
     if (x(i) .ne. 0 .and. x(j) .ne. 0 .and. data(i,j) .gt. 0.01d0*peak) then
        sum1 = sum1 + sqrt(x(i)**2 + x(j)**2) * data(i,j)
        sum2 = sum2 + data(i,j)
     end if
     counter = counter+1
  end do
  print*,'k_r =',sum1/sum2
  
!  sum1 = 0.
!  sum2 = 0.
!  counter = 0
!  do a=1,36864
!     i = mod(counter,192)+1
!     j = counter/192+1
!     if (x(i) .ne. 0 .and. x(j) .ne. 0 .and. data(i,j) .gt. 0.02d0*peak) then
!        sum1 = sum1 + data(i,j)
!        sum2 = sum2 + data(i,j)/sqrt(x(i)**2+x(j)**2)
!     end if
!     counter = counter+1
!  end do
!  print*,'k_r =',sum1/sum2
  
end program SF
