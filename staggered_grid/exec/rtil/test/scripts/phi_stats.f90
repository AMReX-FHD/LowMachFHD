program phi_stats

  character(len=64) fname
  double precision xval(3),temp
  double precision xphi(64,3)
  double precision avg(3), var(3), std(3)
  
  integer i, nFiles

  read(5,*) nFiles
  
  do i=1,nFiles
  
     read(5,*) fname

     open(1,file=fname)
     
     read(1,*) xval(1), xphi(i,1), temp
     read(1,*) xval(2), xphi(i,2), temp
     read(1,*) xval(3), xphi(i,3), temp
     
     close(1)

  end do

  ! compute average at each height
  avg(:) = 0.
  do i=1,nFiles
     avg(1) = avg(1) + xphi(i,1)
     avg(2) = avg(2) + xphi(i,2)
     avg(3) = avg(3) + xphi(i,3)
  end do
  avg(:) = avg(:) / nFiles

  ! compute variance at each height
  var(:) = 0.
  do i=1,nFiles
     var(1) = var(1) + (xphi(i,1)-avg(1))**2
     var(2) = var(2) + (xphi(i,2)-avg(2))**2
     var(3) = var(3) + (xphi(i,3)-avg(3))**2
  end do
  var(:) = var(:) / nFiles

  std(:) = sqrt(var(:))

  
  print*,xval(1),avg(1),var(1),std(1)
  print*,xval(2),avg(2),var(2),std(2)
  print*,xval(3),avg(3),var(3),std(3)

end program phi_stats
