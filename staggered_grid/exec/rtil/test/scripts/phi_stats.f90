program phi_stats

  character(len=64) fname
  double precision xval(194),temp
  double precision xphi(64,194)
  double precision avg(194), var(194), errorbar(194)
  
  integer i,j,nFiles

  read(5,*) nFiles
  
  do i=1,nFiles
  
     read(5,*) fname

     open(1,file=fname)

     do j=1,194
        read(1,*) xval(j), xphi(i,j), temp
     end do
     
     close(1)

  end do

  ! compute average at each height
  avg(:) = 0.
  do i=1,nFiles
     do j=1,194
        avg(j) = avg(j) + xphi(i,j)
     end do
  end do
  avg(:) = avg(:) / nFiles

  ! compute variance at each height
  var(:) = 0.
  do i=1,nFiles
     do j=1,194
        var(j) = var(j) + (xphi(i,j)-avg(j))**2
     end do
  end do
  var(:) = var(:) / nFiles

  errorbar(:) = sqrt(var(:)/nFiles)

  do j=1,194
     print*,xval(j),avg(j),errorbar(j)
  end do

end program phi_stats
