program main

  use BoxLib
  use BoxLibRNGs
  use bl_random_module
  use ParallelRNGs

  implicit none

  integer, parameter :: dp=kind(0.0d0), sp=kind(0.0)

  type(bl_rng_engine) :: rng_eng

  real(dp) :: real_dp
  real(sp) :: real_sp

  integer :: seed
  integer :: i

  call boxlib_initialize()

  seed = 1

  call SeedParallelRNG(seed)

  call bl_rng_build_engine(rng_eng,seed)

  do i=1,10
     call UniformRNG(real_dp,rng_eng%p)
     print*,'BL UniformRNG with dp',real_dp
  end do

  do i=1,10
     call UniformRNG(real_dp)
     print*,'HG UniformRNG with dp',real_dp
  end do

  do i=1,10
     call UniformRNG(real_sp,rng_eng%p)
     print*,'BL UniformRNG with sp',real_sp
  end do

  do i=1,10
     call UniformRNG(real_sp)
     print*,'HG UniformRNG with sp',real_sp
  end do

  call bl_rng_destroy_engine(rng_eng)

  call boxlib_finalize()

end program main
