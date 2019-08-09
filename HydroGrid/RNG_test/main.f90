program main

  use BoxLib
  use parallel
  use BoxLibRNGs
  use ParallelRNGs
  use bl_random_module
  use bl_error_module

  implicit none

  integer, parameter :: dp=kind(0.0d0), sp=kind(0.0)

  ! input parameters (read through command line by read_params)

  integer        :: rng_eng_type = 1               ! 1=built-in HydroGrid RNG engine
                                                   ! 2=c++ RNG engine
  integer        :: seed = 0                       ! random seed
  integer        :: dist_type = 1                  ! 1=uniform
                                                   ! 2=standard normal
                                                   ! 3=Poisson (dist_param1=mean)
                                                   ! 4=binomial (dist_param1=n_trial,dist_param2=sussec_prob)
  integer        :: prec = 1                       ! 1=single precision (sp)
                                                   ! 2=double precision (dp)
  real(sp)       :: param_mean_sp = 0.0            ! mean of Poisson (sp)
  real(dp)       :: param_mean_dp = 0.0d0          ! mean of Poisson (dp)
  integer        :: param_n_trials = 0             ! n_trials of binomial
  real(sp)       :: param_success_prob_sp = 0.0d0  ! success probability of binomial (sp)
  real(dp)       :: param_success_prob_dp = 0.0d0  ! success probability of binomial (dp)
  integer        :: gen_randnum_int = 0            ! (gen_randnum) how many random numbers to be generated
  character(128) :: gen_randnum_file = ''          ! (gen_randnum) output filename
  integer        :: hist_n_samples = 0             ! (hist) how many random numbers to be used for histogram
  integer        :: hist_n_bins = 0                ! (hist) how many bins for histogram
  real(dp)       :: hist_left_end = 0.0d0          ! (hist) value for the left-most bin 
  real(dp)       :: hist_right_end = 0.0d0         ! (hist) value for the right-most bin 
  logical        :: hist_is_center_val = .false.   ! (hist) is hist_left/right_end a center-value?
  character(128) :: hist_file = ''                 ! (hist) output filename

  ! variables

  type(bl_rng_engine) :: rng_eng

  !!!!!!!!
  ! init !
  !!!!!!!!

  call boxlib_initialize()

  if (parallel_nprocs() .ne. 1) then
    call bl_error('this program assumes a single-processor run')
  end if

  call read_params

  select case (rng_eng_type)
  case (1) ! HydroGrid
    call SeedParallelRNG(seed)
  case (2) ! c++
    call bl_rng_build_engine(rng_eng,seed)
  end select

  !!!!!!!!!!!!!!!
  ! gen_randnum !
  !!!!!!!!!!!!!!!

  if (gen_randnum_int .gt. 0) then
    call gen_randnum
  end if

  !!!!!!!!
  ! hist !
  !!!!!!!!
 
  if (hist_n_samples .gt. 0) then
    call hist
  end if

  !!!!!!!!!!!!!!
  ! terminiate !
  !!!!!!!!!!!!!!

  if (rng_eng_type .eq. 2) then  ! c++
    call bl_rng_destroy_engine(rng_eng)
  end if

  call boxlib_finalize()

  contains

    subroutine get_one_random_number(res,res_sp,res_int)

      real(dp), intent(out)           :: res
      real(sp), intent(out), optional :: res_sp
      integer,  intent(out), optional :: res_int 

      real(sp) :: randnum_sp
      real(dp) :: randnum_dp
      integer  :: randnum_int

      select case (dist_type)
      case (1) ! uniform

        select case (prec)
        case (1) ! sp

          select case (rng_eng_type)
          case (1) ! HyroGrid
            call UniformRNG_sp(randnum_sp) 
          case (2) ! c++
            call UniformRNG_sp(randnum_sp,rng_eng%p)
          end select

          ! check
          if (randnum_sp .eq. 1.) then
            call bl_error('UniformRNG_sp returns exact 1.')
          end if

          if (present(res_sp)) then
            res_sp = randnum_sp  ! no datatype conversion
          end if

          res = randnum_sp  ! conversion: sp -> dp 

        case (2) ! dp

          select case (rng_eng_type)
          case (1) ! HyroGrid
            call UniformRNG_dp(randnum_dp) 
          case (2) ! c++
            call UniformRNG_dp(randnum_dp,rng_eng%p)
          end select

          ! check
          if (randnum_dp .eq. 1.d0) then
            call bl_error('UniformRNG_dp returns exact 1.d0')
          end if

          res = randnum_dp  ! no datatype conversion

        end select 

      case (2) ! standard normal

        select case (prec)
        case (1) ! sp

          select case (rng_eng_type)
          case (1) ! HydroGrid
            call NormalRNG_sp(randnum_sp)
          case (2) ! c++
            call NormalRNG_sp(randnum_sp,rng_eng%p)
          end select
         
          if (present(res_sp)) then
            res_sp = randnum_sp  ! no datatype conversion
          end if

          res = randnum_sp  ! conversion: sp -> dp
        
        case (2) ! dp

          select case (rng_eng_type)
          case (1) ! Hydrogrid
            call NormalRNG_dp(randnum_dp)
          case (2) ! c++
            call NormalRNG_dp(randnum_dp,rng_eng%p)
          end select

          res = randnum_dp  ! no datatype conversion

        end select

      case (3) ! Poisson

        select case (prec)
        case (1) ! sp

          select case (rng_eng_type)
          case (1) ! HydroGrid
            call PoissonRNG_sp(randnum_int,param_mean_sp)
          case (2) ! c++
            call PoissonRNG_sp(randnum_int,param_mean_sp,rng_eng%p)
          end select

          ! check
          if (randnum_int .lt. 0) then
            print *,randnum_int
            call bl_error('PoissonRNG_sp returns a negative number')
          end if

        case (2) ! dp

          select case (rng_eng_type)
          case (1) ! HydroGrid
            call PoissonRNG_dp(randnum_int,param_mean_dp)
          case (2) ! c++
            call PoissonRNG_dp(randnum_int,param_mean_dp,rng_eng%p)
          end select

          ! check
          if (randnum_int .lt. 0) then
            print *,randnum_int
            call bl_error('PoissonRNG_dp returns a negative number')
          end if

        end select

        if (present(res_int)) then
          res_int = randnum_int  ! no datatype conversion
        end if

        res = randnum_int ! conversion: int -> dp

      case (4) ! binomial

        select case (prec)
        case (1) ! sp

          select case (rng_eng_type)
          case (1) ! HydroGrid
            call BinomialRNG_sp(randnum_int,param_n_trials,param_success_prob_sp)
          case (2) ! c++
            call BinomialRNG_sp(randnum_int,param_n_trials,param_success_prob_sp,rng_eng%p)
          end select

          ! check
          if ((randnum_int .lt. 0) .or. (randnum_int .gt. param_n_trials)) then
            print *,randnum_int
            call bl_error('BinomialRNG_sp returns a strange number')
          end if

        case (2) ! dp

          select case (rng_eng_type)
          case (1) ! HydroGrid
            call BinomialRNG_dp(randnum_int,param_n_trials,param_success_prob_dp)
          case (2) ! c++
            call BinomialRNG_dp(randnum_int,param_n_trials,param_success_prob_dp,rng_eng%p)
          end select

          ! check
          if ((randnum_int .lt. 0) .or. (randnum_int .gt. param_n_trials)) then
            print *,randnum_int
            call bl_error('BinomialRNG_dp returns a strange number')
          end if

        end select

        if (present(res_int)) then
          res_int = randnum_int  ! no datatype conversion
        end if
 
        res = randnum_int ! conversion: int -> dp

      end select

    end subroutine get_one_random_number

    subroutine gen_randnum

      integer  :: funit = 24
      integer  :: fstat

      real(dp) :: randnum_dp
      real(sp) :: randnum_sp
      integer  :: randnum_int

      integer  :: i

      ! file open
      open (unit=funit,file=gen_randnum_file,action='write',iostat=fstat)
      if (fstat .ne. 0) then
        call bl_error('gen_randnum_file not be opened')
      end if

      ! generate random numbers
      do i=1,gen_randnum_int

        if ((dist_type .eq. 1) .or. (dist_type .eq. 2)) then  ! real RNG

          select case (prec)
          case (1) ! sp
            call get_one_random_number(randnum_dp,res_sp=randnum_sp)
            write (funit,*) randnum_sp
          case (2) ! dp
            call get_one_random_number(randnum_dp)
            write (funit,*) randnum_dp
          end select

        else if ((dist_type .eq. 3) .or. (dist_type .eq. 4)) then ! integer RNG
     
          call get_one_random_number(randnum_dp,res_int=randnum_int)
          write (funit,*) randnum_int

        end if

      end do

      ! file close
      close (unit=funit,iostat=fstat)
      if (fstat .ne. 0) then
        call bl_error('while attempting to close gen_randnum_file')
      end if

    end subroutine gen_randnum
    
    subroutine hist 

      real(dp) :: dx
      real(dp) :: a
      real(dp) :: b
      integer  :: bin(hist_n_bins)
      integer  :: cnt_left = 0    ! count of samples .lt. a
      integer  :: cnt_right = 0   ! count of samples .ge. b

      real(dp) :: sum1 = 0.d0
      real(dp) :: sum2 = 0.d0

      real(dp) :: randnum_dp
      integer  :: n               ! for bin
      integer  :: i               ! for sample

      integer  :: funit = 25
      integer  :: fstat

      character(32) :: tmp

      ! calc dx, a, b
      if (hist_is_center_val) then
        dx = (hist_right_end-hist_left_end)/(hist_n_bins-1)
        a = hist_left_end-0.5d0*dx
        b = hist_right_end+0.5d0*dx
      else
        dx = (hist_right_end-hist_left_end)/hist_n_bins
        a = hist_left_end
        b = hist_right_end
      end if

      ! init bin
      do n=1,hist_n_bins
        bin(n)=0
      end do

      ! main loop
      do i=1,hist_n_samples
        call get_one_random_number(randnum_dp)

        n = floor((randnum_dp-a)/dx)+1
        if (n .lt. 1) then
          cnt_left = cnt_left+1
        else if (n .gt. hist_n_bins) then
          cnt_right = cnt_right+1
        else 
          bin(n) = bin(n)+1
        end if

        sum1 = sum1+randnum_dp
        sum2 = sum2+randnum_dp*randnum_dp
      end do 

      sum1 = sum1/hist_n_samples
      sum2 = sum2/hist_n_samples
      write (tmp,*) sum1
      print *,'** sample mean= ',adjustl(tmp)
      write (tmp,*) sum2-sum1*sum1
      print *,'** sample variance= ',adjustl(tmp)

      print *,''  ! new line
      write (tmp,*) cnt_left
      print *,'** hist: cnt_left= ',adjustl(tmp)
      write (tmp,*) cnt_right
      print *,'** hist: cnt_right= ',adjustl(tmp)

      ! file open
      open (unit=funit,file=hist_file,action='write',iostat=fstat)
      if (fstat .ne. 0) then
        call bl_error('hist_file not be opened')
      end if

      do n=1,hist_n_bins
        write (funit,*) a+(n-0.5)*dx,bin(n)/dx/hist_n_samples
      end do

      ! file close
      close (unit=funit,iostat=fstat)
      if (fstat .ne. 0) then
        call bl_error('while attempting to close hist_file')
      end if

    endsubroutine hist

    subroutine read_params

      integer            :: narg,farg
      character(len=128) :: fname
      character(len=3)   :: prec_tag

      !!!!!!!!
      ! read !
      !!!!!!!!

      narg = command_argument_count()
      farg = 1

      do while (farg<=narg)

        call get_command_argument(farg,value=fname)

        select case (fname)

        case ('--rng_eng_type')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) rng_eng_type

        case ('--seed')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) seed 

        case ('--dist_type')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) dist_type 

        case ('--prec')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) prec 

        case ('--param_mean')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) param_mean_sp 
          call get_command_argument(farg,value=fname)
          read(fname, *) param_mean_dp 

        case ('--param_n_trials')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) param_n_trials 

        case ('--param_success_prob')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) param_success_prob_sp
          call get_command_argument(farg,value=fname)
          read(fname, *) param_success_prob_dp

        case ('--gen_randnum_int')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) gen_randnum_int

        case ('--gen_randnum_file')
          farg = farg+1
          call get_command_argument(farg,value=gen_randnum_file)

        case ('--hist_n_samples')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) hist_n_samples
      
        case ('--hist_n_bins')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) hist_n_bins
  
        case ('--hist_left_end')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) hist_left_end

        case ('--hist_right_end')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) hist_right_end

        case ('--hist_is_center_val')
          farg = farg+1
          call get_command_argument(farg,value=fname)
          read(fname, *) hist_is_center_val

        case ('--hist_file')
          farg = farg+1
          call get_command_argument(farg,value=hist_file)

        case default
          print *,fname
          call bl_error('invalid command-line option')

        end select

        farg = farg+1 
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! check values and report !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! rng_eng_type
      print *,''  ! new line
      select case (rng_eng_type)
      case (1)
        print *,'** built-in HydroGrid RNG engine'
      case (2)
        print *,'** c++ RNG engine'
      case default
        print *,rng_eng_type
        call bl_error('invalid rng_eng_type')
      end select

      ! seed
      write (fname,*) seed
      print *,'** seed= ',adjustl(fname)

      ! dist_type and prec
      select case (dist_type)
      case (1) ! uniform
        fname = 'UniformRNG'
      case (2) ! standard normal
        fname = 'NormalRNG'
      case (3) ! Posson
        fname = 'PoissonRNG'
      case (4) ! binomial
        fname = 'BinomialRNG'
      case default
        print *,dist_type 
        call bl_error('invalid dist_type') 
      end select

      select case (prec)
      case (1)
        prec_tag = '_sp'
      case (2)
        prec_tag = '_dp'
      case default
        print *,prec
        call bl_error('invalid prec')
      end select

      fname = trim(fname) // prec_tag
      print *,'** ',fname

      ! distribution parameters
      select case (dist_type)

      case (1) ! uniform: no params

      case (2) ! standard normal: no params

      case (3) ! Posson

        select case (prec)
        case (1) ! sp
          write (fname,*) param_mean_sp
        case (2) ! dp
          write (fname,*) param_mean_dp
        end select
       
        print *, '** mean= ',adjustl(fname)

      case (4) ! binomial

        write (fname,*) param_n_trials 
        print *, '** n_trials= ',adjustl(fname)

        select case (prec)
        case (1) ! sp
          write (fname,*) param_success_prob_sp
        case (2) ! dp
          write (fname,*) param_success_prob_dp
        end select

        print *, '** success_prob= ',adjustl(fname)

      case default

        print *,dist_type
        call bl_error('invalid dist_type') 

      end select

      ! gen_randnum_int and gen_randnum_file
      if (gen_randnum_int .gt. 0) then
        print *, ''  ! new line

        ! gen_randnum_int
        write (fname,*) gen_randnum_int
        print *, '** gen_randnum_int= ',adjustl(fname)

        ! gen_randnum_file
        if (len(trim(gen_randnum_file)) .eq. 0) then
          call bl_error('set gen_randnum_int=0 or provide gen_randnum_file')
        else
          print *, '** gen_randnum_file= ',trim(gen_randnum_file)
        end if
      end if

      ! hist
      if (hist_n_samples .gt. 0) then
        print *, ''  ! new line

        ! hist_n_samples
        write (fname,*) hist_n_samples
        print *, '** hist_n_samples= ',adjustl(fname)

        ! hist_n_bins
        if (hist_n_bins .gt. 0) then
          write (fname,*) hist_n_bins
          print *, '** hist_n_bins= ',adjustl(fname)
        else
          print *,hist_n_bins
          call bl_error('hist_n_bins should be positive')
        end if

        ! hist_left_end and hist_right_end
        if (hist_left_end .lt. hist_right_end) then
          write (fname,*) hist_left_end
          print *, '** hist_left_end= ',adjustl(fname)
          write (fname,*) hist_right_end
          print *, '** hist_right_end= ',adjustl(fname)
        else
          print *,hist_left_end,hist_right_end
          call bl_error('hist_left_end should be less than hist_right_end')
        end if

        ! hist_is_center_val
        if (hist_is_center_val) then
          print *, '** hist_is_center_val= T'
        else
          print *, '** hist_is_center_val= F'
        end if

        ! hist_file
        if (len(trim(hist_file)) .eq. 0) then
          call bl_error('hist_file is missing')
        else
          print *, '** hist_file= ',trim(hist_file)
        end if
      end if

      ! final new line
      print *,''

    end subroutine read_params

end program main
