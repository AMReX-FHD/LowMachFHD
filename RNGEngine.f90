
module hg_rng_engine_module
  use iso_c_binding
  implicit none

  type, bind(c) :: hg_rng_engine
     type(c_ptr) :: eng = c_null_ptr   ! engine
     type(c_ptr) :: dis = c_null_ptr   ! uniform distribution [0,1)
  end type hg_rng_engine
end module hg_rng_engine_module
