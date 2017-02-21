EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../../inputs_SU_2d

RUNNAME=TEST

#OPT1="--max_step 10000 --plot_int 100 --print_int 10 --hydro_grid_int 0 --n_steps_skip 0 --fixed_dt 1.e-7"
#OPT1="--max_step 10000 --plot_int 0 --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --fixed_dt 1.e-7"
OPT1="--max_step 100000 --plot_int 0 --print_int 1000 --hydro_grid_int 1 --n_steps_skip 10000 --fixed_dt 1.e-7"
#OPT1="--max_step 1000000 --plot_int 0 --print_int 1000 --hydro_grid_int 1 --n_steps_skip 100000 --fixed_dt 1.e-7"

#OPT2="--temporal_integrator 1 --diffusion_type 3 --reaction_type 0 --use_Poisson_rng 2"
OPT2="--temporal_integrator -2 --use_Poisson_rng 1"

#####

if [ ! -f $EXEC ]
then
  echo "ERROR: $EXEC does not exist"
  exit
fi

RUNNAME=RUN_$RUNNAME
if [ -d $RUNNAME ]
then
  echo "ERROR: $RUNNAME exists"
  exit
else
  mkdir $RUNNAME
  cp $0 $RUNNAME
  cp $INPUTS $RUNNAME
  cd $RUNNAME
fi

OPTS="$OPT1 $OPT2"

echo "mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out"
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

#####

grep n_avg scr_out | awk '{print $4, $5}' > res.n_avg

if [ -f SU.S_k.pair=1.Re.dat ]
then
  sed -i 's/0.00000000/#0.00000000/g' SU.S_k.pair=*.Re.dat
fi
