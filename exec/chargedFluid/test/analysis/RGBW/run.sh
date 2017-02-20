EXEC=../../main.Linux.gfortran.mpi.exe
INPUTS=../../inputs_RGBW_react

RUNNAME=plot_a0.01_5T

# OPT1: chose one
OPT1="--algorithm_type 5 --include_reactions T"
#OPT1="--algorithm_type 5 --include_reactions F"
#OPT1="--algorithm_type 0 --include_reactions F"

# OPT2: choose one
OPT2="--max_step 1000 --plot_int 10 --print_int 10 --hydro_grid_int 0 --n_steps_skip 0 --fixed_dt 1.e-9"
#OPT2="--max_step 10000 --plot_int 0 --print_int 100 --hydro_grid_int 1 --n_steps_skip 1000 --fixed_dt 1.e-9"
#OPT2="--max_step 100000 --plot_int 0 --print_int 1000--hydro_grid_int 1 --n_steps_skip 10000 --fixed_dt 1.e-9"

# OPT3: specify dz
DZ=1.e-5
VC=`python -c "print 1./$DZ"`
OPT3="--variance_coef_mom $VC --variance_coef_mass $VC --cross_section $DZ"

# OPT4: specify rate_multiplier
RM=1.e5
OPT4="--rate_multiplier $RM"

#####

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

OPTS="$OPT1 $OPT2 $OPT3 $OPT4"

echo "mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out"
mpiexec -n 4 ../$EXEC ../$INPUTS $OPTS | tee scr_out

#####

grep sum scr_out | grep 'i=           1' | awk '{print $8}' > rho1
grep sum scr_out | grep 'i=           2' | awk '{print $8}' > rho2
grep sum scr_out | grep 'i=           3' | awk '{print $8}' > rho3
grep sum scr_out | grep 'i=           4' | awk '{print $8}' > rho4
grep sum scr_out | grep 'rho=' | awk '{print $5}' > rhotot
paste rho1 rho2 rho3 rho4 rhotot > res.rho
rm rho1 rho2 rho3 rho4 rhotot

sed -i 's/0.00000000/#0.00000000/g' RGBW.S_k.pair=*.Re.dat
