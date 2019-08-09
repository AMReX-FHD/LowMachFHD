# HydroGrid library for computing static and dynamic structure factors (spectra)


## Usage

First try to run the python example of diffusion in a binary mixture of uncorrelated Brownian walkers, which is in the folder/directory example.

* cd example
* cp ../src/MakefileHeader .
* Edit MakefileHeader to edit compilers/paths/etc
* make
* Edit hydroGridOptions.nml where all options for HydroGrid are set, and then run
* python brownian_walkers.py

In this example we have 

num_particles = 1024*16

points diffusing as **independent** Brownian walkers in the plane, with diffusion coefficient 1 in a domain of size 256^2 covered by 128^2 grid cells. The first half the particles is "green" (species 1) and the rest "red" (species 2). In the initial configuration all green particles are in the top half, y>Ly/2, and the red in the bottom half, y<=Ly/2, and they mix diffusively.

Every time step the example pyhon script calls the HydroGrid library to compute structure factors (spectra). Rigth now the code outputs all spectra at the end of the run, but if you change "sample" in brownian_walkers.py to a positive value, say 1000, it will output results every 1000 steps. **Each time HydroGrid outputs statistics it resets all internal variables and starts anew.** Since our configurations evolves dynamically, this is useful to output average statistics over a certain period of time where the statistics don't change much.

## libHydroGrid

The main input file that controls libHydroGrid is the Fortran namelist file:
**hydroGridOptions.nml**
This is basically a text input file with ! for comments.
It contains two namelists, contained between the delimiters:

&hydroAnalysisOptions

   ... ! All options (if excluded they take default values)

/

The first namelist corresponds to doing analysis on the full 2D grid. When initialized it prints:

 Initializing HydroAnalysis for grid of size          128         128           1  nSpecies/isSingleFluid=           3 T

This code will also, however, project the grid along the y axes to get averages that are now just a function of x. When initialized this prints:

 Initializing HydroAnalysis for grid of size          128           1           1  nSpecies/isSingleFluid=           3 T
 
Spectra are then computed also for this 1D grid. Once can separately control the 2D and 1D analysis via the two separate namelists. I have selected the relevant options for you. The output files for the 2D grid have a prefix "run" set by:

   filePrefix = "run"

and the output files for the projected 1D grid have file names with "run_proj" since the second namelists sets:

   filePrefix = "run_proj"

## Variable numbering

HydroGrid numbers the different variables and can compute spectra of one variable or cross-correlations between variables. The numbering for the present 2D setup is:

1=rho (total density), 2-4 (unused = vx, vy, temperature), 5=rho_green, 6=rho_red

where rho here means number density of particles. So for us relevant variables are:

**variable: 1=rho, 5=rho_green and 6=rho_red**

Every time statistics are written out, the code also outputs mean values averaged along the x axes, i.e., as a function of y. This is useful for steady-state runs to very accurately compute the average profiles; for dynamic runs it can be a useful diagnostic but the result is not very useful quantitatively because it is averaged over a number of time samples. For example, to see the average profile of green and red particles at the end of the run, rho_red(y) and rho_green(y), plot columns 6 and 7 of the .means.dat file (first column is y, so variable=5 is column 6 in the file):

xmgrace -block run.means.inst.dat -bxy 1:6 -bxy 1:7

# Structure Factors

The main purpose of HydroGrid is to compute static and dynamic correlations in Fourier space, i.e., static and dynamic structure factors. It does this by computing the FFT of variables and then computing cross-correlations between different variables, either static ones S(k) or time-correlation functions S(k,t).

When computing structure factors, one can select as many **pairs** of variables as desired. For example, for the 1D grid we set:

   nStructureFactors = 4
  
   structureFactorPairs = "0 1 0 1, 0 5 0 5, 0 6 0 6, 0 5 0 6"

Note: The zeros have to be here but mean nothing to us, so (0,5) means variable 5.

This tells HydroGrid to compute S(k) for variable pair (1,1), i.e., to compute rho-rho static correlation in Fourier space, as well as variable pairs (5,5), (6,6) and (5,6). HydroGrid prints:

 nVariablesToFFT=           3  variablesToFFT:            0           1           0           5           0           6

 grid%structureFactorPairs=           1           1           2           2           3           3           2           3

 Min k =          -63           0           0

 Max k =           64           0           0

In HydroGrid wavenumbers here are described by integer indices kappa, so

wavenumber=(kappa_x, kappa_y)

means

k=(2*pi/Lx*kappa_x, 2*pi/Lx*kappa_x)

## Static structure factors

The output file

run.S_k.pair=1.Re.dat

contains the static structure factor S(k) for the first pair of variables, i.e., here for rho-rho, real part. You can just use xmgrace to plot it. Everything is normalized including x and y axes. This spectrum should be flat even for true-2D hydro if the total density of particles is uniform. If the option 

writeSpectrumVTK = T

is set then HydroGrid will also write S(kx,ky) as a vtk file that can be looked at it in visit.

For the projected 1D grid, we also compute separately spectra for green, red, and also correlations red-green:

   nStructureFactors = 4

   structureFactorPairs = "0 1 0 1, 0 5 0 5, 0 6 0 6, 0 5 0 6"

so the file 

run_proj.S_k.pair=2.Re.dat

contains the spectrum for the green particles (variable pair is (5,5)). This is the one that shows giant fluctuations for hydrodynamically-correlated walkers, but is flat here since the walkers are uncorrelated.

## Dynamic Structure Factors

HydroGrid can also compute dynamic structure factors, i.e., S(k,t) or S(k,w). This is done for all pairs of particles for which static factors S(k) are computed, but only for a group of **selected wavenumbers**. This is because HydroGrid keeps a history of S(k) for those selected k's. The length of the history (in terms of number of calls to HydroGrid) is controlled by nSavedSnapshots in the namelist. So if you set nSavedSnapshots=100 this means that every 100 calls to updateHydroGrid the library will compute time correlation functions. These are then averaged over **blocks of nSavedSnapshots spectra** until you call saveHydroGrid to write the average S(k,w/t) to a file.

The output is controlled by the option:

writeDynamicFiles=1 to write S(k,w), 0 to write individual k's, or -1 for one big file

Only if writeDynamicFiles=1 is S(k,w) written to a file (separately for the real and imaginary parts). Otherwise, only time correlation functions S(k,t) (note that these are always real-valued) are written since S(k,t) is computed more accurately than S(k,w) as described below. One can choose whether to write one big file with all k's as columns for each structure factor pair (writeDynamicFiles=-1), or one file per wavenumber k but with all structure factor pairs written as columns (writeDynamicFiles=0).

When computing the time correlation functions, HydroGrid uses the FFT and **doubles the length of the time history** (i.e., of the block) and then padds with zeros to avoid introducing periodic artifacts due to the FFT (this corresponds to the default value of the option paddDynamicFFT=0). This means that the smallest frequency w for which S(k,w) is written (if writeDynamicFiles=1) is:

w_min = 2*Pi/(2*nSavedSnapshots*dtime)

where dtime is the time interval between successive calls to updateHydroGrid. Similarly, the smallest time for which S(k,t) is written is

t_min = dtime

and the largest is

t_max=nSavedSnapshots*dtime

For time correlation functions S(k,t), HydroGrid implements an unbiased estimator that is equivalent to a standard double loop over all time lags and points but using an FFT to accelerate the computation. However, there appears to be no simple way to obtain an unbiased estimator of S(k,w) without using windowing and other complicated tricks. If one wants S(k,w) and not S(k,t) then one can set paddDynamicFFT=-1 which does not padd the sequence with zeros and just does FFTs assuming the sequence is periodic, which will introduce some bias.

### Wavenumber selection

To select how many wavenumbers to compute dynamic factors for, use the input nWavenumbers in the namelist. If this is set to zero, dynamic factors won't be computed. If set to a positive number, then you need to give an explicit list of the wavenumbers you want to track, for example:

selectedWavenumbers="2 16, 32 32"

will track kappa=(2,16) and kappa=(32,32).

If nWavenumbers is negative, then you specify a block of wavenumbers as a subgrid of the full k-space grid to track. This is useful if you want to track, for example, all small wavenumbers but not track the largest wavenumbers to save memory and time. For example

selectedWavenumbers="-4 4, 0 4"

means that we track wavenumbers with kappa_x in the interval (-4:4), and kappa_y in the interval (0,4). This makes for a total of (2*4+1)x(4+1)=9x5=45 wavenumbers to track. Indeed, HydroGrid prints 

 Keeping track of S(k,omega) and S(k,t) for           45  wavenumbers

and outputs a list of the wavenubers it will track:

 Tracking wavenumber index=           7  k_index=           2           0           0
 
This means that output files with name containing k=7 pertain to wavenumber k=(2*2*Pi/Lx, 0*2*Pi/Ly).

For writeDynamicFiles=1, the actual result for S(k,w) and S(k,t) are written to the files:

run.S_k_t.k=7.dat

as a function of time, and

run.S_k_w.k=7.{Re,Im}.dat

as a function of w. Note that the first line of this file shows the wavenumber k so you can make sure you know what k you are looking at. If writeDynamicFiles=-1, then the code writes one big file run.S_k_t.pair=?.dat for each structure factor pair, with all the k's as columns. The order of the k's in the file is the same as the numbering printed in the beginning. The k=2*Pi*(i/Lx,j/Ly,k/Lz) are ordered on a grid and the loop is k outermost, j inner, and i innermost (Fortran order). There is one catch however -- the library skips zero wavenumber k=(0,0,0) since this is time-independent for conservation laws and usually behaves very differently from other wavenumbers.
