SMEx_for (Sample Moments & Extremes fortran package)

Set of simple codes to compute sample central moments and extremes of grid-point daily temperature data
read from netcdf files:

i. Climax-fnc.90 returns day-of-the-year climate normals referenced to preindustrial period;

ii. CMoments-fnc.f90 returns annual/seasonal sample moments of daily anomalies (mean, standard deviation,
skewness coefficient and excess kurtosis) for every grid point and year in the input period ;
 
iii. EXtremes-fnc.f90 returns annual/seasonal extreme values (minimum and maximum) and probabilities of
exceeding high preindustrial quantiles, for every grid point and year as above.

Programs are adjusted for climate model outputs (CMIP5/6).


REQUIREMENTS

Linux/MacOS operating system;
gfortran compiler;
netcdf library (https://github.com/Unidata/netcdf-fortran);
dfftpack library (https://github.com/fortran-lang/fftpack);
quicksort_real_F77.F (https://github.com/jasonanema/Quicksort_Fortran77) 


COMPILATION

The 'Makefile' template in the source directory /src can be used to compile main programs with dependencies,
upon setting user parameters (compilation options, destination directory and library paths).
Programs have been successfully compiled and tested on Linux (Debian 4.19) using gfortran 8.3.


USAGE

Programs take as input a user supplied list of the nc files to be processed together with the list of begin
and end year of each file. Input files are time ordered and have the same variables' and grid structure.
'Execute-fnc' is a shell script template to pass input to main programs and start execution.
Results are stored as single precision real into ascii files for each grid point. 

This set of codes was used to produce results discussed in Simolo, C., Corti, S. 'Quantifying the role of
variability in future intensification of heat extremes'. Nat Commun 13, 7930 (2022).
https://doi.org/10.1038/s41467-022-35571-0

