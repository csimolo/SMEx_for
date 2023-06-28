#
# makefile template
#

# SET your compilation options, destination directory for exec files and library paths
#
fc := gfortran
fflags := 
AR := ar rcs
RM := rm -f
WKD := /..path../wkdir/
#
# netcdf (https://github.com/Unidata/netcdf-fortran)
libnc := -lnetcdf -lnetcdff 
# fftpack (https://github.com/fortran-lang/fftpack) -- used by Climax-fnc.f90
libft := /..path../fftpack/src/libdfftpack.a
# quicksort_real_F77.F (https://github.com/jasonanema/Quicksort_Fortran77) -- used by Calc_quantiles.f90
sort := quicksort_real_F77.F

# list source files
#
src0 := Calc_time.f90 Fix_parameters.f90 Calc_days.f90 Calc_ctmoms.f90 
src1 := compare.f90 $(sort) Calc_quantiles.f90
main1 := Climax-fnc.f90
main2 := CMoments-fnc.f90
main3 := EXtremes-fnc.f90

# define objects and the archive of auxiliary routines
#
obj0 := $(src0:.f90=.o) 
obj1 := compare.o $(sort:.F=.o) Calc_quantiles.o 
libar :=auxlib.a

# define executables
#
Climax := $(WKD)$(main1:.f90=) 
CMoments := $(WKD)$(main2:.f90=) 
EXtremes := $(WKD)$(main3:.f90=) 

.PHONY: all clean 
all: $(Climax) $(CMoments) $(EXtremes) 
#

# Create objects from source and the archive lib
#
%.o: %.f90
	$(fc) -c $< 

%.o: %.F
	$(fc) -c $< 

$(libar): $(obj0) $(obj1)
	$(AR) $@ $^

# Create executables in WKD
#
$(Climax): $(main1) $(libar)
	$(fc) -o $@ $< $(libar) $(libft) $(fflags) $(libnc) 

$(CMoments): $(main2) $(libar)
	$(fc) -o $@ $< $(libar) $(fflags) $(libnc) 

$(EXtremes): $(main3) $(libar)
	$(fc) -o $@ $< $(libar) $(fflags) $(libnc) 

# Cleanup all
#
clean:
	$(RM) *.o *.mod *.a *~
#
