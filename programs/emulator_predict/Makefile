# should contain paths to NLOPT and GPF
MOD = ../../../obj
LIB = ../../../lib

FCO = gfortran-7
F90FLAGS = -std=f2008 -Wall -Wextra -O3 $(addprefix -I,$(MOD))
F90LINKFLAGS =  $(addprefix -L,$(LIB)) -lgpf -lnlopt -llapack -lblas

gp_predict: gp_predict.f90
	$(FCO) $(F90FLAGS) -o gp_predict gp_predict.f90 $(F90LINKFLAGS) -fbounds-check
clean:
	rm gp_predict