# should contain paths to NLOPT and GPF
MOD = ../../obj
LIB = ../../lib
LIBUSER = ../lib_user

FC = gfortran-7
FCFLAGS = -std=f2008 -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbounds-check -g -I$(MOD) -I../../include -I$(LIBUSER)
# for debug replace -O3 with -fbounds-check -g
FCLINKFLAGS = $(LIBUSER)/readio.o -L$(LIB) -lgpf -lnlopt -llapack -lblas

gp_train: gp_train.o
	$(FC) -o $@ $^ $(FCFLAGS) $(FCLINKFLAGS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

.PHONY: clean

clean:
	rm -f gp_train *.o
