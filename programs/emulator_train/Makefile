# should contain paths to NLOPT and GPF
MOD = ../../obj
LIB = ../../lib
LIBUSER = ../lib_user

FCFLAGS = -std=f2008 -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbounds-check -g -I$(MOD) -I../../include -I$(LIBUSER) -I/usr/local/include
# for debug replace -O3 with -fbounds-check -g -ffpe-trap=invalid,zero,overflow,underflow,denormal
FCLINKFLAGS = $(LIBUSER)/readio.o -L$(LIB) -L/usr/local/lib -lgpf -lnlopt -llapack -lblas

gp_train: gp_train.o
	$(FC) -o $@ $^ $(FCFLAGS) $(FCLINKFLAGS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

.PHONY: clean

clean:
	rm -f gp_train *.o
