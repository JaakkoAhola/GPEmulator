# should contain paths to NLOPT and GPF
MOD = ../../obj
LIB = ../../lib
LIBUSER = ../lib_user


ifndef $(RUNTYPE)
	RUNTYPE=normal
endif

$(info $$RUNTYPE is [${RUNTYPE}])

ifeq ($(RUNTYPE),normal)
	RUNTYPEFLAGS = -O3
else
	# debug
	RUNTYPEFLAGS = -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbounds-check -g
endif

FCFLAGS = -std=f2008 -Wall -Wextra $(RUNTYPEFLAGS) -I$(MOD) -I../../include -I$(LIBUSER) -I/usr/local/include

FCLINKFLAGS = $(LIBUSER)/readio.o -L$(LIB) -L/usr/local/lib -lgpf -lnlopt -llapack -lblas

gp_predict: gp_predict.f90
	$(FC) -o $@ $^ $(FCFLAGS) $(FCLINKFLAGS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

.PHONY: clean

clean:
	rm -f gp_predict *.o
