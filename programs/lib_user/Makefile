FCO = gfortran-7
F90FLAGS = -std=f2008 -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbounds-check -g

readio: readio.f90
	$(FCO) -c $(F90FLAGS) readio.f90
clean:
	rm readio.mod readio.o
