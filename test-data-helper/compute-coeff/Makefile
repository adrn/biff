# compilers, flags and libraries
FC  = gfortran
FCFLAGS = -ffixed-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -fdefault-real-8

# source files and objects
SRC = computecoeff.f
PROG = computecoeff.o

all:
	$(FC) $(FCFLAGS) $(SRC) -o $(PROG)

run:
	cp ../../biff/scf/data/hernquist-samples.dat.gz .

	gzip -d hernquist-samples.dat.gz

	./$(PROG) hernquist-samples.dat computed-hernquist.coeff

	rm hernquist-samples.dat

clean:
	rm $(PROG) *-samples.dat.gz *.dat *.coeff
