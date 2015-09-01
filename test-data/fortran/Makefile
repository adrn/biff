# compilers, flags and libraries
FC  = gfortran
FCFLAGS = -ffixed-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -fdefault-real-8

# source files and objects
SRC = bfe.f
PROG = bfe.o
COEFFFILES = $(wildcard *.coeff)

all:
	$(FC) $(FCFLAGS) $(SRC) -o $(PROG)

run:
	./$(PROG) wang-zhao.coeff wang-zhao-accp.dat
	./$(PROG) simple-hernquist.coeff simple-hernquist-accp.dat

clean:
	rm $(PROG) *-accp.dat
