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
	cp ../biff/data/positions.dat.gz .
	cp ../biff/data/plummer-pos.dat.gz .

	gzip -d positions.dat.gz
	gzip -d plummer-pos.dat.gz

	./$(PROG) wang-zhao.coeff wang-zhao-accp.dat
	gzip wang-zhao-accp.dat

	./$(PROG) simple-hernquist.coeff simple-hernquist-accp.dat
	gzip simple-hernquist-accp.dat

	./$(PROG) random.coeff random-accp.dat
	gzip random-accp.dat

	rm positions.dat
	rm plummer-pos.dat

clean:
	rm $(PROG) *-accp.dat.gz
