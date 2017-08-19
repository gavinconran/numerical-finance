
OBJECTS = random_util.o mc_walk.o laplace_mc.o
FFLAGS = 
NUM_PROCS ?= 4   # default if not specified on command line or env variable

.PHONY: plot test clean clobber

test: test.exe
	mpiexec -n $(NUM_PROCS) test.exe

test.exe: $(OBJECTS)
	mpif90 $(FFLAGS) $(OBJECTS) -o test.exe

%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 


mc_laplace_error.txt: test
#	./test.exe

plot: mc_laplace_error.txt
	python plot_mc_laplace_error.py

clean:
	rm -f *.o *.exe *.mod

clobber: clean
	rm -f mc_laplace_error.txt mc_laplace_error.png

