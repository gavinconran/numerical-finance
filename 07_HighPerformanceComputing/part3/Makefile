
OBJECTS = random_util.o mc_walk.o laplace_mc.o
FFLAGS = 

.PHONY: test clean

test: test.exe
	./test.exe

test.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

mc_laplace_error.txt: test.exe
	./test.exe

plot: mc_laplace_error.txt
	python plot_mc_laplace_error.py

clean:
	rm -f *.o *.exe *.mod

clobber: clean
	rm -f mc_laplace_error.txt mc_laplace_error.png

