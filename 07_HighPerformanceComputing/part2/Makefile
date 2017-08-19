
OBJECTS = random_util.o functions.o quadrature_mc.o test_quad_mc.o
FFLAGS = 

.PHONY: test plot clean clobber

test: test.exe
	./test.exe

test.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

mc_quad_error.txt: test.exe
	./test.exe

plot: mc_quad_error.txt
	python plot_mc_quad_error.py

clean:
	rm -f *.o *.exe *.mod

clobber: clean
	rm -f mc_quad_error.txt mc_quad_error.png
