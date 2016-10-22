
OBJECTS = functions.o quadrature.o test.o
OBJECTS2 = functions.o quadrature2.o test2.o
FFLAGS = -fopenmp

.PHONY: test test2 clean 

test: test.exe
	./test.exe

test2: test2.exe
	./test2.exe

test.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test.exe

test2.exe: $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test2.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

