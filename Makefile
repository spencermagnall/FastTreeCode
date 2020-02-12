FC=gfortran
FFLags= -O3
#POTEN= poten.f90
#POTEN= poten-isochrone.f90 # choice of potential used in the code
SRC= octreenode.f90 octree.f90 main.f90 
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

nbody: $(OBJ)
	$(FC) $(FFLAGs) -o $@ $(OBJ)
clean:
	rm *.o *.mod


