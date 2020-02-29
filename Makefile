FC=gfortran
FFLags= -fbounds-check -Wall -ffpe-trap=list
POTEN = potendirect.f90
STEP = step_leapfrog.f90
SETUP = setup_binary.f90
SRC= octreenode.f90 octree.f90  ${POTEN} ${STEP} output.f90 momentum.f90 sphere_dist.f90 ${SETUP}  main.f90 
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

nbody: $(OBJ)
	$(FC) $(FFLAGs) -o $@ $(OBJ)
clean:
	rm *.o *.mod
	rm nbody


