FC=gfortran
FFLAGS= -fdefault-real-8 -fdefault-double-8 -fbounds-check -Wall -O3 
#POTEN = potendirect.f90
POTEN = forcebinary.f90
STEP = step_leapfrog.f90
SETUP = setup_binary.f90
SRC= octreenode.f90 octree.f90 computemass.f90  ${POTEN} ${STEP} output.f90 momentum.f90 sphere_dist.f90 ${SETUP}  main.f90 
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

nbody: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	rm *.o *.mod
	rm nbody
cleanrun:
	rm -f snap_*
	rm -f Momentum


