FC=gfortran
FFLAGS= -fdefault-real-8 -fdefault-double-8 -Wall -funroll-loops -fcheck=bounds -Wno-tabs -g -O3  
POTEN = potendirect.f90
#POTEN = forcebinary.f90
STEP = step_leapfrog.f90
#SETUP = setup_binary.f90
SETUP = plummer_dist.f90
#SETUP = cold_collapse.f90
TREE = contrivedtree.f90
SRC=  read_dump.f90 errors.f90 deltat.f90 octreenode.f90 openingcriterion.f90 octree.f90 contrivedtree.f90 evaluate.f90 taylor.f90 potendirect.f90 interaction.f90   computemass.f90  ${STEP} output.f90 momentum.f90 sphere_dist.f90 ${SETUP} test_direct.f90 test_gravity.f90 main.f90  
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

nbody: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
clean:
	rm *.o *.mod
	rm nbody
cleanruns:
	rm -f snap_*
	rm -f Momentum
	rm -f *.txt


