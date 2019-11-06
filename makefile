MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 				 = $(wildcard *.cpp Solid/*.cpp Solid/Fiber/*.cpp)
FCOMPILER        = gfortran -O2
CXXFLAGS        += -w
CURRENT_DIR      = $(shell pwd)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

s:$(GCH_FILES:.h=.gch) $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -lboost_system -std=c++0x

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

clear:
	@$ rm *.o *~ s *.vtu mirror* domain* *.mod *.dat ma26* tensao* esforc* saida omega.txt *.geo *.msh $(CURRENT_DIR)/src/*.gch

run1:
	@$ mpirun -np 1 ./s

run2:
	@$ mpirun -np 2 ./s

run3:
	@$ mpirun -np 3 ./s

run4:
	@$ export OMP_NUM_THREADS=1 && mpirun -n 4 ./s

run5:
	@$ mpirun -np 5 ./s

run6:
	@$ mpirun -np 6 ./s

run7:
	@$ mpirun -np 7 ./s

run8:
	@$ mpirun -np 8 ./s -pc_factor_nonzeros_along_diagonal 1.e-8

run12:
	@$ mpirun -np 12 ./s
#	@$ mpirun -np 16 ./f -pc_type jacobi -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000
