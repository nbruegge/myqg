
F90 = gfortran
F90FLAGS = 

#$(F90) $(F90FLAGS) myqg_module.o solve_poisson_cg.o myqg.o
all: clean \
		 myqg_module.o \
     solve_poisson_cg.o \
     myqg.o
		 $(F90) $(F90FLAGS) myqg_module.o solve_poisson_cg.o myqg.o

clean: 
	@echo "Cleaning everythin up!"
	rm -f *.o *.mod *.out

myqg_module.mod: myqg_module.o myqg_module.f90
	$(F90) $(F90FLAGS) -c myqg_module.f90

myqg_module.o: myqg_module.f90
	$(F90) $(F90FLAGS) -c myqg_module.f90

myqg.o: myqg_module.mod myqg.f90 
	$(F90) $(F90FLAGS) -c myqg.f90

solve_poisson_cg.o: solve_poisson_cg.f90
	$(F90) $(F90FLAGS) -c solve_poisson_cg.f90


