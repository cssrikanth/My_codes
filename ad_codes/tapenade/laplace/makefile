laplace: main.o solver.o solver_d.o
	f95 -o laplace main.o solver.o solver_d.o

	



main.o: main.f95
	f95 -c main.f95

solver.o : solver.f95
	f95 -c solver.f95

solver_d.o: solver_d.f95
	f95 -c solver_d.f95


	    
	
		

clean:

	rm -rf *.o laplace
	


	
