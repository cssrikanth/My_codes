test: main.o  glob.o 
	f95 -o test main.o  glob.o 
	



glob.mod:glob.o glob.f95
	f95 -c glob.f95

glob.o:glob.f95
	f95 -c glob.f95

main.o: glob.mod main.f95
	f95 -c main.f95





	    
	
		

clean:

	rm -rf *.o test  *.mod
	


	
