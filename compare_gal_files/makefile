compare_gal_files: compare_gal_files.o file_operations.o
	gcc -Wall -o compare_gal_files compare_gal_files.o file_operations.o -lm

compare_gal_files.o: compare_gal_files.c
	gcc -Wall -o compare_gal_files.o -c compare_gal_files.c -lm

file_operations.o: ../file_operations/file_operations.c ../file_operations/file_operations.h
	gcc -Wall -o file_operations.o -c ../file_operations/file_operations.c 

clean: 
	rm -f *.o
	rm -f compare_gal_files
