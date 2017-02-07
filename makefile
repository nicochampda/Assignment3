galsim: galsim.o file_operations.o
	gcc -Wall -o galsim file_operations.o galsim.o

galsim.o: galsim.c
	gcc -Wall -o galsim.o -c galsim.c -lm

file_operations.o: file_operations/file_operations.c
	gcc -Wall -o file_operations.o -c file_operations/file_operations.c

clean:
	rm -f *.o
	rm -f galsim
