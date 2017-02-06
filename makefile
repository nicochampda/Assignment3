galsim.o: galsim.c
	gcc -Wall -o galsim.o galsim.c

clean:
	rm -f galsim.o
