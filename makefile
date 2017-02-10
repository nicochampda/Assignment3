CFLAGS=-Wall -O3
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: galsim.o file_operations.o graphics.o
	gcc $(CFLAGS) $(INCLUDES) -o galsim file_operations.o graphics.o galsim.o $(LDFLAGS)

galsim.o: galsim.c
	gcc $(CFLAGS) $(INCLUDES) -o galsim.o -c galsim.c

file_operations.o: file_operations/file_operations.c
	gcc $(CFLAGS) -o file_operations.o -c file_operations/file_operations.c

graphics.o: graphics/graphics.c graphics/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics/graphics.c

clean:
	rm -f *.o
	rm -f galsim
	rm -f result.gal
