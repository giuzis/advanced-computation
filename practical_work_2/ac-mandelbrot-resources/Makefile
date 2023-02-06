CC=gcc
CC2=mpicc
MFLAGS=-lm
GFLAGS=-lGL -lGLU  -lglut

#use the 1st line in alternative to the 2nd (1st is for valgrind, 2nd is for benchmarking)
CFLAGS=-Wall -O0 -g
# CFLAGS=-Wall -O0
# CFLAGS=-Wall -O2
MPIFLAGS=-I/share/apps/openmpi-4.1.4/include

all: mandelbrot-gui-serial.exe mandelbrot-gui-mpi.exe mandelbrot-gui-test.exe

mandelbrot-gui-serial.exe: mandelbrot-gui-serial.c
	$(CC) $(CFLAGS) mandelbrot-gui-serial.c $(MFLAGS) $(GFLAGS) -o mandelbrot-gui-serial.exe 

mandelbrot-gui-test.exe: mandelbrot-gui-test.c
	$(CC) $(CFLAGS) mandelbrot-gui-test.c $(MFLAGS) $(GFLAGS) -o mandelbrot-gui-test.exe 

mandelbrot-gui-mpi.exe: mandelbrot-gui-mpi.c
	$(CC2) $(CFLAGS) $(MPIFLAGS) mandelbrot-gui-mpi.c $(MFLAGS) $(GFLAGS) -o mandelbrot-gui-mpi.exe 

clean:
	rm -f *.exe *.o a.out

