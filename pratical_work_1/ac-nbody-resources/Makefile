CC=gcc
MFLAGS=-lm
TFLAGS=-pthread

#uncomment only if SDL bgi is installed
GFLAGS=-lSDL_bgi -lSDL2

#use the 1st line in alternative to the 2nd (1st is for valgrind, 2nd is for benchmarking)
CFLAGS=-Wall -O0
#CFLAGS=-Wall -O2

#use the 1st line in alternative to the 2nd (use 2nd only if SDL bgi installed)
#all: nbody-serial.exe nbody-threads.exe 
all: nbody-serial.exe nbody-serial-opt.exe nbody-threads.exe  nbody-serial-gui.exe calc_chunk.exe

calc_chunk.exe: calc_chunk.c
	$(CC) $(CFLAGS) calc_chunk.c $(MFLAGS) -o calc_chunk.exe 

nbody-serial.exe: nbody-serial.c
	$(CC) $(CFLAGS) nbody-serial.c $(MFLAGS) -o nbody-serial.exe 

nbody-serial-opt.exe: nbody-serial-opt.c
	$(CC) $(CFLAGS) nbody-serial-opt.c $(MFLAGS) -o nbody-serial-opt.exe 

#nbody-threads.exe: nbody-threads.c
#	$(CC) $(CFLAGS) nbody-threads.c $(MFLAGS) $(TFLAGS) -o nbody-threads.exe

#uncomment the following two lines only if SDL bgi is installed
#nbody-serial-gui.exe: nbody-serial-gui.c
#	$(CC) $(CFLAGS) nbody-serial-gui.c $(MFLAGS) $(GFLAGS) -o nbody-serial-gui.exe

clean:
	rm -f *.exe *.o a.out
