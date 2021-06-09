# Makefile

all: mainv timing

timing: timing.c rns.o rnsv.o
	gcc -march=native -Wno-overflow -o timing timing.c rns.o rnsv.o -lgmp

mainv: mainv.c rns.o rnsv.o
	gcc -O3 -march=native -Wno-overflow -o mainv mainv.c rns.o rnsv.o -lgmp

test: test_rns.c rns.o tests.c
	gcc -o rns test_rns.c rns.o -lgmp

rns.o: rns.c rns.h structs_data.h 
	gcc -O3 -Wno-overflow -c rns.c -lgmp 
	
rnsv.o: rnsv.c rns.h structs_data.h 
	gcc -O3 -march=native -Wno-overflow -c rnsv.c -lgmp 

clean:
	rm *.o mainv timing
