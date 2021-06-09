# Makefile
main: main.c rns.o
	gcc -O3 -o main main.c rns.o -lgmp

test: test_rns.c rns.o tests.c
	gcc -O3 -o rns test_rns.c rns.o -lgmp

rns.o: rns.c rns.h structs_data.h 
	gcc -O3 -c rns.c -lgmp 

clean:
	rm *.o main

