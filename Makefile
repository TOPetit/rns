# Makefile

all: timing testing fulltest main

main: main.c rns.o rns.o
	gcc -g -march=native -Wno-overflow -o main main.c rns.o rnsv.o -lgmp

fulltest: fulltest.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o fulltest fulltest.c rns.o rnsv.o -lgmp

testing: testing.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o testing testing.c rns.o rnsv.o -lgmp

timing: timing.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o timing timing.c rns.o rnsv.o -lgmp

test: test_rns.c rns.o tests.c
	gcc -o rns test_rns.c rns.o -lgmp

rns.o: rns.c rns.h structs_data.h 
	gcc -O3 -Wno-overflow -c rns.c -lgmp 
	
rnsv.o: rnsv.c rns.h structs_data.h 
	gcc -g -O3 -march=native -Wno-overflow -c rnsv.c -lgmp 

clean:
	rm *.o timing testing fulltest main

