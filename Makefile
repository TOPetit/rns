# Makefile

all: mainv timing testing fulltest

fulltest: fulltest.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o fulltest fulltest.c rns.o rnsv.o -lgmp

testing: testing.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o testing testing.c rns.o rnsv.o -lgmp

timing: timing.c rns.o rnsv.o
	gcc -g -march=native -Wno-overflow -o timing timing.c rns.o rnsv.o -lgmp

mainv: mainv.c rns.o rnsv.o
<<<<<<< HEAD
	gcc -g -march=native -Wno-overflow -o mainv mainv.c rns.o rnsv.o -lgmp
=======
	gcc -march=native -Wno-overflow -o mainv mainv.c rns.o rnsv.o -lgmp
>>>>>>> fix_bad_timings

test: test_rns.c rns.o tests.c
	gcc -o rns test_rns.c rns.o -lgmp

rns.o: rns.c rns.h structs_data.h 
	gcc -O3 -Wno-overflow -c rns.c -lgmp 
	
rnsv.o: rnsv.c rns.h structs_data.h 
<<<<<<< HEAD
	gcc -g -march=native -Wno-overflow -c rnsv.c -lgmp 
=======
	gcc -g -O3 -march=native -Wno-overflow -c rnsv.c -lgmp 
>>>>>>> fix_bad_timings

clean:
	rm *.o mainv timing testing fulltest

doc :
	pandoc theory.md -V fontsize=12pt -o theory.pdf

cleandoc:
	rm theory.pdf

