all : modular timing

modular :
	gcc -g -march=native -Wno-overflow -o modular modular_mult.c -lgmp

timing: timing.c rns.o
	gcc -g -march=native -Wno-overflow -o timing timing.c rns.o -lgmp

rns.o: rns.c rns.h structs_data.h 
	gcc -Wno-overflow -c rns.c

clean :
	rm modular timing *.o

doc :
	pandoc theory.md -o theory.pdf

cleandoc:
	rm theory.pdf