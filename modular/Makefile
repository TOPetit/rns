all : modular

modular :
	gcc -g -march=native -Wno-overflow -o modular modular_mult.c -lgmp

clean :
	rm modular

doc :
	pandoc theory.md -o theory.pdf