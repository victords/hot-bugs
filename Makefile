ep: main.c
	gcc -Wall -O3 -o ep main.c -pthread -lm

debug:
	gcc -Wall -g -o ep main.c -pthread -lm
