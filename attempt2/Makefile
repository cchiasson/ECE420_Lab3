
CC = gcc

all: main.c serialtester.c datagen.c
	$(CC) serialtester.c Lab3IO.c -o serialtester -lm
	$(CC) datagen.c Lab3IO.c -o datagen
	$(CC) -g -Wall main.c -o main -lpthread -lm -fopenmp

clean:
	-rm main
	-rm serialtester
	-rm datagen
