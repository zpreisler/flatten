CC=gcc
CLIBS=-lm
EXEC=flatten
MAIN=main.c

all: $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) -o $(EXEC) $(CLIBS)
