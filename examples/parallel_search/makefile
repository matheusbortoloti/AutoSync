.PHONY: all build test clean


all: build test

build: 
	gcc main_posix.c -o MainPosix.o -lpthread

test: build
	./MainPosix.o	

helgrind:
	valgrind --tool=helgrind ./MainPosix.o

clean:
	rm -rf *o *out
