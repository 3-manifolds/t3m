#CFLAGS = -O3 -funroll-loops -march=i686 -mcpu=i686 -I/usr/local/include/python2.2 -I/usr/include/python2.2
#MMX = mmx.o

CFLAGS = -O3 -funroll-loops -I/usr/local/include/python2.2 -I/usr/include/python2.2 
MMX = no_mmx.o

#CC = gcc-3.2
CC = gcc

solvermodule.so: solver.o $(MMX) vertex.o
	$(CC) -shared -o solvermodule.so solver.o $(MMX) vertex.o

test: test0 test1 test2 test3

test0:  test0.c vertex.o $(MMX) vertex.h
	$(CC) -g $(CFLAGS) -o test0 test0.c vertex.o $(MMX) -static -lc

test1:  test1.c vertex.o $(MMX) vertex.h
	$(CC) $(CFLAGS) -o test1 test1.c vertex.o $(MMX) -static -lc

test2:  test2.c vertex.o $(MMX) vertex.h
	$(CC) $(CFLAGS) -o test2 test2.c vertex.o $(MMX) -static -lc

magic5:  magic5.c vertex.o $(MMX) vertex.h
	$(CC) $(CFLAGS) -o magic5 magic5.c vertex.o $(MMX) -static -lc

solver.o: solver.c 

vertex.o: vertex.c vertex.h

mmx.o: mmx.c vertex.h

no_mmx: no_mmx.c vertex.h

clean:
	-rm *.o *.so test0 test1 test2 magic5
