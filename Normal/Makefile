#CFLAGS = -O2 -g
CFLAGS = -O3 -march=i686 -mcpu=i686 -I/usr/local/include/python2.2

_normalmodule.so: _normal.o mmx.o vertex.o
	gcc -shared -o _normalmodule.so _normal.o mmx.o vertex.o

test: test0 test1 test2 test3

test0:  test0.c vertex.o mmx.o vertex.h
	gcc $(CFLAGS) -o test0 test0.c vertex.o mmx.o -static -lc

test1:  test1.c vertex.o mmx.o vertex.h
	gcc $(CFLAGS) -o test1 test1.c vertex.o mmx.o -static -lc

test2:  test2.c vertex.o mmx.o vertex.h
	gcc $(CFLAGS) -o test2 test2.c vertex.o mmx.o -static -lc

test3:  test3.c vertex.o mmx.o vertex.h
	gcc $(CFLAGS) -o test3 test3.c vertex.o mmx.o -static -lc

_normal.o: _normal.c 

vertex.o: vertex.c vertex.h

mmx.o: mmx.c vertex.h

clean:
	-rm *.o *.so test0 test1 test2 test3
