CFLAGS = -O2 -Wall -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations 
PYTHON_VERSION = 2.2
KERNEL_HEADERS = -I../headers/ -I../unix_kit/ 
HEADERS = -Iheaders -Iunix_kit -I/usr/local/include/python$(PYTHON_VERSION)  -I/usr/include/python$(PYTHON_VERSION) 

all: KernelObjects/BuildDate SnapPeaC.so

install: SnapPeaC.so
	cp SnapPeaC.so ..

KernelObjects/BuildDate:	code/*.c headers/*.h unix_kit/*
	mkdir -p KernelObjects
	cd KernelObjects ;\
	date > BuildDate ;\
	gcc -c $(CFLAGS) $(KERNEL_HEADERS) -DDEBUG_MALLOC=0 \
		../code/*.c				\
		../unix_kit/unix_UI.c			\
		../unix_kit/unix_cusped_census.c	\
		../unix_kit/unix_file_io.c              \
#		../unix_kit/OldFileFormat.hack.c        

SnapPeaC.so:    SnapPeaC.o KernelObjects/BuildDate decode_new_DT.o
	gcc -shared -o SnapPeaC.so $(CFLAGS)    \
	SnapPeaC.o		\
	KernelObjects/*.o	\
	decode_new_DT.o 	\
	-lm

SnapPeaC.o: SnapPeaC.c
	gcc -c $(CFLAGS) $(HEADERS) SnapPeaC.c

decode_new_DT.o: decode_new_DT.c
	gcc -c $(CFLAGS) $(HEADERS) decode_new_DT.c

clean:
	rm -rf KernelObjects
	rm SnapPeaC.o SnapPeaC.so decode_new_DT.o