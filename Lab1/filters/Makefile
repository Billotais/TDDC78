
CFLAGS = 

LFLAGS= -lpthread -lrt -lm 

all: blurc thresc

clean:
	-$(RM) *.o blur thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	icc -Nmpi -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	cc -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cf - *.c *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
	
blur:
	mpicc -Nmpi -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)
	
thres: 
	mpicc -Nmpi -o thres thresmain.c thresfilter.c ppmio.c
	
