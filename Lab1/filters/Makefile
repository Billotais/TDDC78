
CFLAGS =

LFLAGS= -lpthread -lrt -lm

all: blurc thresc
blur: blurfilter.c blurmain.c gaussw.c ppmio.c
	mpicc -Nmpi -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)	
thres: thresmain.c thresfilter.c ppmio.c
	mpicc -Nmpi -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)
blur_home: blurfilter.c blurmain.c gaussw.c ppmio.c
	mpicc -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)
thres_home: thresmain.c thresfilter.c ppmio.c
	mpicc -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)

