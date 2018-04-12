
CFLAGS =

LFLAGS= -lpthread -lrt -lm

all: blurc thresc
blur:
	mpicc -Nmpi -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)	
thres:
	mpicc -Nmpi -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)
blur_home:
	mpicc -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)
thres_home:
	mpicc -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)

