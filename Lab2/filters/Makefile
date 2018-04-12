
CFLAGS =

LFLAGS= -lpthread -lrt -lm

all: blurc thresc
blur: blurfilter.c blurmain.c gaussw.c ppmio.c
	gcc -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)	
thres: thresmain.c thresfilter.c ppmio.c
	gcc -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)

