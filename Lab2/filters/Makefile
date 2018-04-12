
CFLAGS =

LFLAGS= -lpthread -lrt -lm

all: blurc thresc
blur:
	rm -f blur
	gcc -o blur blurfilter.c blurmain.c gaussw.c ppmio.c $(LFLAGS)	
thres:
	gcc -o thres thresmain.c thresfilter.c ppmio.c $(LFLAGS)

