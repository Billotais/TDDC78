#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include <mpi.h>

int main (int argc, char ** argv) {
   int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double stime, etime;
#define MAX_RAD 1000

    double w[MAX_RAD];

    /* Take care of the arguments */

    if (argc != 4) {
	fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
	exit(1);
    }
    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) {
	fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
	exit(1);
    }

	MPI_Init(NULL, NULL);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    /* read file */
    if (myid == 0)
    {
		if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
			exit(1);
        if (colmax > 255) {
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}
	}
   
	MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&colmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    

    printf("Has read the image, generating coefficients\n");

    /* filter */
    get_gauss_weights(radius, w);
    
    printf("Calling filter\n");
 
    //clock_gettime(CLOCK_REALTIME, &stime);
    stime = MPI_Wtime();

    

    if (myid==0) {
      blurfilter(xsize, ysize, src, radius, w);
    } else {
      blurfilter(xsize, ysize, src, radius, w);
    }
	etime = MPI_Wtime();
    printf("Filtering took: %g secs\n", (etime  - stime)) ;

    /* write result */
    

    if (myid == 0) {
        printf("Writing output file\n");
        if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
          exit(1);
    }

    
    MPI_Finalize();

    return(0);
}
