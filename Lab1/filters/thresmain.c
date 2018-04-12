#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include <mpi.h>

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double stime, etime;

    /* Take care of the arguments */

    if (argc != 3) {
	fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
	exit(1);
    }
    
	MPI_Init(NULL, NULL);
	
	int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    /* read file */
    if (myid == 0)
    {
		if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
			exit(1);
		if (colmax > 255) {
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}
	}
	MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&colmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

  

    printf("Has read the image, calling filter\n");
   
    
    stime = MPI_Wtime();

    
    if (myid==0) {
      thresfilter(xsize, ysize, src);
    } else {
      //pixel dummy_src[1];
      thresfilter(xsize, ysize, src);
    }

 
    //clock_gettime(CLOCK_REALTIME, &etime);
	
	etime = MPI_Wtime();
    printf("Filtering took: %g secs\n", (etime - stime)) ;

    /* write result */
    printf("Writing output file\n");
    
    if (myid == 0) {
        printf("Writing output file\n");
        if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
          exit(1);
    }
    MPI_Finalize();


    return(0);
}
