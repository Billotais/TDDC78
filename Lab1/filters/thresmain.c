#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include <mpi.h>

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel* src;
    double stime, etime;

    /* Take care of the arguments */

    if (argc != 3) {
	fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
	exit(1);
    }
    
	MPI_Init(NULL, NULL);
	
    int myid, size_mpi;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);
    
    /* read file */
    if (myid == 0)
    {
        // Allocate memory only on worker 0 
        src = (pixel*)calloc(MAX_PIXELS, sizeof(pixel));
		if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
			exit(1);
		if (colmax > 255) {
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}
	} else {
        // dummy source
        src = (pixel*)calloc(1, sizeof(pixel));
    } 

	MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&colmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("Has read the image, calling filter\n");

    stime = MPI_Wtime();

    thresfilter(xsize, ysize, src);

	etime = MPI_Wtime();
    // printf("Filtering took: %g secs\n", (etime - stime)) ;

    double total_time = etime - stime;
    double total_time_all;

    // Get time from all workers
    MPI_Reduce(&total_time, &total_time_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    total_time_all = total_time_all / (float)size_mpi;

    /* write result */
    if (myid == 0) {
        printf("Writing output file\n");
        if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
          exit(1);

        printf("Average filtering time: %g secs\n", total_time_all) ;
    }

    free(src);
    MPI_Finalize();

    return(0);
}
