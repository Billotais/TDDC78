#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include <pthread.h>

#include <math.h>


int main (int argc, char ** argv) {
   int radius;
    int xsize, ysize, colmax;
    pixel *src = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);
    pixel *dst = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);
    struct timespec stime, etime;
	#define MAX_RAD 1000

    double w[MAX_RAD];

    /* Take care of the arguments */

    if (argc != 5) {
        fprintf(stderr, "Usage: %s radius infile outfile nb_threads\n", argv[0]);
        exit(1);
    }
    radius = atoi(argv[1]);
    int nthreads = atoi(argv[4]);
    if((radius > MAX_RAD) || (radius < 1)) {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        exit(1);
    }

    /* read file */
    if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }

    printf("Has read the image, generating coefficients\n");

    /* filter */
    get_gauss_weights(radius, w);

    printf("Calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);
    
    // Splitting => nb threads
    
    
    pthread_t* threads = calloc(nthreads, sizeof(pthread_t)); // All threads
    thread_data* thread_datas = calloc(nthreads, sizeof(thread_data)); // All data given to threads

    // Init barrier
    pthread_barrier_t barr;
    pthread_barrier_init(&barr, NULL, nthreads);

    // Create all threads
    int i;
    for (i = 0; i < nthreads; ++i)
    {
        // Decide what each thread will work on
		int num_rows_per_job = ceil((float)ysize/(float)nthreads);
		int from_row = (i*num_rows_per_job);
		int to_row = from_row + num_rows_per_job < ysize ? from_row + num_rows_per_job : ysize;
		
        // Give specific data to each thread
        thread_datas[i] = (thread_data){xsize, ysize, from_row, to_row, src, dst, radius, w, &barr};
		
        // Create the thread itself
		pthread_create(threads + i, NULL, blurfilter, thread_datas + i);
		
	}
    // We now wait until all threads are done executing, using join
	for (i = 0; i < nthreads; ++i)
	{
		pthread_join(threads[i], NULL);
	}
    
    // Don't need the barrier anymore
	pthread_barrier_destroy(&barr);


    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    
    free(threads);
    free(thread_datas);

    printf("Writing output file\n");
    if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
      exit(1);


    return(0);
}
