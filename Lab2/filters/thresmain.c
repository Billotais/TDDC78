#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include <pthread.h>
#include <math.h>

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel *src = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);
    struct timespec stime, etime;

    /* Take care of the arguments */

    if (argc != 4) {
        fprintf(stderr, "Usage: %s infile outfile nthreads\n", argv[0]);
        exit(1);
    }
    int nthreads = atoi(argv[3]);

    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }

    printf("Has read the image, calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    


    pthread_t* threads = calloc(nthreads, sizeof(pthread_t)); // All threads
    thread_data* thread_datas = calloc(nthreads, sizeof(thread_data)); // All data given to threads

    // Init barrier
    pthread_barrier_t barr;
    pthread_barrier_init(&barr, NULL, nthreads);
    // Init lock
    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);

    // Gobal value use to compute the average value of pixels
    int sum = 0;

    // Create all threads
    int i;
    for (i = 0; i < nthreads; ++i)
    {
        // Decide what each thread will work on
		int num_rows_per_job = ceil((float)ysize/(float)nthreads);
		int from_row = (i*num_rows_per_job);
		int to_row = from_row + num_rows_per_job < ysize ? from_row + num_rows_per_job : ysize;
		
        // Give specific data to each thread
        thread_datas[i] = (thread_data){xsize, ysize, from_row, to_row, &sum, src, &barr, &lock};
		
         // Create the thread itself
		pthread_create(threads + i, NULL, thresfilter, thread_datas + i);
		
	}
    // We now wait until all threads are done executing, using join
	for (i = 0; i < nthreads; ++i)
	{
		pthread_join(threads[i], NULL);
	}
    
    // Don't need the barrier anymore
	pthread_barrier_destroy(&barr);

    clock_gettime(CLOCK_REALTIME, &etime);


    printf("Average filtering time: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");
    
    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      exit(1);

    return(0);
}
