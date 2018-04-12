/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

#include <pthread.h>

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

typedef struct {
    int xsize;
    int ysize;
    int start;
    int end;
    pixel* src;
    pixel* dst;
    int radius;
    const double *w;
    pthread_barrier_t* barr;
} thread_data;

void* blurfilter(void * args);

#endif
