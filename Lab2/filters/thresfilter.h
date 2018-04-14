/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

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
    int* sum;
    pixel* src;
    pthread_barrier_t* barr;
    pthread_mutex_t* lock;
} thread_data;

void* thresfilter(void * args);


#endif
