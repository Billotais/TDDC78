#include "thresfilter.h"
#include <stdio.h>

void* thresfilter(void * args){
#define uint unsigned int 

  uint sum_local, i, psum;

  // Get data 
	thread_data* data = (thread_data*) args;
	int xsize = data->xsize;
	int ysize = data->ysize;
	int start = data->start;
	int end = data->end;
  int* sum = data->sum;
	pixel* src = data->src;

  // Compute sum for region
  for(i = start*xsize, sum_local = 0; i < end*xsize; i++) {
    sum_local += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  // Add it to common sum
  pthread_mutex_lock(data->lock);
  *sum += sum_local;
  pthread_mutex_unlock(data->lock);

  // wait until everyone has added its sum 
  pthread_barrier_wait(data->barr);

  // Get the average value
  int mean = (*data->sum)/(xsize*ysize);

  // Apply the threshold
  for(i = start*xsize; i < end*xsize; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(mean > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}
