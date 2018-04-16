/*
  File: blurfilter.c

  Implementation of blurfilter function.    
 */
#include <stdio.h>
#include "blurfilter.h"




pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

void* blurfilter(void * args){
  int x,y,x2,y2, wi;
  double r,g,b,n, wc;
  
	// "Rename" all the data given to the thread for easier use
	thread_data* data = (thread_data*) args;
	int xsize = data->xsize;
	int ysize = data->ysize;
	int start = data->start;
	int end = data->end;
	pixel* src = data->src;
	pixel* dst = data->dst;
	int radius = data->radius;
	const double *w = data->w;

	// Compute horizontal blur
  for (y=start; y<end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(src, x, y, xsize)->r;
      g = w[0] * pix(src, x, y, xsize)->g;
      b = w[0] * pix(src, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
				wc = w[wi];
				x2 = x - wi;
				if(x2 >= 0) {
					r += wc * pix(src, x2, y, xsize)->r;
					g += wc * pix(src, x2, y, xsize)->g;
					b += wc * pix(src, x2, y, xsize)->b;
					n += wc;
				}
				x2 = x + wi;
				if(x2 < xsize) {
					r += wc * pix(src, x2, y, xsize)->r;
					g += wc * pix(src, x2, y, xsize)->g;
					b += wc * pix(src, x2, y, xsize)->b;
					n += wc;
				}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }

	// Wait until all threads are at this point
	pthread_barrier_wait(data->barr);

	// Do vertical blur
  for (y=start; y<end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(dst, x, y, xsize)->r;
      g = w[0] * pix(dst, x, y, xsize)->g;
      b = w[0] * pix(dst, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
				wc = w[wi];
				y2 = y - wi;
				if(y2 >= 0) {
					r += wc * pix(dst, x, y2, xsize)->r;
					g += wc * pix(dst, x, y2, xsize)->g;
					b += wc * pix(dst, x, y2, xsize)->b;
					n += wc;
				}
				y2 = y + wi;
				if(y2 < ysize) {
					r += wc * pix(dst, x, y2, xsize)->r;
					g += wc * pix(dst, x, y2, xsize)->g;
					b += wc * pix(dst, x, y2, xsize)->b;
					n += wc;
				}
      }
      pix(src,x,y, xsize)->r = r/n;
      pix(src,x,y, xsize)->g = g/n;
      pix(src,x,y, xsize)->b = b/n;
    }
  }
}



