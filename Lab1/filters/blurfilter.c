/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdio.h>
#include "blurfilter.h"
#include "ppmio.h"
#include <mpi.h>



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

void blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w){
  int x,y,x2,y2, wi;
  double r,g,b,n, wc;
  pixel dst[MAX_PIXELS];
  
  int myid;
  int size_mpi;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);
  
  
  
  
  
  pixel item;
  MPI_Datatype pixel_mpi;
  int block_length[] = {1, 1, 1};
  MPI_Datatype block_types [] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
  MPI_Aint start, displ[3];
  
  MPI_Address(&item, &start);
  MPI_Address(&item.r, &displ[0]);
  MPI_Address(&item.g, &displ[1]);
  MPI_Address(&item.b, &displ[2]);
  
  displ[0] -= start;
  displ[1] -= start;
  displ[2] -= start;

  MPI_Type_struct(3, block_length, displ, block_types, &pixel_mpi);
  MPI_Type_commit(&pixel_mpi);


  // Scatter horizontally
  
  
  
  
  int num_rows_per_job = ceil(ysize/size_mpi);
  int from_row = (myid*num_rows_per_job);
  int to_row = from_row + num_rows_per_job < ysize ? from_row + num_rows_per_job : ysize;
  
  pixel* rcv_buf = pixel*[num_rows_per_job*xsize];
   
  MPI_Scatterv(src, NULL, NULL, pixel_mpi, num_rows_per_job*xsize, pixel_mpi, 0, MPI_COMM_WORLD);
  
  // TODO Correct indexing
  for (y=from_row; y<to_row; y++) { // For each row
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(rcv_buf, x, y, xsize)->r;
      g = w[0] * pix(rcv_buf, x, y, xsize)->g;
      b = w[0] * pix(rcv_buf, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
		wc = w[wi];
		x2 = x - wi;
		if(x2 >= 0) {
		  r += wc * pix(rcv_buf, x2, y, xsize)->r;
		  g += wc * pix(rcv_buf, x2, y, xsize)->g;
		  b += wc * pix(rcv_buf, x2, y, xsize)->b;
		  n += wc;
		}
		x2 = x + wi;
		if(x2 < xsize) {
		  r += wc * pix(rcv_buf, x2, y, xsize)->r;
		  g += wc * pix(rcv_buf, x2, y, xsize)->g;
		  b += wc * pix(rcv_buf, x2, y, xsize)->b;
		  n += wc;
		}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }
  
  // Gather
  MPI_Scatterv(rcv_buf, num_rows_per_job*xsize, pixel_mpi, src, NULL, NULL, pixel_mpi, 0, MPI_COMM_WORLD);
  
  // Scatter verticaly
  
  int num_cols_per_job = ceil(xsize/size_mpi);
  int from_col = (myid*num_cols_per_job);
  int to_col = from_row + num_col_per_job < xsize ? from_col+ num_cols_per_job : xsize;
  for (y=0; y<ysize; y++)  {
    for (x=from_col; x<to_col; x++) {
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
  
  // Gather
  MPI_Finalize();

}



