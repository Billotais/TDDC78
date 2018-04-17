/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdlib.h>
#include <stdio.h>
#include "blurfilter.h"
#include "ppmio.h"
#include <mpi.h>
#include <math.h>

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

  
  int myid;
  int size_mpi;

  // Get mpi info
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &size_mpi);
  
  // mpi datatype
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
  MPI_Status status;

  // The image is split horizontally into blocks, and each worker works
  // on one block

  int num_rows_per_job = ceil((float)ysize/(float)size_mpi);
  int from_row = (myid*num_rows_per_job);
  int to_row = from_row + num_rows_per_job < ysize ? from_row + num_rows_per_job : ysize;
  
  pixel* rcv_buf = (pixel*)calloc(num_rows_per_job*xsize, sizeof(pixel));
  
  int counts[size_mpi];
  int offsets[size_mpi];
  
  pixel* prov_dst = (pixel*)calloc(xsize*(num_rows_per_job  + 2*radius), sizeof(pixel));

  int i = 0;
  for (i = 0; i < size_mpi-1; ++i)
  {
	  counts[i] = num_rows_per_job*xsize;
  } 
  counts[size_mpi-1] = (ysize - (size_mpi-1)*num_rows_per_job)*xsize;
   
  for (i = 0; i < size_mpi; ++i)
  {
	  offsets[i] = (num_rows_per_job*i)*xsize;
  } 
  
  // Scatter
  MPI_Scatterv(src, counts, offsets, pixel_mpi, rcv_buf, counts[myid], pixel_mpi, 0, MPI_COMM_WORLD);
  
  for (y=0; y<counts[myid]/xsize; y++) { // For each row
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
      pix(prov_dst,x,y+radius, xsize)->r = r/n;
      pix(prov_dst,x,y+radius, xsize)->g = g/n;
      pix(prov_dst,x,y+radius, xsize)->b = b/n;
    }
  }
  
  // When applying the blur filter along y, we need to get the processed data from
  // other workers
  // Send top rows to previous section, and next section
  
  int num_blocks_to_send =  ceil((float)radius / (float)num_rows_per_job);

  for(i=0; i<num_blocks_to_send; i++){
    int num_elements_to_send;
    if (i == num_blocks_to_send -1) {
	  num_elements_to_send = xsize*(radius - i*num_rows_per_job);
    } else {
      num_elements_to_send = xsize*num_rows_per_job;
    }
 
    // Send to block above you
    if (myid > i) {
      MPI_Send(prov_dst+xsize*radius, num_elements_to_send, pixel_mpi, myid-1-i, 1000+i, MPI_COMM_WORLD);
    }

    int num_elements_to_rcv;
    
    // Receive from block below you
    if (myid < size_mpi-1-i) {
      MPI_Recv(prov_dst+counts[myid]+radius*xsize+i*num_rows_per_job*xsize, num_elements_to_send, pixel_mpi, myid+1+i, 1000+i, MPI_COMM_WORLD, &status);
    }
  }
 
  for(i=0; i<num_blocks_to_send; i++){
    int num_elements_to_send;
    if (i == num_blocks_to_send -1) {
	  num_elements_to_send = xsize*(radius - i*num_rows_per_job);
    } else {
      num_elements_to_send = xsize*num_rows_per_job;
    }
 
    // Send to block below you
    if (myid < size_mpi-1-i) {
      MPI_Send(prov_dst+xsize*radius + num_rows_per_job*xsize - num_elements_to_send, num_elements_to_send, pixel_mpi, myid+1+i, 2000+i, MPI_COMM_WORLD);
    }

    // Receive from block above you
    if (myid > i) {
      MPI_Recv(prov_dst+radius*xsize-i*num_rows_per_job*xsize-num_elements_to_send, num_elements_to_send, pixel_mpi, myid-1-i, 2000+i, MPI_COMM_WORLD, &status);
    }
  } 

  // Apply vertical blur
  for (y=radius; y<counts[myid]/xsize+radius; y++)  {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(prov_dst, x, y, xsize)->r;
      g = w[0] * pix(prov_dst, x, y, xsize)->g;
      b = w[0] * pix(prov_dst, x, y, xsize)->b;
      n = w[0];

      for ( wi=1; wi <= radius; wi++) {
        wc = w[wi];
        y2 = y - wi;
        if(y2 >=0) {
          r += wc * pix(prov_dst, x, y2, xsize)->r;
          g += wc * pix(prov_dst, x, y2, xsize)->g;
          b += wc * pix(prov_dst, x, y2, xsize)->b;
          n += wc;
        }
        y2 = y + wi;
        if(y2 < counts[myid]/xsize+2*radius) {
          r += wc * pix(prov_dst, x, y2, xsize)->r;
          g += wc * pix(prov_dst, x, y2, xsize)->g;
          b += wc * pix(prov_dst, x, y2, xsize)->b;
          n += wc;
        }
      }

      pix(rcv_buf,x,y-radius, xsize)->r = r/n;
      pix(rcv_buf,x,y-radius, xsize)->g = g/n;
      pix(rcv_buf,x,y-radius, xsize)->b = b/n;
    }
  }


  MPI_Gatherv(rcv_buf, counts[myid], pixel_mpi, src, counts, offsets, pixel_mpi, 0, MPI_COMM_WORLD);

  free(rcv_buf);
  free(prov_dst);
}



