#include "thresfilter.h"
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
void thresfilter(const int xsize, const int ysize, pixel* src){
#define uint unsigned int 

  int myid;
  int size_mpi;

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
  MPI_Status status;


  // The image is split horizontally into blocks, and each worker works
  // on one block
  int num_rows_per_job = ceil((float)ysize/(float)size_mpi);
  pixel* rcv_buf = (pixel*)calloc(num_rows_per_job*xsize, sizeof(pixel));
  
  int counts[size_mpi];
  int offsets[size_mpi];

  uint i;
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
  
  uint sum, psum, nump;

  nump = xsize * ysize;

  // Calculate the sum
  for(i = 0, sum = 0; i < counts[myid]; i++) {
    sum += (uint)rcv_buf[i].r + (uint)rcv_buf[i].g + (uint)rcv_buf[i].b;
  }

  // Reduce
  uint sum_all;
  MPI_Reduce(&sum, &sum_all, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  // Calculate average
  sum = sum_all;
  if (myid==0) {
      sum /= nump;
  }

  // Broadcast average to all workers
  MPI_Bcast(&sum, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  // Threshold
  for(i = 0; i < counts[myid]; i++) {
    psum = (uint)rcv_buf[i].r + (uint)rcv_buf[i].g + (uint)rcv_buf[i].b;
    if(sum > psum) {
      rcv_buf[i].r = rcv_buf[i].g = rcv_buf[i].b = 0;
    }
    else {
      rcv_buf[i].r = rcv_buf[i].g = rcv_buf[i].b = 255;
    }
  }

  MPI_Gatherv(rcv_buf, counts[myid], pixel_mpi, src, counts, offsets, pixel_mpi, 0, MPI_COMM_WORLD);
  free(rcv_buf);
}
