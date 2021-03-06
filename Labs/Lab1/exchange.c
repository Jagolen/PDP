/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  MPI_Status status;
  MPI_Request send,rec;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  int prev = rank-1;
  int next = rank+1;
  if(next > size-1) next = 0;
  if(prev < 0) prev = size-1;

  if(rank == 0)
    if(size <= 1) {
      printf("There needs to be at least 2 processes\n");
      MPI_Finalize();
      exit(0);
    }
  a = 100.0 + (double) rank;  /* Different a on different processors */

/* Exchange variable a, notice the send-recv order */
  MPI_Irecv(&b, 1, MPI_DOUBLE, prev, rank, MPI_COMM_WORLD, &rec);
  MPI_Isend(&a, 1, MPI_DOUBLE, next, next, MPI_COMM_WORLD, &send);
  MPI_Wait(&rec, &status);
  printf("Processor %d got %f from processor %d\n",rank, b, prev);

  MPI_Finalize(); 

  return 0;
}
