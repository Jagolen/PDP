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

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  if(rank == 0)
    if(size <= 1) {
      printf("There needs to be at least 2 processes\n");
      MPI_Finalize();
      exit(0);
    }
  a = 100.0 + (double) rank;  /* Different a on different processors */

  /* Exchange variable a, notice the send-recv order */
  if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, 1, rank+1, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, size-1, rank, MPI_COMM_WORLD, &status);
    printf("Processor 0 got %f from processor %d\n", b,size-1);
  } else if (rank==size-1) {
    MPI_Recv(&b, 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD, &status);
    MPI_Send(&a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    printf("Processor %d got %f from processor %d\n",rank, b,rank-1);
  } else{
    MPI_Recv(&b, 1, MPI_DOUBLE, rank-1, rank, MPI_COMM_WORLD, &status);
    MPI_Send(&a, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD);
    printf("Processor %d got %f from processor %d\n",rank, b, rank-1);
  }

  MPI_Finalize(); 

  return 0;
}
