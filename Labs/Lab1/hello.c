/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int size,rank;
  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  printf("Hello World! from process %d\n",rank);             /* Print a message              */
  if(rank == 0)
    printf("Process %d say there are %d processes\n",rank,size);

  MPI_Finalize();                       /* Shut down and clean up MPI   */

  return 0;
}
