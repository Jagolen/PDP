/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a;
  MPI_Status status;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  int maxsend = 0;
  int sender = 0;
  int reciever = 1;
  int tag = rank;
  /* Processor 0 send to all others */
  if (rank == 0) {
    a=999.999;
  }

  while((int)pow(2, maxsend) < size){
    if(rank < (int)pow(2, maxsend)){
      printf("process %d sending to process %d\n", rank, reciever);
      MPI_Send(&a, 1, MPI_DOUBLE, reciever,reciever,MPI_COMM_WORLD);
    }
    else if(rank < pow(2, maxsend+1) && rank > maxsend+1 ) {
      MPI_Recv(&a, 1, MPI_DOUBLE, sender, rank, MPI_COMM_WORLD, &status);
      printf("Process %d recieving %f from %d\n", rank, a, sender);
    }
    maxsend += 1;
    reciever = (int)pow(2, maxsend) + rank;
    sender = rank % (int)pow(2, maxsend);
  }
  MPI_Finalize(); 

  return 0;
}
