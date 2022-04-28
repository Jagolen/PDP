#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    int left,right,rank,size, prev, next;
	// Check arguments
	if (3 != argc) {
		printf("Usage: ./a.out input1 input2, where input1 is the number of rows and input2 is the number of columns in the local array.\n");
		return -3;
	}
    MPI_Init(&argc, &argv);
	char *nt = argv[1];
	char *mt = argv[2];

    int n = atoi(nt);
    int m = atoi(mt);

    MPI_Status status;



    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    next = (rank + 1) % size;
    prev = (size + rank - 1) % size;

    int *matrix = (int*)malloc(n*m*sizeof(int));

    for(int i = 0; i < m*n; i++){
        matrix[i] =  rank * m * n + i;
    }
    if(rank==1){
        printf("BEFORE: RANK %d, matrix \n", rank);
        for(int i = 0; i<n; i++){
            for(int j = 0; j<m; j++){
                printf("%d \t", matrix[i*n+j]);
            }
            printf("\n");
        }
    }
    MPI_Datatype contig;
    MPI_Type_contiguous(2*m, MPI_INT, &contig);
    MPI_Type_commit(&contig);

    if(rank == 0){
        printf(" %d sending to %d with tag %d\n", rank, next, next);
        MPI_Send(&matrix[n*m - 2 * m], 1, contig, next, next, MPI_COMM_WORLD);
        MPI_Recv(&matrix[n*m - 2 * m - 1], 1, contig, prev, rank, MPI_COMM_WORLD, &status);
    }
    else{
    printf("rank %d recieving from, %d with tag %d\n", rank, prev, rank);
    MPI_Recv(&matrix[n*m - 2 * m], 1, contig, prev, rank, MPI_COMM_WORLD, &status);
    MPI_Send(&matrix[n*m - 2 * m], 1, contig, next, next, MPI_COMM_WORLD);
    }


    if(rank==1){
        printf("AFTER: RANK %d, matrix \n", rank);
        for(int i = 0; i<n; i++){
            for(int j = 0; j<m; j++){
                printf("%d \t", matrix[i*n+j]);
            }
            printf("\n");
        }
    }




    MPI_Type_free(&contig);
    MPI_Finalize();
    return 0;
}