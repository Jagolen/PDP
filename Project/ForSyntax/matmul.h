#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int read_input(const char *file_name, double **A, double **B);

int write_output(const char *file_name, double *C, int matrix_size);