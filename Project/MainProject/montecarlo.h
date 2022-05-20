#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * Compute propensities for the Malaria model.
 * @param x State vector Should be of length 7!
 * @param w Result vector (propensities). Should be of length 15!
 *
 */
void prop(int *x, double *w);