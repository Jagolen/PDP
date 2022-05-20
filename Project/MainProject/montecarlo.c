/* 
Monte Carlo with the Stochastic Simulation Algorithm
By: Jakob Gölén
*/

#include "montecarlo.h"

int main(int argc, char **argv){
    if (3 != argc) {
        printf("Usage: num_experiments output_file\n");
        return 1;
	}

	//Parameters

	//Inputs
    int num_experiments = atoi(argv[1]);
    char *output_name = argv[2];

	//To store timings
	double local_time[4];
	double *times;
	double start_final, start_sub, end_final;

	//Given constants
	const int R = 15;
	const double T = 100;

	//Parameters used in the time loop
	double u1, u2, a0, tau, tt;
	int quarter_done, half_done, three_quarter_done, all_done;

	//State and result vectors
	double *w;
	int *x, *p;
	double *result;

	//For MPI
	int rank, size;

	//Initializing MPI
	MPI_Init(&argc, &argv);

	//Getting rank and size
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//number of experiments per process
	int n = num_experiments/size;
	if(num_experiments%size != 0){
		if(rank == 0) printf("Number of experiments (%d) is not divisible by number of processes (%d)\n", num_experiments, size);
		MPI_Abort(MPI_COMM_WORLD, -1);
		exit(-1);
	}

	//Initializing x, w and p
	x = (int*)malloc(7*sizeof(int));
	w = (double*)malloc(15*sizeof(double));
	p = (int*)calloc(15*7, sizeof(int));

	//the initial values of x was given, the values are initialized
	const int x0[7] = {900, 900, 30, 330, 50, 270, 20};

	//P is a constant and most values are zero which was initialized with calloc. The non zero elements are initialized
	p[0] = 1; p[7] = -1; p[14] = -1; p[16] = 1; p[22] = 1; p[29] = -1; p[36] = -1; p[38] = 1; p[44] = -1; p[51] = -1; p[53] = 1; 
	p[59] = -1; p[66] = -1; p[68] = 1; p[74] = -1; p[81] = -1; p[83] = 1; p[89] = -1; p[91] = 1; p[97] = -1; p[104] = -1;

	//To generate different numbers every time the program is ran and gives different random numbers on the different processes
	srand(time(0)+rank);


	//Main loop
	for(int iter = 0; iter<n; iter++){
		//Time is reset and x is set to initial value
		tt = 0;
		for(int i = 0; i<7; i++) x[i] = x0[i];
		quarter_done = 0;
		half_done = 0;
		three_quarter_done = 0;
		all_done = 0;

		//Time loop for every iteration
		while(tt<T){
			a0 = 0;
			prop(x,w);
			for(int i = 0; i<15;i++) a0 += w[i];
			u1 = ((double)(rand()%(RAND_MAX-1)+1))/(RAND_MAX);
			u2 = ((double)(rand()%(RAND_MAX-1)+1))/(RAND_MAX);
			tau = -log(u1/a0);
			tt += 10;
		}

	
	}
	double hi = MPI_Wtime();
	printf("u1 = %lf, u2 = %lf, tau = %lf\n", u1, u2, tau);

	//Freeing memory
	free(x);
	free(w);
	free(p);

	//Finalizing MPI and end the program
	MPI_Finalize();
    return 0;
}

//Prop function used to find w
void prop(int *x, double *w) {
	// Birth number, humans
	const double LAMBDA_H = 20;
	// Birth number, mosquitoes
	const double LAMBDA_M = 0.5;
	// Biting rate of mosquitoes
	const double B = 0.075;
	/* Probability that a bite by an infectious mosquito results in transmission
	   of disease to human*/
	const double BETA_H = 0.3;
	/* Probability that a bite results in transmission of parasite to a
	   susceptible mosquito*/
	const double BETA_M = 0.5;
	// Human mortality rate
	const double MU_H = 0.015;
	// Mosquito mortality rate
	const double MU_M = 0.02;
	// Disease induced death rate, humans
	const double DELTA_H = 0.05;
	// Disease induced death rate, mosquitoes
	const double DELTA_M = 0.15;
	// Rate of progression from exposed to infectious state, humans
	const double ALFA_H = 0.6;
	// Rate of progression from exposed to infectious state, mosquitoes
	const double ALFA_M = 0.6;
	// Recovery rate, humans
	const double R = 0.05;
	// Loss of immunity rate, humans
	const double OMEGA = 0.02;
	/* Proportion of an antibody produced by human in response to the incidence
	   of infection caused by mosquito. */
	const double NU_H = 0.5;
	/* Proportion of an antibody produced by mosquito in response to the
	   incidence of infection caused by human. */
	const double NU_M = 0.15;

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}