//Assignment4.c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void read_input(double *L_p, int *N_p, double *t_F_p, double *t_D_p, double *K_p, double *A_p, double *B_p);

double max(double A, double B, double K);

int main(void) {
  // **********
  // Parameters
  // **********

	//Wrote to file to allow for a speedy graphical transition without using command prompt.
	FILE *outfile;
	outfile= fopen("outputAssignment4.txt", "w");
	if (outfile == NULL){
	    exit(1);
    }

  // Length of domain
  double L;

  // Number of grid points
  int N;

  // Initial condition parameters.
  double K;
  double A;
  double B;


  // Length of time to run simulation.
  double t_F;

  // Timestep for diagnostic output.
  double t_D;

  // Read values in from file.
  read_input(&L, &N, &t_F, &t_D, &K, &A, &B);

  // Grid spacing
  double dx = L / N;

  //Our simulation time step - Calculated for the standard heat equation.
  double dt;
  double maxCoeff=max(A, B, K);
  if(A > 1.0 || B > 1.0 || K > 1.0) {
    dt=(dx * dx) / (4 * maxCoeff * maxCoeff * maxCoeff);
  }
  else{
    dt = (dx * dx) / 4;
  }

  int j; //Declare index for current gridpoint globally as opposed to at the start of every for loop.
  double x; //Represents our x-value on output.
  double ctime = 0.0; //Starts simulation at t=0.

  //Allocates memory for our function variables U and V.
  double *U;
  U= (double*) malloc((N) * sizeof(double));
  double *V;
  V= (double*) malloc((N) * sizeof(double));
  double *U_next;
  U_next= (double*) malloc((N) * sizeof(double));
  double *V_next;
  V_next= (double*) malloc((N) * sizeof(double));

  //Checks to see if memory was correctly allocated for U and V, and exiting if not.
  if(U == NULL || V == NULL || U_next == NULL || V_next == NULL) {
        printf("Memory allocation didn't work!");
        exit(1);
  }

  //Initialising U and V as per the brief
  for(j = 0; j < N; ++j) {
      x= j * dx;
      U[j] = K + A * cos((2 * M_PI * x) / L);
      V[j]=B * sin((2 * M_PI * x) / L);
  }

  double output_time = 0.0;
  //Works in a similar manner as before
  while (ctime < t_F){
    double dt_0 = dt;
    int output = 0;
    if(ctime + dt_0 > output_time) {
        dt_0 = output_time - ctime;
        output = 1;
    }
    //Uses numerical method to construct an approximation to the next timestep
    for(j = 1; j < N - 1; ++j){
        double deriv_U = (U[j - 1] + U[j + 1] - 2 * U[j]) / (dx * dx);
        double deriv_V = (V[j - 1] + V[j + 1] - 2 * V[j]) / (dx * dx);
        U_next[j] = U[j] + deriv_U * dt_0 + 2 * dt_0 * U[j] - 4 * U[j] * dt_0 * (U[j] * U[j] + V[j] * V[j]);
        V_next[j] = V[j] + deriv_V * dt_0 + 2 * dt_0 * V[j] - 4 * V[j] * dt_0 * (U[j] * U[j] + V[j] * V[j]);
    }

    //Set our boundary conditions

    //For the left boundary:
    double deriv_U = (U[N - 1] + U[1] - 2 * U[0]) / (dx * dx); //uses the fact that U,V are periodic
    double deriv_V = (V[N - 1] + V[1] - 2 * V[0]) / (dx * dx);
    U_next[0] = U[0] + deriv_U * dt_0 + 2 * dt_0 * U[0] - 4 * U[0] * dt_0 * (U[0] * U[0] + V[0] * V[0]);
    V_next[0] = V[0] + deriv_V * dt_0 + 2 * dt_0 * V[0] - 4 * V[0] * dt_0 * (U[0] * U[0] + V[0] * V[0]);

    //And now for the right boundary, again using the fact that we have periodicity:
    deriv_U = (U[N - 2] + U[0] - 2 * U[N - 1]) / (dx * dx);
    deriv_V = (V[N - 2] + V[0] - 2 * V[N - 1]) / (dx * dx);
    U_next[N - 1] = U[N - 1] + deriv_U * dt_0 + 2 * dt_0 * U[N - 1] - 4 * U[N - 1] * dt_0 * (U[N - 1] * U[N - 1] + V[N - 1] * V[N - 1]);
    V_next[N - 1] = V[N - 1] + deriv_V * dt_0 + 2 * dt_0 * V[N - 1] - 4 * V[N - 1] * dt_0 * (U[N - 1] * U[N - 1] + V[N - 1] * V[N - 1]);


    //Copy next timestep info into current timestep ready for next iteration.
    for(j = 0; j < N; ++j) {
        U[j] = U_next[j];
        V[j] = V_next[j];
    }
    ctime += dt_0;

    //Outputs values calculated
    if(output) {
        for (j = 0; j < N; ++j){ 
                double x = j * dx;
                printf("%g %g %g %g\n", ctime, x, U[j], V[j]);          // Added a print to terminal statement so the user can check the output doesn't blow up while the program runs.
                fprintf(outfile, "%g %g %g %g\n", ctime, x, U[j], V[j]); 
        }
        output_time += t_D;
    }
  }
  //Freeing pointers now we're finished with them, and closing our output file
  //fclose(outfile);
  free(U);
  free(V);
  free(U_next);
  free(V_next);
  return 0;
}

void read_input(double *L_p, int *N_p, double *t_F_p, double *t_D_p, double *K_p, double *A_p, double *B_p) {
   FILE *infile;
   if(!(infile=fopen("inputAssignment4.txt", "r"))) {
       printf("Error opening file.\n");
       exit(1);
   }
   else{
   }
   if(7 != fscanf(infile,"%lf %d %lf %lf %lf %lf %lf", L_p, N_p, t_F_p, t_D_p, K_p, A_p, B_p)) {
       printf("Error reading parameters from file.\n");
       exit(1);
   }
   //added in order to make sure we get a file that has no more and no less than 7 inputs
   if(fscanf(infile, "%lf %d %lf %lf %lf %lf %lf", L_p, N_p, t_F_p, t_D_p, K_p, A_p, B_p) > 7) {
       printf("Too many inputs in this file.\n");
       exit(1);
   }
   else{
   }

   fclose(infile);
}

double max(double A, double B, double K) {
  double max = A;
  if(max < B) {
      max = B;
  }
  if(max < K) {
     max = K;
  }
  return max;
}
