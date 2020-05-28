//assignment3.c

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void read_input(double *C, double *L, int *nx, double *t_F,double *t_out);

int main(void) {
  // **********
  // Parameters
  // **********

	//Wrote to file to allow for a speedy graphical transition without using command prompt.
	FILE *file2;
	file2= fopen("outputAssignment3.txt", "w");
	  if (file2 == NULL) {
      return 1;
	  }

  // Number of grid points
  int *nx_p;
  nx_p= (int*) malloc(sizeof(int));


  // Length of domain
  double *L_p;
  L_p= (double*) malloc(sizeof(double));


  // Equation coefficients
  double *C_p;
  C_p= (double*) malloc(sizeof(double));


  // Length of time to run simulation.
  double *t_F_p;
  t_F_p= (double*) malloc(sizeof(double));


  // How frequently in time to output.
  double *output_timestep_p;
  output_timestep_p= (double*) malloc(sizeof(double));


  // Read in from file;
  read_input(C_p, L_p, nx_p, t_F_p, output_timestep_p);


  // Copy values from pointers to variables.
  int nx = *nx_p;
  double L = *L_p;
  double C = *C_p;
  double t_F = *t_F_p;
  double output_timestep = *output_timestep_p;


  // Grid spacing, now fixed
  double dx = L / (nx - 1);


  // Time step
  double dt;


  // Set a timestep so the code is stable
  if(C == 0.0){
    dt = 0.001;
  }
  else{
  dt = dx / (2 * C);
  } //means we don't get a stupidly huge timestep at c=0, but gives us a stable timestep in both cases.
  // ************
  // Grid Storage
  // ************
  double *U;
  double *U_next;  //x at current and next timestep


  double *V;
  double *V_next;  //y at current and next timestep

  /* Allocate memory according to size of nx */
  U       = (double*) malloc(nx * sizeof(double));
  U_next  = (double*) malloc(nx * sizeof(double));
  V       = (double*) malloc(nx * sizeof(double));
  V_next  = (double*) malloc(nx * sizeof(double));

//Ensures we do in fact assign memory appropriately before carrying on with the program, and terminating if we do not.
  if(U == NULL || U_next == NULL || V == NULL || V_next == NULL){
	printf("Memory allocation failed \n");
	return 1;
  }

  int j; //Index for current gridpoint
  double x; //Denotes the x-value.
  double ctime;
  ctime = 0.0; // Current time.

  // **************
  // initialisation
  // **************
  for(j = 0;j < nx; ++j) {
    U[j]  = 1.0; //Initialises U on all timesteps to 1.
    V[j]  = 0.0; //Initialises V on all timesteps to 0.
  }

  double next_output_time = 0.0; //Initialises the timestep we will first look at, then acts as the next timestep in subsequent iterations - misleading title on first iteration, but then okay for subsequent iterations.

  //loop over timesteps

  while (ctime < t_F) {
		double dt0 = dt;
		int output = 0;
		// If we would go past the next output step, reduce the timestep.
		if (ctime + dt0 > next_output_time) {
			dt0 = next_output_time - ctime;
			output = 1;
		}

		//Rotation factors for time-splitting scheme.
		double cfac = cos(dt0);
		double sfac = sin(dt0);

		//loop over points: first rotation
		for (j = 0; j < nx; ++j) {
			// Apply rotation matrix to U and V (in other words, solve
			// the differential equation without the spatial derivatives).
                U[j] = cfac * U[j] + sfac * V[j];
                V[j] = - sfac * U[j] + cfac * V[j];
		}

		//loop over points: then advection
		for (j = 1; j < nx; ++j) {
			signed int jp = j-1;

			// Need upwinding for stability
			double dUdx = (U[j] - U[jp]) / dx;
			double dVdx = (V[j] - V[jp]) / dx;

			U_next[j] = U[j] + V[j] * dt0 - dt0 * C * dUdx; //adds a damping term as per the lecture notes.
			V_next[j] = V[j] - U[j] * dt0 - dt0 * C * dVdx; // ""  ""   ""     "" "" ""   ""    ""    "".
		}

		// Set boundary values
    U_next[0] = 0.0;
    V_next[0] = 0.0;

		// Copy next values at timestep to U,V arrays. Now actually does this element by element.
    for(j = 0; j < nx; ++j){
		  U[j] = U_next[j];
		  V[j] = V_next[j];
    }
		//Now increment our time up by dt0.
		ctime += dt0;
		if (output) {
			for (j = 0; j < nx; ++j) {
				x = j * dx;
				printf("%g %g %g %g\n", ctime, x, U[j], V[j]);
				fprintf(file2, "%g %g %g %g\n", ctime, x, U[j], V[j]);
            }
			next_output_time += output_timestep;
        }
        else{
        }
	}
  free(U);
  free(U_next);
  free(V);
  free(V_next);
  free(nx_p);
  free(L_p);
  free(C_p);
  free(t_F_p);
  free(output_timestep_p); //frees everything properly - no pointers left in the code unfreed.
  fclose(file2); //closes file used for data transfer.
  return 0; //should exit with 0 if we're lucky!
}

// The lines below don't contain any bugs! Don't modify them
void read_input(double *C, double *L, int *nx, double *t_F, double *t_out) {
   FILE *infile;
   if(!(infile=fopen("inputAssignment3.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(5!=fscanf(infile,"%lf %lf %d %lf %lf", C, L, nx, t_F, t_out)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
