// compile with -lm for math library link

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXORD  8 // the highest order system (largest # if dependent variables
#define MAXPAR 16 // largest number of parameters

char *prog; // the will be set to the program name. it is a global identifier and
            // as such is known to ALL functions

int
main(int argc, char *argv[]) // main() with arguments!
{
	double nr0, nf0, p[MAXPAR], tmax, t, dt, r[MAXORD], rho, x_err;
	void rk4(int, double*, double, double, double*);
	int n, nsteps;
	FILE *fp;

	prog = argv[0];

	printf("Using the 4th-order Runge-Kutta (RK4) method with the fox-rabbit system of equations.\n");

	printf("Enter initial number of rabbits and foxes 'nr0' and 'nf0':\n");

	if (scanf("%le %le", &nr0, &nf0) !=2) {
		printf("Bad input, program will terminate.\n");
		return 1;
	}

	printf("Enter the ratio of the starvation rate of foxes to the growth rate of rabbits rR/rF 'rho':\n");

	if (scanf("%le", &rho) != 1) {
		printf("Bad input, program will terminate.\n");
		return 1;
	}

	printf("Enter the maximum time 'tmax' and total number of time steps 'M':\n");

	if (scanf("%le %d", &tmax, &nsteps) !=2) {
	        printf("Bad input, program will terminate.\n");
       		return 1;
	}


/*
 * need to put stuff into arrays for passing to rk4() b/c rk4()
 * works for ANY # of dependent variables, at least up to MAXORD
 */
	r[0] = nr0; // use r[0] for rabbit numbers
	r[1] = nf0; // use r[1] for fox numbers
	p[0] = rho;  // use p[0] for rho

	dt = tmax/nsteps;
	x_err = 0.0;
	t = 0.0;

	if ((fp=fopen("gbp2b.txt", "w")) == NULL) {
		fprintf(stderr, "%s: Can't open gbp2b.txt for writing ", prog);
		perror("because");
		return 2;
	}
/*
 * want initial condition in input file so write them here before looping
 */
	fprintf(fp, "%le %le %le\n", t, r[0], r[1]);

	for (n = 0; n < nsteps; n++) { // N.B.: start counting at n=0 so that ...

		t = n*dt; // ... "t" is the time at the BEGINNING of the interval!

		/*
		 * advance solution from t --> t+dt:
		 */
		rk4(2, r, t, dt, p);

		t += dt;

		fprintf(fp, "%le %le %le\n", t, r[0], r[1]);
	}

	fclose(fp);

	printf("Rabbits and foxes data 'nR(t)' and 'nF(t)' versus time 't' is sent to the file gbp2b.txt\n");

	return 0;
}

void
rk4(int ord, double *r, double t, double dt, double *p)
{
	double drdt[MAXORD]; // array for RHS derivatives
	int i;
	void model(int, double*, double*, double, double, double*);

	model(ord, r, drdt, t, dt, p);

	for (i = 0; i < ord; i++){
	  r[i] += dt/6*(drdt[i] + 2*drdt[i+2] + 2*drdt[i+4] + drdt[i+6]);

}

	return;
}

void
model(int ord, double *r, double *drdt, double t, double dt, double *p)
{
/*
 * SHO:
 */
   int i;

	drdt[0] = r[0]*(1-r[1]);
	drdt[1] = p[0]*r[1]*(r[0]-1);
        drdt[2] = (r[0] + dt/2*drdt[1])*(1-(r[1]+dt/2*drdt[1]));
        drdt[3] = p[0]*(r[1]+dt/2*drdt[0])*(r[0]+dt/2*drdt[0]-1);
	drdt[4] = (r[0] + dt/2*drdt[3])*(1-(r[1]+dt/2*drdt[3]));
        drdt[5] = p[0]*(r[1]+dt/2*drdt[2])*(r[0]+dt/2*drdt[2]-1);
        drdt[6] = (r[0] + dt*drdt[5])*(1-(r[1]+dt*drdt[5]));
        drdt[7] = p[0]*(r[1]+dt*drdt[4])*(r[0]+dt*drdt[4]-1);

	return;
}
