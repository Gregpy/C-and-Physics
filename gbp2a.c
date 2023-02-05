// compile with -lm for math library link

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXORD  8 // the highest order system (largest # if dependent variables
#define MAXPAR 16 // largest number of parameters

char *prog; // the will be set to the program name. it is a global identifier and
            // as such is known to ALL functions

int
main(int argc, char *argv[])
{
	double x0, v0, p[MAXPAR], tmax, t, dt, xe, ve, r[MAXORD], w, x_err;
	void rk4(int, double*, double, double, double*);
	int n, nsteps;
	FILE *fp;

	prog = argv[0];
/*
 * open input file and read input data:
 */
	if ((fp=fopen("gbp2a.inp", "r")) == NULL) {
		fprintf(stderr, "%s: Cant't open gbp2a.inp for reading ", prog);
		perror("because");
		return 2;
	}

	printf("Using the 4th-order Runge-Kutta (RK4) method to solve the simple harmonic oscillator.\n");

	printf("Enter initial conditions 'x0' and 'v0':\n");

	if (scanf("%le %le", &x0, &v0) !=2) {
		printf("Bad input\n");
		return 1;
	}

	printf("Enter the maximum time 'tmax' and total number of time steps 'M':\n");

	if (scanf("%le %d", &tmax, &nsteps) !=2) {
	        printf("Bad input\n");
       		return 1;
	}

	if (fscanf(fp, "%le", &w) != 1) {
		fprintf(stderr, "%s: bad input on line 1 of gbp2a.inp.\n", prog);
		return 1;
	}

	fclose(fp);

/*
 * need to put stuff into arrays for passing to rk4() b/c rk4()
 * works for ANY # of dependent variables, at least up to MAXORD
 */
	r[0] = x0; // use r[0] for positions
	r[1] = v0; // use r[1] for velocities
	p[0] = w;  // use p[0] for frequency

	dt = tmax/nsteps;
	x_err = 0.0;
	t = 0.0;

	if ((fp=fopen("gbp2a.txt", "w")) == NULL) {
		fprintf(stderr, "%s: Can't open gbp2a.txt for writing ", prog);
		perror("because");
		return 2;
	}
/*
 * want initial condition in input file so write them here before looping
 */
	fprintf(fp, "%le %le %le %le %le\n", t, r[0], x0, r[1], v0);

	for (n = 0; n < nsteps; n++) { // N.B.: start counting at n=0 so that ...

		t = n*dt; // ... "t" is the time at the BEGINNING of the interval!

		/*
		 * advance solution from t --> t+dt:
		 */
		rk4(2, r, t, dt, p);

		/*
		 * NOW update t to END of interval, compute exact results:
		 */
		t += dt;
		xe =    x0*cos(w*t) + v0*sin(w*t)/w;
		ve = -w*x0*sin(w*t) + v0*cos(w*t);

		/*
		 * do trapeziodal rule integration of [Xem(t)-Xex(t)]^2
		 * accumulate the running sum of the integrand
		 * !SUBTLEY! : Xem(0)-Xex(0)=0 BY DEFINITION!
		 */
		if (n == nsteps-1) // last time around n=nsteps-1, NOT nsteps
			x_err += pow(r[0]-xe, 2.0)/2; // last step weight = 1/2 !!!!
		else
			x_err += pow(r[0]-xe, 2.0);

		fprintf(fp, "%le %le %le %le %le\n", t, r[0], xe, r[1], ve);
	}

	fclose(fp);

	x_err *= dt; // multiply by width of trapezoid
	x_err = sqrt(x_err/tmax); // compute RMS average from integral

	printf("Global position error = %le\n", x_err);

	printf("Data for 4th-order Runge-Kutta calculated position 'xRK', exact position 'xExact', 4th-order Runge-Kutta calculated velocity 'vRK', and exact velocity 'vExact' versus time 't' is sent to the file gbp2a.txt.\n");

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

	drdt[0] = r[1];
	drdt[1] = -p[0]*p[0]*r[0];
        drdt[2] = r[1] + dt/2*drdt[1];
        drdt[3] = -p[0]*p[0]*(r[0] + dt/2*drdt[0]);
	drdt[4] = r[1] + dt/2*drdt[3];
        drdt[5] = -p[0]*p[0]*(r[0] + dt/2*drdt[2]);
        drdt[6] = r[1] + dt*drdt[5];
        drdt[7] = -p[0]*p[0]*(r[0] + dt*drdt[4]);

	return;
}
