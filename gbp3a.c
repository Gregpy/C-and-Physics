// compile with -lm for math library link

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int

main(int argc, char *argv[])
{

	double *p, *l, *d, *u, *r, x, dx, xL, xR, lamk, lamp;

	double exact(double x, double lamk, double lamp);

	int n, i;

	int td_solve(double *l, double *d, double *u, double *r, double *y, int n);

	FILE *fp;

	printf("Using FEBS to compute a two-point boundary value problem, and comparing it to the exact solution, part 1.\n");

	printf("Enter lambda-k and lambda-p: ");

	if(scanf("%le %le", &lamk, &lamp) != 2){

		printf("Error, program will terminate.\n");

		exit(2);
	}

	printf("Enter N: ");

	if(scanf("%d", &n) != 1)
	{
		printf("Error, program will terminate.\n");

		exit(2);
	}

	p = (double*)malloc(n*sizeof(double));

	l = (double*)malloc(n*sizeof(double));

	d = (double*)malloc(n*sizeof(double));

	u = (double*)malloc(n*sizeof(double));

	r = (double*)malloc(n*sizeof(double));

	xL = -0.5;

	xR = 0.5;

	dx = (xR - xL)/(n-1);

	for (i = 0; i < n; i++) {

		u[i] = 1/(dx*dx);

		l[i] = 1/(dx*dx);

		d[i] = -2/(dx*dx) - 1.0/(lamk*lamk);

		r[i] = 0.0;
}

	u[0] += l[0];

	d[0] -= 2*dx*l[0]/lamp;

	r[0] -= 2*dx*l[0]/lamp;

	l[n-1] += u[n-1];

	d[n-1] -= 2*dx*u[n-1]/lamp;

	r[n-1] -= 2*dx*u[n-1]/lamp;


	if (td_solve(&l[0], &d[0], &u[0], &r[0], &p[0], n) != 0) {

		fprintf(stderr, "td_solve() failed!\n");

		exit(2);
}
	fp = fopen("gbp3a.txt", "w");
	printf("Results sent to gbp3a.txt\n");
	for (i = 0; i < n; i++) {

		x = xL + i*dx;

		fprintf(fp, "%le %le %le\n", x, p[i], exact(x, lamk, lamp));
}
	fclose(fp);

	exit(0);
}

	double
	exact(double x, double lamk, double lamp)
{

	return cosh(x/lamk)/(cosh(1/(2*lamk))+(lamp/lamk)*sinh(1/(2*lamk)));
}

	int
	td_solve(double *l, double *d, double *u, double *r, double *y, int n)
{

	int i;
/*
* do Gaussian elimination:
*/

	for (i = 1; i < n; i++) {

		if (d[i-1] == 0.0) return 1;

		d[i] -= l[i]*u[i-1]/d[i-1];

		r[i] -= l[i]*r[i-1]/d[i-1];
}
/*
* do back-substitution:
*/
	if (d[n-1] == 0.0) return 2;

	y[n-1] = r[n-1]/d[n-1];

	for (i = n-2; i >= 0; i--) {

		if(d[i] == 0.0) return 2;

		y[i] = (r[i] - u[i]*y[i+1])/d[i];
}
	return 0;
}

