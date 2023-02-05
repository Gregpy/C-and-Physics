// compile with -lm for math library link

/*The values the variables correspond to are: 'n' is the number of discretization points; 'i' is the counter used in the for loop; 'a' is the lower limit of the integral; 'b' is the upper limit of the integral; 'd' is dx or delta x, the width of the trapezoids; 's' is the sum of the area of the trapezoids. The sum of the first and last functions are calculated separately from the rest since they are divided by 2.*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double fun(double x);

int main()
{

  	int n,i;
	double a,b,d,s=0;

	printf("Using the Trapezoidal Rule to determine the definite integral of 1+x^2+x^3.\n");

	printf("Enter the lower limit 'a':\n");
	scanf("%le,\n",&a);

	printf("Enter the upper limit 'b':\n");
	scanf("%le,\n",&b);

	printf("Enter number of discretization points 'n', greater than 1:\n");
	scanf("%d",&n);

	if(n<=1){

	printf("Invalid input, program will terminate.\n");
	exit(0);

	}

	d=(b-a)/(n-1);

	s=d/2*(fun(a) + fun(b));

	for(i=1;i<=n-2;i++){

	s=s+d*fun(a+i*d);

	}

	printf("The signed area under the curve is equal to: %.6lf\n",s);

return (0);
}

double fun(double x)

{

double form;

	form = 1 + x*x + x*x*x;

return (form);

}




