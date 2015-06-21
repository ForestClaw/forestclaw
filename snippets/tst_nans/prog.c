#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double compute_slopes_(double *sl, double *sr);

static
void set_qnan(double* f)
{
    /*
     The quiet nan from math.h
    */
    *f = NAN;
}


static
void set_snan(double* f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)f) = 0x7ff0000000000001LL;
}

int main()
{
    double x,y,z;

    set_snan(&x);

    y = 1.0;

    z = compute_slopes_(&x,&y);
    printf("slope is %24.15f\n",z);

    z = compute_slopes_(&y,&x);
    printf("slope is %24.15f\n",z);

}
