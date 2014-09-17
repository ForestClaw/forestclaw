#include <stdio.h>
#include <stdlib.h>

void mult(int *n, double* x);
void assign(int* n, double x[]);

int main()
{
    int n,i;
    double *x;

    n = 10;
    x = (double*) malloc(sizeof(double)*n);
    assign(&n,x);

    for (i = 0; i < n; i++)
    {
        printf("x is %4.0f\n",x[i]);
    }

    return 0;
}

void mult(int *n, double* x)
{
    int i;

    for (i = 0; i < *n; i++)
    {
        x[i] = 2*x[i];
    }
}
