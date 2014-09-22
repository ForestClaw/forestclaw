#include <stdio.h>
#include <stdlib.h>

typedef struct simple
{
    int multiplier;
} simple_t;


void mult(simple_t *s, int *n, double x[]);
void assign(simple_t *s, int* n, double x[]);

int main()
{
    int n,i;
    double *x;
    simple_t s;

    n = 3;
    s.multiplier = 23;
    x = (double*) malloc(sizeof(double)*n);
    assign(&s,&n,x);

    for (i = 0; i < n; i++)
    {
        printf("x(%d) is %4.0f\n",i,x[i]);
    }

    return 0;
}

void mult(simple_t *s, int *n, double* x)
{
    int i;

    for (i = 0; i < *n; i++)
    {
        x[i] = s->multiplier*x[i];
    }
}
