#include <stdio.h>
#include <stdlib.h>

void assign_array_(double x[], int* n);
void print_array_();

int main()
{
    int n,i;
    double *x;

    n = 5;
    x = (double*) malloc(sizeof(double)*n);
    if (x == NULL)
    {
        printf("Malloc did not success\n");
        exit(0);
    }
    else
    {
        printf("Malloc succeeded\n");
    }

    for (i = 0; i < n; i++)
        x[i] = (double) i;

    assign_array_(x,&n);
    print_array_();

    return 0;
}
