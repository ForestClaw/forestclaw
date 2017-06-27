#include <stdio.h>
#include <stdlib.h>

void assign_array_();
void print_array_();
void set_ptr_(double **x, int* n);

int main()
{
    int n,i;
    double *concen;

    n = 5;
    concen = (double*) malloc(sizeof(double)*n*n);

    set_ptr_(&concen,&n);
    assign_array_();
    print_array_();

    free(concen);

    return 0;
}
