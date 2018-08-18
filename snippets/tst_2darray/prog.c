#include <stdlib.h>

typedef struct array array_t;

struct array
{
    double *data[2];
};

void assign_array_(int* m, double* array, double* val);

int main()
{
    int n, k;
    double val;
    array_t *a;
    double *b;

    n = 8;
    val = 3.14159;

    a = (array_t*) malloc(sizeof(array_t));    
    for(k = 0; k < 2; k++)
    {
        a->data[k] = (double*) malloc(sizeof(double)*n);            
    }

    for(k = 0; k < 2; k++)
    {
        /* Can produce gdb "a=<error reading variable:..." error */
        assign_array_(&n,a->data[k],&val);  

        /* Doesn't seem to produce gdb  error */
        b = a->data[k];    
        assign_array_(&n,b,&val);
    }

    for(k = 0; k < 2; k++)
    {
        free(a->data[k]);        
    }
    free(a);    

    return 0;
}
