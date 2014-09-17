#include <iostream>
#include <cstdio>

extern "C"
{
    void mult(int *n, double* x);
    void assign(int* n, double x[]);
}

int main()
{
    int n = 10;
    double *x = new double[n];
    assign(&n,x);

    for (int i = 0; i < n; i++)
    {
        std::cout << x[i] << std::endl;
    }

    return 0;
}

void mult(int *n, double* x)
{
    for (int i = 0; i < *n; i++)
    {
        x[i] = 2*x[i];
    }
}
