#include <iostream>
#include <cstdio>

typedef struct simple
{
    int multiplier;
} simple_t;

extern "C"
{
    void mult(simple_t *s, const int& n, double x[]);
    void f_assign(simple_t *s, const int& n, double x[]);
}

int main()
{
    int n = 10;
    double *x = new double[n];
    simple_t s, *s_ptr = &s;
    s.multiplier = 3;
    f_assign(s_ptr,n,x);

    for (int i = 0; i < n; i++)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}

void mult(simple_t *s, const int& n, double x[])
{
    for (int i = 0; i < n; i++)
    {
        x[i] = s->multiplier*x[i];
    }
}
