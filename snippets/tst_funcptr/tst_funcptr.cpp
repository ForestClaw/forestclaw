
#include <iostream>
#include <stdlib.h>

using namespace std;

typedef void (*set_value_t)(int& a);
typedef void (*math_t)(const int& a, int& b);

extern "C"
{
    // We can have the external definitions here, even if the file is not
    // compiled in.
    void set_value_to_2_(int& a);
    void set_value_to_5_(int& a);
    void math_value_(math_t f, const int& a, int& b);
    void square(const int& a, int& b);
    void cube(const int& a, int& b);
}

int main ()
{
    int a;
    set_value_t f_set;
    math_t f_math;

    // Only these files need to be compiled into executable.  See Makefile.inc
    f_set = &set_value_to_5_; // Fortran function called directly from C
    f_math = &square; // C function that is passed to and called from Fortran

    f_set(a);
    cout << "a = " << a << endl;

    int b;
    math_value_(f_math,a,b);
    cout << "f(a) = " << b << endl;
}
