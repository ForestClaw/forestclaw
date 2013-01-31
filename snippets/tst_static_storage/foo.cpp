#include "foo.H"
#include <iostream>
#include <cstdlib>

// This acts as a default value (??)
double foo::x;

foo::foo()
{
    foo::x = 1;
}

void foo::print()
{
    printf("x = %f\n",foo::x);
}
