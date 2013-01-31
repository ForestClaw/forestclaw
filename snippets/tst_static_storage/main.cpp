#include <iostream>
#include <cstdlib>

#include "foo.H"

int main()
{
    foo f,g;

    f.print();
    g.print();

    // This call sets foo::x to the desired value.
    foo::x = 5;

    f.print();
    g.print();

}
