#include <iostream>
using namespace std;

typedef void (*set_value_t)(int& a);
typedef void (*square_value_t)(int& a);

extern "C"
{
    void set_to_2_(int& a);
    void set_to_5_(int& a);
    void fortran_cannot_square_(square_value_t f, int& a);
}

void square(int& a)
{
    a *= a;
}


int main ()
{
    int a;
    int ichoice;
    set_value_t f;

    cout << "Input choice (0 or 1) : " << endl;
    cin >> ichoice;

    if (ichoice == 0)
    {
        f = &set_to_2_;
    }
    else if (ichoice == 1)
    {
        f = &set_to_5_;
    }
    else
    {
        cout << "(:-((" << endl;
        exit(1);
    }

    f(a);
    cout << "Success! a = " << a << endl;

    fortran_cannot_square_(square,a);
    cout << "Squared value : a = " << a << endl;
}
