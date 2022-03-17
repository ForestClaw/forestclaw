#include <signal.h>      // for signal()                                         
#include <fenv.h>       // for fegetenv(), fesetenv()        
#include <math.h>        // for sqrt()                                           

#include <stdio.h>
#include <stdlib.h>

/* Needed for print_hexval, below */
#include <string.h>    /* Needed for memcpy */
#include <inttypes.h>  /* Needed for PRIx64 */

/* This doesn't seem to do anything */
void fpe_signal_handler(int sig) {
  printf("Floating point exception!\n");
  return;
}

void enable_floating_point_exceptions() 
{
 fenv_t env;
 fegetenv(&env);
 env.__fpcr = env.__fpcr | __fpcr_trap_invalid | __fpcr_trap_overflow 
     | __fpcr_trap_divbyzero | __fpcr_trap_underflow;
 fesetenv(&env);
 signal(SIGFPE, fpe_signal_handler);
}


static
void set_qnan(double *f)
{
    /*
     The quiet nan from math.h
    */
    *f = NAN;
}

static
void set_snan(double* f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)f) = 0x7ff0000000000001LL;
}

static
void print_hexval(double x)
{
    /* Kind of cool - didn't know it was possible to do this */
    uint64_t xn; 
    memcpy(&xn, &x, sizeof x);
    printf("x = %16f (%" PRIx64 ")\n", x, xn);
}

static
double create_exception(int choice)
{
    double x,y;
    switch(choice)
    {
        case 1:
        x = -1;
        y = sqrt(x);  // Or acos(x-1), log(x)
        print_hexval(y);
        printf("0: INVALID : y = sqrt(-1)\n");
        break;

        case 2:
        x = 1000;
        y = exp(x);
        print_hexval(y);
        printf("1: OVERFLOW : y = exp(1000)\n");
        break;

        case 3:
        y = pow(2,-1074);
        print_hexval(y);
        printf("1: UNDERFLOW : y = 2^(-1074)\n");
        break;

        case 4:
        x = 0;
        y = 1/x;
        print_hexval(y);
        printf("2: DIVBYZERO : y = 1/0\n");
        break;

        case 5:
        set_qnan(&x);
        print_hexval(x);
        y = x + 1;
        printf("3: Quiet NAN (not trapped) : x = NAN; y = x + 1;\n");
        break;

        case 6:
        set_snan(&x);
        print_hexval(x);
        y = x + 1;
        printf("3: Signaling NAN (trapped) : x = snan; y = x + 1;\n");
        break;

        default:
        y = -1;
        printf("Usage : tst_fpe <fpe_id>\n");
        printf("        fpe_id in [1-5]\n\n");
        exit(0);
    }
    return y;
}

int main(int argc, char** argv) 
{
    int fpe_id = 1;

    if (argc > 1)        
        fpe_id = atoi(argv[1]);

    printf("Exception handling disabled : \n");

    double y = create_exception(fpe_id);

    printf("Result : y = %g\n",y);

    printf("\n");
    printf("Enabling exception handling\n");  
    enable_floating_point_exceptions();

    y = create_exception(fpe_id);
    printf("Result : y = %g\n",y);
}

