#define _GNU_SOURCE 1
#include <fenv.h>

static void __attribute__ ((constructor))
trapfpe ()
{
    /* Enable some exceptions.  At startup all exceptions are masked.  */
     
    feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}