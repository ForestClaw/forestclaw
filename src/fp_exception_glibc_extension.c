/*   Title:  Floating-point exception handling example
    Author:  David N. Williams
      File:  fe-handlng-example.c
   License:  Public Domain
   Version:  0.5.0
   Started:  21-Sep-09
   Revised:  22-Sep-09
   Revised:  30-Sep-09 (comment typo)
   Revised:  18 Oct-12 (chnaged char* to const char * on line 228, by Richard Booth)

This code is an example of alternate, nondefault handling of
IEEE 754 floating-point exceptions in OS X and Linux, based on
the GNU functions feenableexcept(), fedisableeexcept(), and
fegetexcept() [in libm], plus POSIX sigaction().

The GNU functions above are not implemented in OS X Leopard,
gcc 4.x, but are present in Linux.  We implement them here for
OS X, at least until the underlying mechanism is no longer
supported by Apple.

The mechanism is to use the POSIX functions fegetenv() and
fesetenv(), which *are* present in OS X, to manipulate the ppc
and intel floating-point control registers, after changing bits
in fields corresponding to those registers in the fenv_t data
type.

Assembly language code to directly access the floating-point
status and control registers for ppc and intel is also included.

This example grew out of an update to legacy code for Apple
ppc's.  The original legacy code is in Listing 7-1 in "PowerPC
Numerics", 2004:

http://lists.apple.com/archives/unix-porting/2003/May/msg00026.html

Another version of the ppc legacy code is here:

http://developer.apple.com/documentation/Performance/Conceptual/Mac_OSX_Numerics/Mac_OSX_Numerics.pdf

Terry Lambert pointed out that our naive update of the legacy
example to Mac OS X Leopard made egregious unsupported use of
system context structures in the handler.  See his reply to

http://lists.apple.com/archives/Darwin-dev/2009/Sep/msg00091.html

The example in this file is more plain vanilla, and aims at
alternate handling that does not return to the application, but
rather aborts with a diagnostic message.

To compile it under Mac OS X, execute:

  cc -o fe-handling fe-handling-example.c

To compile it under Linux, execute:

  cc -DLINUX -lm -o fe-handling fe-handling-example.c
*/

/* Need this to evaluate the FCLAW_HAVE_ definitions */
#include <fclaw_config.h>

/* We only do something if we need the feenableexcept replacement.
   Otherwise this file noops. */
#ifndef FCLAW_HAVE_FEENABLEEXCEPT

/* Added this ifdef to compile the file conditionally */
/* Question: why does this file not include fp_exception_glibc_extension.h? */
/* Answer: now it does -- is this how we like it? */
#include <fp_exception_glibc_extension.h>

#ifdef FCLAW_HAVE_SIGNAL_H
#include <signal.h>
#endif

#ifdef LINUX
/* BEGIN quote
http://graphviz.sourcearchive.com/documentation/2.16/gvrender__pango_8c-source.html
*/
/* _GNU_SOURCE is needed (supposedly) for the feenableexcept
 * prototype to be defined in fenv.h on GNU systems.
 * Presumably it will do no harm on other systems.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/* We are not supposed to need __USE_GNU, but I can't see
 * how to get the prototype for fedisableexcept from
 * /usr/include/fenv.h without it.
 */
#ifndef __USE_GNU
#define __USE_GNU
#endif
/* END quote */
#endif // LINUX

#ifdef FCLAW_HAVE_FENV_H
#include <fenv.h>
#endif

#define DEFINED_PPC      (defined(__ppc__) || defined(__ppc64__))
#define DEFINED_INTEL    (defined(__i386__) || defined(__x86_64__))

#ifndef LINUX
#if DEFINED_PPC

#define FE_EXCEPT_SHIFT 22  // shift flags right to get masks
#define FM_ALL_EXCEPT    FE_ALL_EXCEPT >> FE_EXCEPT_SHIFT

/* GNU C Library:
http://www.gnu.org/software/libc/manual/html_node/Control-Functions.html

     - Function: int fegetexcept (int excepts)

       The function returns a bitmask of all currently enabled
       exceptions.  It returns -1 in case of failure.

   The excepts argument appears in other functions in fenv.h,
   and corresponds to the FE_xxx exception flag constants.  It
   is unclear whether the bitmask is for the flags or the masks.
   We return that for the flags, which corresponds to the
   excepts argument in feenableexcept(excepts) and
   fedisableexcept(excepts).  In GNU/Linux the argument is void,
   and that's what we implement.  Linux "man fegetenv" appears
   to suggest that it's the mask corresponding to bits in
   excepts that is returned.
*/
int
fegetexcept (void)
{
  static fenv_t fenv;

  return ( fegetenv (&fenv) ? -1 :
    (
      ( fenv & (FM_ALL_EXCEPT) ) << FE_EXCEPT_SHIFT )
    );
}

int
feenableexcept (int excepts)
{
  static fenv_t fenv;
  int new_excepts = (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;

  fenv = (fenv & ~new_excepts) | new_excepts;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

int
fedisableexcept (int excepts)
{
  static fenv_t fenv;
  int still_on = ~( (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT ),
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;

  fenv &= still_on;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#elif DEFINED_INTEL

int
fegetexcept (void)
{
  static fenv_t fenv;

  return fegetenv (&fenv) ? -1 : (fenv.__control & FE_ALL_EXCEPT);
}

int
feenableexcept (int excepts)
{
  static fenv_t fenv;
  int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

int
fedisableexcept (int excepts)
{
  static fenv_t fenv;
  int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // mask
  fenv.__control |= new_excepts;
  fenv.__mxcsr   |= new_excepts << 7;

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#endif  // PPC or INTEL enabling
#endif  // not LINUX

#if DEFINED_PPC

#define getfpscr(x)    asm volatile ("mffs %0" : "=f" (x));
#define setfpscr(x)    asm volatile ("mtfsf 255,%0" : : "f" (x));

typedef union {
    struct {
        unsigned long hi;
        unsigned long lo;
    } i;
    double d;
} hexdouble;

#endif  // DEFINED_PPC

#if DEFINED_INTEL

// x87 fpu
/* DAC (9/26/2016) : Changed 'asm' to '__asm__ volatile' so that macro expansion
   passed -std=c99 -pedentic compiler flags */

#define getx87cr(x)    __asm__ volatile ("fnstcw %0" : "=m" (x));
#define setx87cr(x)    __asm__ volatile ("fldcw %0" : "=m" (x));
#define getx87sr(x)    __asm__ volatile ("fnstsw %0" : "=m" (x));

// SIMD, gcc with Intel Core 2 Duo uses SSE2(4)
#define getmxcsr(x)    __asm__ volatile ("stmxcsr %0" : "=m" (x));
#define setmxcsr(x)    __asm__ volatile ("ldmxcsr %0" : "=m" (x));

#endif  // DEFINED_INTEL

#ifdef FCLAW_HAVE_SIGNAL_H
#include <signal.h>
#endif

#include <stdio.h>   // printf()
#include <stdlib.h>  // abort(), exit()

#if 0

static const char *fe_code_name[] = {
  "FPE_NOOP",
  "FPE_FLTDIV", "FPE_FLTINV", "FPE_FLTOVF", "FPE_FLTUND",
  "FPE_FLTRES", "FPE_FLTSUB", "FPE_INTDIV", "FPE_INTOVF"
  "FPE_UNKNOWN"
};

/* SAMPLE ALTERNATE FP EXCEPTION HANDLER

   The sample handler just reports information about the
   exception that invoked it, and aborts.  It makes no attempt
   to restore state and return to the application.

   More sophisticated handling would have to confront at least
   these issues:

     * interface to the system context for restoring state
     * imprecision of interrupts from hardware for the intel x87
       fpu (but not the SIMD unit, nor the ppc)
     * imprecision of interrupts from system software
*/
void
fhdl ( int sig, siginfo_t *sip, ucontext_t *scp )
{
  int fe_code = sip->si_code;
  int excepts = fetestexcept (FE_ALL_EXCEPT);

  switch (fe_code)
  {
#ifdef FPE_NOOP  // occurs in OS X
    case FPE_NOOP:   fe_code = 0; break;
#endif
    case FPE_FLTDIV: fe_code = 1; break; // divideByZero
    case FPE_FLTINV: fe_code = 2; break; // invalid
    case FPE_FLTOVF: fe_code = 3; break; // overflow
    case FPE_FLTUND: fe_code = 4; break; // underflow
    case FPE_FLTRES: fe_code = 5; break; // inexact
    case FPE_FLTSUB: fe_code = 6; break; // invalid
    case FPE_INTDIV: fe_code = 7; break; // overflow
    case FPE_INTOVF: fe_code = 8; break; // underflow
            default: fe_code = 9;
   }

  if ( sig == SIGFPE )
  {
#if DEFINED_INTEL
    unsigned short x87cr,x87sr;
    unsigned int mxcsr;

    getx87cr (x87cr);
    getx87sr (x87sr);
    getmxcsr (mxcsr);
    printf ("X87CR:   0x%04X\n", x87cr);
    printf ("X87SR:   0x%04X\n", x87sr);
    printf ("MXCSR:   0x%08X\n", mxcsr);
#endif

#if DEFINED_PPC
   hexdouble t;

   getfpscr (t.d);
   printf ("FPSCR:   0x%08X\n", t.i.lo);
#endif

    printf ("signal:  SIGFPE with code %s\n", fe_code_name[fe_code]);
    printf ("invalid flag:    0x%04X\n", excepts & FE_INVALID);
    printf ("divByZero flag:  0x%04X\n", excepts & FE_DIVBYZERO);
  }
  else printf ("Signal is not SIGFPE, it's %i.\n", sig);

  abort();
}
#endif

#if 0
int main (int argc, char **argv)
{
    double s;
    struct sigaction act;

    act.sa_sigaction = (void(*))fhdl;
    sigemptyset (&act.sa_mask);
    act.sa_flags = SA_SIGINFO;


//  printf ("Old divByZero exception: 0x%08X\n", feenableexcept (FE_DIVBYZERO));
    printf ("Old invalid exception:   0x%08X\n", feenableexcept (FE_INVALID));
    printf ("New fp exception:        0x%08X\n", fegetexcept ());

    // set handler
    if (sigaction(SIGFPE, &act, (struct sigaction *)0) != 0)
    {
        perror("Yikes");
        exit(-1);
    }

//  s = 1.0 / 0.0;  // FE_DIVBYZERO
    s = 0.0 / 0.0;  // FE_INVALID
    return 0;
}
#endif

#endif /* !FCLAW_HAVE_FEENABLEEXCEPT */
