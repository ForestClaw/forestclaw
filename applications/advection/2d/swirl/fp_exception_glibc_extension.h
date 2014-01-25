/* 
   glibc floating point extension replacement for OS X. 
   
   Title:  Floating-point exception handling example
   Author:  David N. Williams
   File:  fe-handlng-example.c
   License:  Public Domain
   Version:  0.5.0

*/

#ifndef __FENV_H_LINUX_REPLACEMENT__
#define __FENV_H_LINUX_REPLACEMENT__

//#include <fenv.h>
#include <signal.h>
#ifdef __cplusplus
extern "C" {
#endif

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
  int fegetexcept (void) ;
  int feenableexcept (unsigned int excepts) ;
  int fedisableexcept (unsigned int excepts) ;
  
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
  void fhdl ( int sig, siginfo_t *sip, ucontext_t *scp ) ;

#ifdef __cplusplus
}
#endif


#endif // __FENV_H_LINUX_REPLACEMENT__
