/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <fclaw_fpe.h>
#include <fclaw2d_global.h>

#include <signal.h>    /* for signal()  */
#include <fenv.h>      /* for fegetenv(), fesetenv() */

void fclaw_fpe_signal_handler_default(int sig) 
{
    fclaw_global_essentialf("Floating point exception encountered\n");
    switch (sig)
    {
        case FE_INVALID:
            fclaw_global_essentialf("Invalid value detected\n");
            exit(0);
            break;
        case FE_DIVBYZERO:
            fclaw_global_essentialf("Division by zero detected\n");
            exit(0);
            break;
        case FE_OVERFLOW:
            fclaw_global_essentialf("Overflow detected\n");
            exit(0);
            break;
        case FE_UNDERFLOW:
        case FE_INEXACT:
            /* Ignore these - too many exceptions detected */
            break;
        default : 
            break;
    }
    return;
}


void fclaw_fpe_handling_default(fclaw2d_global_t* glob)
{
    fclaw_global_essentialf("Floating point exception handling is disabled.\n");
}
