/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_DIAGNOSTICS_H
#define FCLAW2D_DIAGNOSTICS_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
/* Now stored in vtable */
typedef void (*fclaw2d_run_diagnostics_t)(fclaw2d_domain_t *domain,
                                          const double t);
#endif

/* See forestclaw2d.h for the maximum version of this function */
double fclaw2d_domain_global_minimum (fclaw2d_domain_t* domain, double d);

#if 0
void link_run_diagnostics(fclaw2d_domain_t* domain,
                          fclaw2d_run_diagnostics_t f_run_diagnostics);

void run_diagnostics_default(fclaw2d_domain_t* domain, const double t);
#endif

void fclaw2d_run_diagnostics(fclaw2d_domain_t *domain);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
