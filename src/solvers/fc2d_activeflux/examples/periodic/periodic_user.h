/*
Copyright (c) 2012-2020 Carsten Burstedde, Donna Calhoun, Erik Chudzik, Christiane Helzel
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

#ifndef PERIODIC_USER_H
#define PERIODIC_USER_H

#include <fclaw2d_include_all.h>

#include <fc2d_activeflux.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* ------------------------------------- Options ---------------------------------------*/
typedef struct user_options
{
    int example; 
    
    double ubar;
    double vbar;

    int refinement_strategy;

    int claw_version;
    int is_registered;

} user_options_t;


user_options_t* periodic_options_register (fclaw_app_t * app, const char *configfile);

void periodic_options_store (fclaw2d_global_t* glob, user_options_t* user);

const user_options_t* periodic_get_options(fclaw2d_global_t* glob);

/* ------------------- Public functions defined in periodic_user -----------------------*/

void periodic_link_solvers(fclaw2d_global_t *glob);

/* ------------------------------------ Fortran ----------------------------------------*/
#define PERIODIC_SETPROB FCLAW_F77_FUNC(periodic_setprob, PERIODIC_SETPROB)
void PERIODIC_SETPROB();

#define PERIODIC_COMPUTE_ERROR FCLAW_F77_FUNC(periodic_compute_error,PERIODIC_COMPUTE_ERROR)

void PERIODIC_COMPUTE_ERROR(int* blockno, int *mx, int *my, int* mbc, int* meqn,
                           double *dx, double *dy, double *xlower,
                           double *ylower, double *t, double q[],
                           double error[], double soln[]);

#define RPN2CONS_UPDATE FCLAW_F77_FUNC(rpn2cons_update,RPN2CONS_UPDATE)

void RPN2CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
                      double q[], double aux_center[], double aux_edge[], double flux[]);



/* ------------------------------- Other Fortran headers -------------------------------*/

/* Other headers (e.g. activeflux_rpn2, clawpack2_rpt2) are  provided by default in 
   src/solvers/fc2d_clawpack4.6/activeflux_fort_user.h */


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
