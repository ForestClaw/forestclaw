/*
  Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef WAVETANK_USER_H
#define WAVETANK_USER_H

#include <fclaw2d_include_all.h>
#include <fc2d_thunderegg_fort.h>  /* For virtualized functions */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
    /* Fix syntax highlighting */
#endif    

typedef struct user_options
{
    int example;
    double g;

    double dry_tolerance;
    double sea_level;

    int claw_version;

    int is_registered;
} user_options_t;


#define WAVETANK_SETPROB FCLAW_F77_FUNC(wavetank_setprob, WAVETANK_SETPROB)
void WAVETANK_SETPROB();


//void wavetank_problem_setup(fclaw2d_global_t *glob);
void wavetank_link_solvers(fclaw2d_global_t *glob);

user_options_t* wavetank_options_register(fclaw_app_t * app,
                                          const char *configfile);

void wavetank_options_store(fclaw2d_global_t* glob, user_options_t* user);

user_options_t* wavetank_get_options(fclaw2d_global_t* glob);

void wavetank_run(fclaw2d_global_t *glob);

fclaw2d_map_context_t* fclaw2d_map_new_nomap();


void fc2d_geoclaw_output_ascii(fclaw2d_global_t* glob,int iframe);


/* ***************************** FORTRAN - Riemann solvers **************************** */

#define RPN2_TSUNAMI    FCLAW_F77_FUNC(rpn2_tsunami, RPN2_TSUNAMI)
void RPN2_TSUNAMI(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc,const int* mx, double ql[], double qr[],
                  double auxl[], double auxr[], double wave[],
                  double s[], double amdq[], double apdq[]);

#define RPT2_TSUNAMI    FCLAW_F77_FUNC(rpt2_tsunami, RPT2_TSUNAMI)
void RPT2_TSUNAMI(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc, const int* mx, double ql[], double qr[],
                  double aux1[], double aux2[], double aux3[], const int* imp,
                  double dsdq[], double bmasdq[], double bpasdq[]);

#define RPN2_GEOCLAW    FCLAW_F77_FUNC(rpn2_geoclaw,   RPN2_GEOCLAW)
void RPN2_GEOCLAW(const int* ixy,const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc,const int* mx, double ql[], double qr[],
                  double auxl[], double auxr[], double wave[],
                  double s[], double amdq[], double apdq[]);

#define RPT2_GEOCLAW    FCLAW_F77_FUNC(rpt2_geoclaw,   RPT2_GEOCLAW)
void RPT2_GEOCLAW(const int* ixy, const int* maxm, const int* meqn, const int* mwaves,
                  const int* mbc, const int* mx, double ql[], double qr[],
                  double aux1[], double aux2[], double aux3[], const int* imp,
                  double dsdq[], double bmasdq[], double bpasdq[]);


/* --------------------------------- Output functions ----------------------------------- */

#define FC2D_GEOCLAW_FORT_WRITE_HEADER FCLAW_F77_FUNC(fc2d_geoclaw_fort_write_header,\
                                                      FC2D_GEOCLAW_FORT_WRITE_HEADER)
void FC2D_GEOCLAW_FORT_WRITE_HEADER(int* iframe, double* time, int* meqn, 
                                    int* maux, int* ngrids);

#define FC2D_GEOCLAW_FORT_WRITE_FILE FCLAW_F77_FUNC(fc2d_geoclaw_fort_write_file, \
                                                    FC2D_GEOCLAW_FORT_WRITE_FILE)
void FC2D_GEOCLAW_FORT_WRITE_FILE(int* mx,int* my,int* meqn,int* maux, 
                                  int* mbathy,int* mbc, 
                                  double* xlower, double* ylower, 
                                  double* dx,double* dy,
                                  double q[], double aux[], int* iframe, 
                                  int* patch_num,int* level,
                                  int* blockno,int* mpirank);

#define WAVETANK_FORT_TAG4REFINEMENT FCLAW_F77_FUNC(wavetank_fort_tag4refinement, \
                                                 WAVETANK_FORT_TAG4REFINEMENT)

void WAVETANK_FORT_TAG4REFINEMENT(const int* mx,const int* my,
                                  const int* mbc,const int* meqn,
                                  const double* xlower, const double* ylower,
                                  const double* dx, const double* dy,
                                  double* time,
                                  const int* blockno,
                                  double q[],
                                  const double* tag_threshold,
                                  int* level,
                                  const int* init_flag,
                                  int* tag_patch);

#define WAVETANK_FORT_TAG4COARSENING FCLAW_F77_FUNC(wavetank_fort_tag4coarsening, \
                                                WAVETANK_FORT_TAG4COARSENING)

void WAVETANK_FORT_TAG4COARSENING(const int* mx, const int* my,
                                  const int* mbc, const int* meqn,
                                  double xlower[], double ylower[],
                                  const double* dx, const double* dy,
                                  double* time,
                                  const int* blockno,
                                  double q0[],double q1[],
                                  double q2[],double q3[],
                                  const double* tag_threshold,
                                  const int* initflag,
                                  int* tag_patch);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
