/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#ifndef RADIALDAM_USER_H
#define RADIALDAM_USER_H

#include <fclaw_include_all.h>

#include <fclaw_include_all.h>

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>

#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack5_options.h>

#include <fc2d_clawpack46.h>
#include <fc2d_clawpack5.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int example;

    double alpha;

    double g;
    double x0;
    double y0;
    double r0;
    double hin;
    double hout;

    int claw_version;

    int is_registered;
} user_options_t;


void radialdam_link_solvers(fclaw_global_t *glob);

user_options_t* radialdam_options_register (fclaw_app_t * app,
                                          const char *configfile);

void radialdam_options_store (fclaw_global_t* glob, user_options_t* user);

user_options_t* radialdam_get_options(fclaw_global_t* glob);


void radialdam_global_post_process(fclaw_options_t *fclaw_opt,
                                   fclaw_clawpatch_options_t *clawpatch_opt,
                                   user_options_t *user_opt);


/* ------------------------------- Mapping functions ---------------------------------- */

fclaw_map_context_t* fclaw2d_map_new_nomap();

fclaw_map_context_t* fclaw2d_map_new_pillowdisk(const double scale[],
                                                  const double shift[],
                                                  const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha);
/* ----------------------------- Conservative update ---------------------------------- */

#define RPN2_CONS_UPDATE FCLAW_F77_FUNC(rpn2_cons_update,RPN2_CONS_UPDATE)

void RPN2_CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
                      double q[], double aux_center[], double aux_edge[], double flux[]);


#define RPN2_CONS_UPDATE_MANIFOLD FCLAW_F77_FUNC(rpn2_cons_update_manifold, \
                                                 RPN2_CONS_UPDATE_MANIFOLD)

void RPN2_CONS_UPDATE_MANIFOLD(const int* meqn, const int* maux, const int* idir,
                               const int* iface,
                               double q[], double aux_center[], double aux_edge[],
                               double flux[]);

#define RPN2_CONS_UPDATE_ZERO FCLAW_F77_FUNC(rpn2_cons_update_zero, \
                                                 RPN2_CONS_UPDATE_ZERO)

void RPN2_CONS_UPDATE_ZERO(const int* meqn, const int* maux, const int* idir,
                           const int* iface,
                           double q[], double aux_center[], double aux_edge[],
                           double flux[]);

/* ----------------------------- Setaux for manifolds --------------------------------- */

#define USER46_SETAUX_MANIFOLD FCLAW_F77_FUNC(user46_setaux_manifold, \
                                             USER46_SETAUX_MANIFOLD)

void USER46_SETAUX_MANIFOLD(const int* mbc,
                            const int* mx, const int* my,
                            const double* xlower, const double* ylower,
                            const double* dx, const double* dy,
                            const int* maux, double aux[],
                            double xnormals[], double xtangents[],
                            double ynormals[], double ytangents[],
                            double surfnormals[],
                            double area[]);

#define USER5_SETAUX_MANIFOLD FCLAW_F77_FUNC(user5_setaux_manifold, \
                                             USER5_SETAUX_MANIFOLD)

void USER5_SETAUX_MANIFOLD(const int* mbc,
                           const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           const int* maux, double aux[],
                           double xnormals[], double xtangents[],
                           double ynormals[], double ytangents[],
                           double surfnormals[],
                           double area[]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
