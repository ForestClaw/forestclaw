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

#ifndef METRIC_USER_H
#define METRIC_USER_H

#include <fclaw2d_forestclaw.h>
#include <fc2d_clawpack46.h>

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
    double beta;
    int claw_version;

    int is_registered;

} user_options_t;

#define METRIC_SETPROB FCLAW_F77_FUNC(metric_setprob, METRIC_SETPROB)
void METRIC_SETPROB(const double* beta);

const user_options_t* metric_user_get_options(fclaw2d_domain_t* domain);


void metric_problem_setup(fclaw2d_domain_t* domain);



void metric_check_(const int& mbc,const int& mx, const int& my,
                   const double& xlower, const double& ylower,
                   const double& dx,const double& dy,
                   const int& blockno,
                   const double xp[], double yp[],const double zp[],
                   const double xd[], double yd[],const double zd[],
                   const double area[], const int& mpirank);

void compute_error(int* blockno, int* mx, int* my,
                   int* mbc, int* meqn,
                   double* dx, double* dy,
                   double* xlower, double* ylower,
                   double *t, double q[],double error[]);

void initialize(const int& mx, const int& my,
                const int& meqn, const int& mbc,
                const double& xlower, const double& ylower,
                const double& dx, const double& dy,
                double q[],
                const double curvature[], const double area[]);

double total_area_(const int& mx, const int& my, const int& mbc,
                   const double area[]);

void min_grid_cell_area_(const int& mx, const int& my, const int& mbc,
                         const double& dx, const double& dy,
                         const double area[],double *minvalue);

void max_grid_cell_area_(const int& mx, const int& my, const int& mbc,
                         const double& dx, const double& dy,
                         const double area[], double *maxvalue);

void metric_patch_initialize(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx);

void metric_link_patch(fclaw2d_domain_t *domain);

void metric_diagnostics(fclaw2d_domain_t *domain, const double t);


fclaw_map_context_t* fclaw2d_map_new_nomap();

fclaw_map_context_t* fclaw2d_map_new_cart(const double scale[],
                                            const double shift[],
                                            const double rotate[]);


fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double rotate[],
                                                 const double alpha);

fclaw_map_context_t* fclaw2d_map_new_pillowdisk(const double scale[],
                                                  const double shift[],
                                                  const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_squareddisk(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

fclaw_map_context_t* fclaw2d_map_new_pillowsphere(const double scale[],
                                                    const double shift[],
                                                    const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_cubedsphere(const double scale[],
                                                   const double shift[],
                                                   const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

fclaw_map_context_t* fclaw2d_map_new_torus(const double scale[],
                                             const double shift[],
                                             const double rotate[],
                                             const double alpha);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
