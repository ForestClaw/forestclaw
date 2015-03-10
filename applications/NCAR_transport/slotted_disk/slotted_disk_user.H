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

#ifndef SLOTTED_DISK_USER_H
#define SLOTTED_DISK_USER_H

#include <fclaw2d_clawpatch.H>
#include <fc2d_clawpack46.H>

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
    int vflag;       /* specify edge normal velocities */
    int init_choice; /* for slotted disk */
    int is_registered;

} user_options_t;

void slotted_disk_problem_setup(fclaw2d_domain_t *domain);

#define SETPROB_TRANSPORT FCLAW_F77_FUNC(setprob_transport,SETPROB_TRANSPORT)
void SETPROB_TRANSPORT(const int* vflag, const int* ichoice);

/* Initialization */
void slotted_disk_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx);

#define QINIT_TRANSPORT FCLAW_F77_FUNC(qinit_transport,QINIT_TRANSPORT)
void QINIT_TRANSPORT(const int* mx, const int* my, const int* meqn,
                     const int* mbc,
                     const double* xlower, const double* ylower,
                     const double* dx, const double* dy,
                     double q[], const int* maux, double aux[],
                     const int* blockno,
                     double xp[], double yp[], double zp[]);

void slotted_disk_patch_setup(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx);

#define SETAUX_TRANSPORT FCLAW_F77_FUNC(setaux_transport,SETAUX_TRANSPORT)
void SETAUX_TRANSPORT(const int* mx, const int* my, const int* mbc,
                      const double* xlower, const double* ylower,
                      const double* dx, const double* dy,
                      const int* maux, double aux[],
                      const int* blockno,
                      double xd[], double yd[], double zd[],double area[]);

void slotted_disk_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt);

#define B4STEP2_TRANSPORT FCLAW_F77_FUNC(b4step2_transport,B4STEP2_TRANSPORT)
void B4STEP2_TRANSPORT(const int* mx, const int* my, const int* mbc,
                       const double* dx, const double* dy,
                       const double* t, const int* maux, double aux[],
                       const int* blockno,
                       double xd[], double yd[], double zd[]);

double slotted_disk_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt);

void slotted_disk_link_solvers(fclaw2d_domain_t *domain);


fclaw2d_map_context_t * fclaw2d_map_new_cubedsphere (const double scale[],
                                                     const double shift[],
                                                     const double rotate[]);

fclaw2d_map_context_t * fclaw2d_map_new_pillowsphere (const double scale[],
                                                      const double shfit[],
                                                      const double rotate[]);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
