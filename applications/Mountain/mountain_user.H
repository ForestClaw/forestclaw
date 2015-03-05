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

#ifndef MOUNTAIN_USER_H
#define MOUNTAIN_USER_H

#include <fc2d_clawpack46.H>
#include <amr_forestclaw.H>

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
    int is_registered;

} user_options_t;

void mountain_patch_setup(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_blockno,
                          int this_patchno);


#define SET_CAPACITY FCLAW_F77_FUNC(set_capacity,SET_CAPACITY)
void SET_CAPACITY(int* mx, int *my, int *mbc, double *dx,double *dy,
                  double area[],int* maux, double *aux);


#if 0
    void mountain_tag4refinement_(const int& mx,const int& my,const int& mbc,
                               const int& meqn, const int& blockno,
                               const double& xlower, const double& ylower,
                               const double& dx, const double& dy,
                               double q[], const int& init_flag, int& tag_patch);

    void mountain_tag4coarsening_(const int& mx,const int& my,const int& mbc,
                               const int& meqn,const double& xlower,
                               const double& ylower, const double& dx, const double& dy,
                               double qcoarsened[], int& tag_patch);

    void mountain_write_tfile_(const int& iframe,const double& time,
                            const int& mfields,const int& ngrids,
                            const int& maux);

    void mountain_write_qfile_(const int& meqn,
                            const int& mbc, const int& mx,const int& my,
                            const double& xlower,const double& ylower,
                            const double& dx, const double& dy,
                            double q[], const int& iframe,const int& patch_num,
                            const int& level,const int& blockno,
                            const int& mpirank);
#endif

#if 0
    void qinit_manifold_(const int& meqn, const int& mbc,
                         const int& mx, const int& my,
                         const double& xlower, const double& ylower,
                         const double& dx, const double& dy,
                         const int& blockno,double q[],
                         const int& maux, double aux[]);

    void setaux_manifold_(const int& mbc,const int& mx, const int& my,
                          const double& xlower, const double& ylower,
                          const double& dx, const double& dy,
                          const int& maux, double aux[],
                          const int& blockno,
                          double xp[], double yp[], double zp[],
                          double xd[], double yd[], double zd[],
                          double area[]);

    void b4step2_manifold_(const int& mbc,const int& mx,const int& my,
                           const int& meqn, double q[],
                           const double& dx,const double& dy,
                           const int& blockno,
                           double xd[], double yd[], double zd[],
                           const double& time,
                           const double& dt, const int& maux, double aux[]);
#endif


void mountain_link_solvers(fclaw2d_domain_t *domain);


double mountain_patch_update(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt);

/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_mountain(fclaw2d_map_context_t* brick,
                                                const double scale[],
                                                const double shift[]);

fclaw2d_map_context_t* fclaw2d_map_new_identity(fclaw2d_map_context_t* brick,
                                                const double scale[],
                                                const double shift[]);

#define MAPC2M_MOUNTAIN FCLAW_F77_FUNC (mapc2m_mountain,MAPC2M_MOUNTAIN)
void MAPC2M_MOUNTAIN (int* blockno, double *xc, double *yc,
                     double *xp, double *yp, double *zp, double *scale);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
