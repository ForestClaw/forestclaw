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

#ifndef AMR_MOL_H
#define AMR_MOL_H

#include "forestclaw2d.h"
#include "amr_options.h"

// Needed to avoid circular typedef refs.  This struct is defined in
// 'fclaw_typedefs.h', which includes this file.
struct fclaw2d_level_time_data;

typedef struct fclaw_level_mol_data fclaw_level_mol_data_t;
typedef struct fclaw2d_f_rhs_data_fort fclaw2d_f_rhs_data_fort_t;

// These are used to define domain data in 'fclaw_typedefs.h'
typedef void (*fclaw_mol_rhs_patch_t)(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx,
                                      int this_patch_idx,
                                      double t,
                                      double *rhs);

// Interface to MOL solver
typedef void (*fclaw_mol_solver_t)(int neqn,double q[],
                                   double t, double dt);


struct fclaw_level_mol_data
{
    int count;
    int level;
    int patch_size;
    fclaw2d_patch_t **patches;   // Store patches at a single level.
    int *patch_indices;
    int *block_indices;
    double *patch_data;
    double *patch_data_restore;
};


struct fclaw2d_f_rhs_data_fort
{
    fclaw2d_domain_t* domain;
    struct fclaw2d_level_time_data *time_data;
    struct fclaw_level_mol_data *mol_data;
};

extern "C"
{

    /* Fortran routine that does actual time stepping */
    void feuler_(void* f, const int& neqn, double q[],
                 const double& t, const double& dt,
                 double rhs_work[]);

    void rkc_(int& neqn, void* f, double q_amr[], double& ti,
              double& tip1, double& rtol, double& atol,
              int info[], double rhs_work[], int& idid);
}

// This is called from Fortran, so needs to be in this 'extern C' scope.
void fclaw_mol_rhs(const double& t_inner, double *q, double *rhs_work);

void fclaw_mol_step(fclaw2d_domain_t *domain,
                    int level, struct fclaw2d_level_time_data *time_data);

void parabolic_step_feuler(int neqn, double q[], double t, double dt);
void parabolic_step_rkc(int neqn, double q[], double t, double dt);


#endif
