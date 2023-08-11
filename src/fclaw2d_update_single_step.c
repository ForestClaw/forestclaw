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

#include <fclaw2d_update_single_step.h>

#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>

static
void cb_single_step_count(fclaw_domain_t *domain,
                          fclaw_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{    
    fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;
    int *count = (int*) g->user;
    (*count)++;
}

static
void cb_single_step(fclaw_domain_t *domain,
                    fclaw_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user)
{
    fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;

    double maxcfl;
    
    fclaw2d_single_step_data_t *ss_data = (fclaw2d_single_step_data_t *) g->user;
    double dt = ss_data->dt;
    double t = ss_data->t;
    
    maxcfl = fclaw_patch_single_step_update(g->glob,this_patch,
                                              this_block_idx,
                                              this_patch_idx,t,dt,
                                              &ss_data->buffer_data);

    ss_data->buffer_data.iter++;  /* Used for patch buffer */
    g->glob->count_single_step++;

    /* Note : In openMP, this reduction is probably not correct */
    ss_data->maxcfl = fmax(maxcfl,ss_data->maxcfl);

}


double fclaw2d_update_single_step(fclaw_global_t *glob,
                                  int level,
                                  double t, double dt)
{

    /* Iterate over every patch at this level */
    fclaw2d_single_step_data_t ss_data;
    ss_data.t = t;
    ss_data.dt = dt;
    ss_data.maxcfl = 0;
    ss_data.buffer_data.total_count = 0;
    ss_data.buffer_data.iter = 0;

    /* If there are not grids at this level, we return CFL = 0 */
#if defined(_OPENMP)        
    fclaw_global_iterate_level_mthread(glob, level, 
                                         cb_single_step,(void *) &ss_data);
#else
    /* Count number of grids to be updated in this call;  not sure how this
       works in OpenMP.  */
    int count = 0;
    fclaw_global_iterate_level(glob, level, cb_single_step_count,&count);
    ss_data.buffer_data.total_count = count;

    fclaw_global_iterate_level(glob, level, 
                                 cb_single_step,(void *) &ss_data);
#endif   
 

    return ss_data.maxcfl;
}
