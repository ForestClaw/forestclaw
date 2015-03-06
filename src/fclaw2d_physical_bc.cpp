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

#include "forestclaw2d.h"
#include "fclaw2d_physical_bc.h"
#include "fclaw2d_vtable.h"


void fclaw2d_physbc_default(fclaw2d_domain *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            double t,
                            double dt,
                            fclaw_bool intersects_phys_bdry[],
                            fclaw_bool time_interp)
{
    /* This can be used when no BCs are to be called */
}

static
void cb_set_phys_bc(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    struct time_info_t {double level_time; fclaw_bool time_interp; } *t_info;
    t_info = (time_info_t*) user;

    fclaw_bool intersects_bc[NumFaces];
    double dt = 1e20;
    fclaw2d_get_physical_bc(domain,this_block_idx,this_patch_idx,intersects_bc);

    vt.patch_physical_bc(domain,
                         this_patch,
                         this_block_idx,
                         this_patch_idx,
                         t_info->level_time,dt,
                         intersects_bc,
                         t_info->time_interp);
}

/* This is needed by other routines, so we don't set it to static. */
void fclaw2d_get_physical_bc(fclaw2d_domain_t *domain,
                             int this_block_idx,
                             int this_patch_idx,
                             fclaw_bool *intersects_bdry)
{
    // const int numfaces = get_faces_per_patch(domain);
    int bdry[NumFaces];
    fclaw2d_patch_boundary_type(domain,this_block_idx,this_patch_idx,bdry);
    for(int i = 0; i < NumFaces; i++)
    {
        // Physical boundary conditions
        intersects_bdry[i] = bdry[i] == 1;
    }
}


/* -----------------------------------------------------------------------------
   Set physical boundary conditions on a patch
   ----------------------------------------------------------------------------- */

void fclaw2d_set_physical_bc(fclaw2d_domain_t *domain, int a_level,
                             double a_level_time, fclaw_bool time_interp)
{
    struct time_info_t {double level_time; fclaw_bool time_interp; } t_info;
    t_info.level_time = a_level_time;
    t_info.time_interp = time_interp;
    fclaw2d_domain_iterate_level(domain, a_level,cb_set_phys_bc,(void *) &t_info);
}
