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

#include "metric_user.h"

#include <fclaw2d_diagnostics.h>
#include <fclaw_clawpatch.hpp>

static
double metric_surface_area(fclaw2d_map_context_t* cont)
{
    double exact_area;

    double alpha = 0.4;
    if (FCLAW2D_MAP_IS_DISK(&cont))
    {
        exact_area = M_PI;
    }
    else if (FCLAW2D_MAP_IS_SPHERE(&cont))
    {
        exact_area = 4.0*M_PI;
    }
    else if (FCLAW2D_MAP_IS_TORUS(&cont))
    {
        exact_area = 4.0*alpha*(M_PI)*(M_PI);
    }
    else
    {
        exact_area = 4.0;
    }
    return exact_area;
}


static
void cb_total_area(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   void *user)
{
    double *sum = (double*) user;

    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
#if 0
    ClawPatch *cp = fclaw2d_clawpatch_cp(domain,this_patch);
#endif
    ClawPatch* cp = fclaw2d_clawpatch_get_cp(this_patch);
    double *area = cp->area();

    *sum += total_area_(mx,my,mbc,area);
}

static
void cb_min_cell_area(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    double *minvalue = (double*) user;

    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    double *area = cp->area();
    double dx = cp->dx();
    double dy = cp->dy();

    min_grid_cell_area_(mx,my,mbc,dx,dy,area,minvalue);
}

static
void cb_max_cell_area(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    double *maxvalue = (double*) user;
    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    double *area = cp->area();
    double dx = cp->dx();
    double dy = cp->dy();

    max_grid_cell_area_(mx,my,mbc,dx,dy,area,maxvalue);
}

void metric_diagnostics(fclaw2d_domain_t *domain, const double t)
{
    /* Compute a global sum */
    double sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_total_area,(void *) &sum);
    sum = fclaw2d_domain_global_sum (domain, sum);

    fclaw2d_map_context_t* cont = fclaw2d_domain_get_map_context(domain);
    double exact_area = metric_surface_area(cont);

    fclaw_global_productionf("%30s %24.16f\n","Total area",sum);
    fclaw_global_productionf("%30s %24.16f\n","Exact area",exact_area);
    fclaw_global_productionf("%30s %24.4e\n","Error ",fabs(exact_area-sum));
    fclaw_global_productionf("\n");

    const fclaw_options_t *gparms = get_domain_parms(domain);
    if (gparms->minlevel == gparms->maxlevel)
    {
        /* Only compare ratio of smallest grid cell to largest if the grid is
           uniformly refined */

        double minvalue = 100;
        fclaw2d_domain_iterate_patches(domain,cb_min_cell_area,(void *) &minvalue);
        minvalue = fclaw2d_domain_global_minimum (domain, minvalue);

        double maxvalue = 0;
        fclaw2d_domain_iterate_patches(domain,cb_max_cell_area,(void *) &maxvalue);
        maxvalue = fclaw2d_domain_global_maximum (domain, maxvalue);

        fclaw_global_productionf("%30s %24.8e\n","Minimum value (kappa)",minvalue);
        fclaw_global_productionf("%30s %24.8e\n","Maximum value (kappa)",maxvalue);
        fclaw_global_productionf("%30s %24.8f\n","Ratio of largest to smallest",
                                 maxvalue/minvalue);
        fclaw_global_productionf("\n\n");
    }
}
