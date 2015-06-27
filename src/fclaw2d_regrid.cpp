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

#include <fclaw2d_global.h>

#include "fclaw2d_regrid.h"
#include "fclaw2d_ghost_fill.h"

#include <forestclaw2d.h>
#include <amr_utils.H>
#include "fclaw2d_partition.h"
#include <fclaw2d_vtable.h>
#include <fclaw2d_clawpatch.hpp>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


void fclaw2d_regrid_tag4refinement(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *this_patch,
                                   int this_block_idx,
                                   int this_patch_idx,
                                   void *user)
{
    fclaw2d_vtable_t vt;
    int refine_patch, maxlevel, level;
    const amr_options_t* gparms;

    int domain_init = *((int*) user);

    vt = fclaw2d_get_vtable(domain);
    gparms = get_domain_parms(domain);

    maxlevel = gparms->maxlevel;
    level = this_patch->level;

    if (level < maxlevel)
    {
        refine_patch  =
            vt.patch_tag4refinement(domain,this_patch,this_block_idx,
                                    this_patch_idx, domain_init);
        if (refine_patch == 1)
        {
            fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
        }
    }
}

/* Tag family for coarsening */
static
void cb_tag4coarsening(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *fine_patches,
                       int blockno, int fine0_patchno,
                       void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(domain);

    int minlevel = gparms->minlevel;

    int level = fine_patches[0].level;

    if (level > minlevel)
    {
        int family_coarsened = 1;
        family_coarsened = vt.patch_tag4coarsening(domain,&fine_patches[0],
                                                  blockno, fine0_patchno);
        if (family_coarsened == 1)
        {
            for (int igrid = 0; igrid < NumSiblings; igrid++)
            {
                int fine_patchno = fine0_patchno + igrid;
                fclaw2d_patch_mark_coarsen(domain,blockno, fine_patchno);
            }
        }
    }
}


/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */

void fclaw2d_regrid_repopulate(fclaw2d_domain_t * old_domain,
                               fclaw2d_patch_t * old_patch,
                               fclaw2d_domain_t * new_domain,
                               fclaw2d_patch_t * new_patch,
                               fclaw2d_patch_relation_t newsize,
                               int blockno,
                               int old_patchno,
                               int new_patchno,
                               void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(new_domain);

    int domain_init = *((int*) user);

    fclaw2d_domain_data_t *ddata_old = fclaw2d_domain_get_data (old_domain);
    fclaw2d_domain_data_t *ddata_new = fclaw2d_domain_get_data (new_domain);

    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        new_patch->user = old_patch->user;
        old_patch->user = NULL;
        ++ddata_old->count_delete_clawpatch;
        ++ddata_new->count_set_clawpatch;
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        fclaw2d_patch_t *fine_siblings = new_patch;
        fclaw2d_patch_t *coarse_patch = old_patch;
        for (int i = 0; i < NumSiblings; i++)
        {
            fclaw2d_patch_t *fine_patch = &fine_siblings[i];
            int fine_patchno = new_patchno + i;
            fclaw2d_patch_data_new(new_domain,fine_patch);
            fclaw2d_clawpatch_build_cb(new_domain,fine_patch,blockno,
                                       fine_patchno,(void*) NULL);
            if (domain_init)
            {
                vt.patch_initialize(new_domain,fine_patch,blockno,fine_patchno);
            }
        }

        if (!domain_init)
        {
            int coarse_patchno = old_patchno;
            int fine_patchno = new_patchno;

            vt.patch_interpolate2fine(new_domain,coarse_patch,fine_siblings,
                                      blockno,coarse_patchno,fine_patchno);
        }
        fclaw2d_patch_data_delete(old_domain,coarse_patch);
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        /* Old grids are the finer grids;  new grid is the coarsened grid */
        fclaw2d_patch_t *fine_siblings = old_patch;
        int fine_patchno = old_patchno;

        fclaw2d_patch_t *coarse_patch = new_patch;
        int coarse_patchno = new_patchno;
        fclaw2d_patch_data_new(new_domain,coarse_patch);
        fclaw2d_clawpatch_build_cb(new_domain,coarse_patch,blockno,
                                   coarse_patchno,(void*) NULL);

        if (domain_init)
        {
            fclaw_debugf("fclaw2d_regrid.cpp (repopulate): We shouldn't end up here\n");
            exit(0);
        }
        else
        {
            vt.patch_average2coarse(new_domain,fine_siblings,coarse_patch,
                                    blockno,coarse_patchno, fine_patchno);
            for(int i = 0; i < 4; i++)
            {
                fclaw2d_patch_t* fine_patch = &fine_siblings[i];
                fclaw2d_patch_data_delete(old_domain,fine_patch);
            }
        }
    }
    else
    {
        fclaw_global_essentialf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
}

void fclaw2d_regrid_new_domain_setup(fclaw2d_domain_t* old_domain,
                                     fclaw2d_domain_t* new_domain)
{
    const amr_options_t *gparms;
    double t;

    if (old_domain == NULL)
    {
        fclaw_global_infof("Building initial domain\n");
        t = 0;
        fclaw2d_domain_set_time(new_domain,t);

    }
    else
    {
        fclaw_global_infof("Rebuilding  domain\n");
        /* Allocate memory for user data types (but they don't get set) */
        fclaw2d_domain_data_new(new_domain);
        fclaw2d_domain_data_copy(old_domain,new_domain);
    }

    gparms = get_domain_parms(new_domain);

    fclaw2d_block_data_new(new_domain);

    int num = new_domain->num_blocks;
    for (int i = 0; i < num; i++)
    {
        fclaw2d_block_t *block = &new_domain->blocks[i];
        /* This will work for rectangular domains ... */
        fclaw2d_block_set_data(block,gparms->mthbc);
    }

    /* Set up the parallel ghost patch data structure. */
    fclaw_global_infof("  -- Setting up parallel ghost exchange ... \n");

    fclaw2d_partition_setup(new_domain);

    fclaw_global_infof("Done\n");
}

/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */
void fclaw2d_regrid(fclaw2d_domain_t **domain)
{
    fclaw2d_domain_data_t* ddata = fclaw2d_domain_get_data(*domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_REGRID]);

    fclaw_global_infof("Regridding domain\n");


    /* First determine which families should be coarsened. */
    fclaw2d_domain_iterate_families(*domain, cb_tag4coarsening,
                                    (void*) NULL);

    int domain_init = 0;
    fclaw2d_domain_iterate_patches(*domain, fclaw2d_regrid_tag4refinement,
                                   (void *) &domain_init);

    /* Rebuild domain if necessary */
    /* Will return be NULL if no refining was done */
    fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
    fclaw_bool have_new_refinement = new_domain != NULL;

    /* Domain data may go out of scope now. */
    ddata = NULL;

    if (have_new_refinement)
    {
        fclaw_global_infof(" -- Have new refinement\n");

        /* allocate memory for user patch data and user domain data in the new
           domain;  copy data from the old to new the domain. */
        fclaw2d_regrid_new_domain_setup(*domain, new_domain);

        /* Average to new coarse grids and interpolate to new fine grids */
        fclaw2d_domain_iterate_adapted(*domain, new_domain,
                                       fclaw2d_regrid_repopulate,
                                       (void *) &domain_init);

        /* free memory associated with old domain */
        fclaw2d_domain_reset(domain);
        *domain = new_domain;
        new_domain = NULL;

        /* We have a newly created mesh;  We need to get all of the ghost cells
           filled */

        fclaw2d_partition_domain(domain, -1);

        /* fclaw2d_ghost_update_all_levels (*domain,FCLAW2D_TIMER_REGRID); */
        int minlevel = (*domain)->global_minlevel;
        int maxlevel = (*domain)->global_maxlevel;
        int time_interp = 0;
        fclaw2d_ghost_update(*domain,minlevel,maxlevel,time_interp,FCLAW2D_TIMER_REGRID);
        fclaw_global_infof ("Global minlevel %d maxlevel %d\n",
                            minlevel,maxlevel);
    }

    /* Stop timer */
    ddata = fclaw2d_domain_get_data(*domain);
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_REGRID]);

    /* Count calls to this function */
    ++ddata->count_amr_regrid;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
