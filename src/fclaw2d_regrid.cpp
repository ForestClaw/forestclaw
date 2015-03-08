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

#include <forestclaw2d.h>
#include <amr_utils.H>
#include <fclaw2d_vtable.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


static
void cb_tag4refinement(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    fclaw2d_vtable_t vt;
    int initflag, refine_patch, maxlevel, level;
    const amr_options_t* gparms;

    vt = fclaw2d_get_vtable(domain);
    gparms = get_domain_parms(domain);

    maxlevel = gparms->maxlevel;
    level = this_patch->level;

    /* The initial mesh is built in a separate call (see amrinit) */
    initflag = 0;

    if (level < maxlevel)
    {
        refine_patch  =
            vt.patch_tag4refinement(domain,this_patch,this_block_idx,
                                    this_patch_idx, initflag);
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
        double xlower[4], ylower[4];
        int family_coarsened = 1;
        for (int igrid = 0; igrid < NumSiblings; igrid++)
        {
            xlower[igrid] = fine_patches[igrid].xlower;
            ylower[igrid] = fine_patches[igrid].ylower;
        }
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


static
void cb_rebuild_patches(fclaw2d_domain_t * old_domain,
                        fclaw2d_patch_t * old_patch,
                        fclaw2d_domain_t * new_domain,
                        fclaw2d_patch_t * new_patch,
                        fclaw2d_patch_relation_t newsize,
                        int blockno, int old_patchno,
                        int new_patchno,  void *user)
{
    fclaw2d_vtable_t vt;
    vt = fclaw2d_get_vtable(new_domain);

    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        /* We should be able to just pass patch pointer here */
        vt.patch_copy2samesize(new_domain,old_patch,new_patch,
                               blockno,old_patchno, new_patchno);
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        fclaw2d_patch_t *fine_siblings = new_patch;
        int fine_patchno = old_patchno;

        fclaw2d_patch_t *coarse_patch = old_patch;
        int coarse_patchno = old_patchno;

        vt.patch_interpolate2fine(new_domain,coarse_patch,fine_siblings,
                                  blockno,coarse_patchno,fine_patchno);
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        /* Old grids are the finer grids;  new grid is the coarsened grid */
        fclaw2d_patch_t *fine_siblings = old_patch;
        int fine_patchno = old_patchno;

        fclaw2d_patch_t *coarse_patch = new_patch;
        int coarse_patchno = new_patchno;

        vt.patch_average2coarse(new_domain,fine_siblings,coarse_patch,
                                blockno,coarse_patchno, fine_patchno);
    }
    else
    {
        fclaw_global_essentialf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
}


void fclaw2d_regrid(fclaw2d_domain_t **domain)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(*domain);
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_REGRID]);

    /* First determine which families should be coarsened. */
    fclaw2d_domain_iterate_families(*domain, cb_tag4coarsening,
                                    (void*) NULL);

    /* Then refine. */
    fclaw2d_domain_iterate_patches(*domain, cb_tag4refinement,
                                   (void *) NULL);

    /* Rebuild domain if necessary */
    /* Will return be NULL if no refining was done */
    fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
    fclaw_bool have_new_refinement = new_domain != NULL;

    /* Domain data may go out of scope now. */
    ddata = NULL;

    if (have_new_refinement)
    {
        /* allocate memory for user patch data and user domain data in the new
           domain;  copy data from the old to new the domain. */
        rebuild_domain(*domain, new_domain);

        /* Average to new coarse grids and interpolate to new fine grids */
        fclaw2d_domain_iterate_adapted(*domain, new_domain,cb_rebuild_patches,
                                       (void *) NULL);

        /* TODO: can we use global min/maxlevels here? */
        /* DAC : I am not sure - The coarse/fine exchanges will
           probably break if passed levels on which there are no
           patches.  This could of course be fixed, but it doesn't
           seem necessary at this point.
        */
        /* Do a level exchange to make sure all data is valid;
           Coarse and fine grid exchanges will happen during time stepping
           as needed. */
#if 0
        /* currently unused */
        minlevel = new_domain->global_minlevel;
        maxlevel = new_domain->global_maxlevel;
#endif

        /* free memory associated with old domain */
        amrreset(domain);
        *domain = new_domain;
        new_domain = NULL;

        repartition_domain(domain, -1);
        update_ghost_all_levels (*domain,FCLAW2D_TIMER_REGRID);
        fclaw_global_infof ("Global minlevel %d maxlevel %d\n",
                            (*domain)->global_minlevel, (*domain)->global_maxlevel);
    }

    /* Stop timer */
    ddata = get_domain_data(*domain);
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
