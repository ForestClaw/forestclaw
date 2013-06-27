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

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw2d_typedefs.h"
#include "amr_regrid.H"

/* -----------------------------------------------------------------
   Callback routine for tagging
   ----------------------------------------------------------------- */
static
void cb_tag4refinement(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    int maxlevel = gparms->maxlevel;
    int level = this_patch->level;

    if (level < maxlevel)
    {
        fclaw2d_regrid_functions_t* rf = get_regrid_functions(domain);

        int initflag = 0;
        fclaw_bool refine_patch =
            (rf->f_patch_tag4refinement)(domain,this_patch,this_block_idx,
                                         this_patch_idx,initflag);
        if (refine_patch)
        {
            fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
        }
    }
}

/* Tag family for coarsening */
static
void cb_tag4coarsening(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *fine_patches,
                       int blockno,
                       int fine0_patchno,
                       void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    int minlevel = gparms->minlevel;

    int level = fine_patches[0].level;

    if (level > minlevel)
    {
        /* --------------------------------------------------------------
           Create temporary patch to hold coarsened data that we will use
           to determine if we need to coarsen data
          ----------------------------------------------------------------- */

        fclaw2d_patch_t *temp_coarse_patch = new fclaw2d_patch_t;
        init_patch_data(temp_coarse_patch);

        temp_coarse_patch->xlower = fine_patches[0].xlower;
        temp_coarse_patch->ylower = fine_patches[0].ylower;
        temp_coarse_patch->xupper = fine_patches[NumSiblings-1].xupper;
        temp_coarse_patch->yupper = fine_patches[NumSiblings-1].yupper;
        temp_coarse_patch->level  = level;
        int coarse_patchno = -1;

        set_clawpatch(domain,temp_coarse_patch,blockno,coarse_patchno);

        // One-time setup of patch
        fclaw2d_solver_functions_t *sf = get_solver_functions(domain);
        (sf->f_patch_setup)(domain,temp_coarse_patch,blockno,coarse_patchno);

        fclaw2d_regrid_functions_t *rf = get_regrid_functions(domain);
        (rf->f_patch_average2coarse)(domain,fine_patches,temp_coarse_patch,
                                     blockno,coarse_patchno, fine0_patchno);

        /* --------------------------------------------------------------
           Test to see if temporary patch needs refining.  If so, then we
           shouldn't coarsen it.
          ----------------------------------------------------------------- */
        fclaw_bool patch_coarsened =
            (rf->f_patch_tag4coarsening)(domain, temp_coarse_patch, blockno,
                                         coarse_patchno);

        // not deleting ClawPatch - need to fix this!
        delete_clawpatch(temp_coarse_patch);
        delete temp_coarse_patch;

        if (patch_coarsened)
        {
            for (int igrid = 0; igrid < NumSiblings; igrid++)
            {
                int fine_patchno = fine0_patchno + igrid;
                fclaw2d_patch_mark_coarsen(domain, blockno, fine_patchno);
            }
        }
    }
}


static
void cb_domain_adapt(fclaw2d_domain_t * old_domain,
                     fclaw2d_patch_t * old_patch,
                     fclaw2d_domain_t * new_domain,
                     fclaw2d_patch_t * new_patch,
                     fclaw2d_patch_relation_t newsize,
                     int blockno,
                     int old_patchno, int new_patchno,
                     void *user)
{
    if (newsize == FCLAW2D_PATCH_SAMESIZE)
    {
        // Grid doesn't change
        set_clawpatch(new_domain,new_patch,blockno,new_patchno);

        // Setup new patch using solver specific routine
        fclaw2d_solver_functions_t *sf = get_solver_functions(old_domain);
        (sf->f_patch_setup)(new_domain,new_patch,blockno,new_patchno);

        // Need a copy function in regrid_functions
        fclaw2d_regrid_functions_t *rf = get_regrid_functions(old_domain);
        (rf->f_patch_copy2samesize)(new_domain,old_patch,new_patch,blockno,old_patchno,
                                    new_patchno);
    }
    else if (newsize == FCLAW2D_PATCH_HALFSIZE)
    {
        // Old grid is the coarse patch; new grids are the finer patches

        fclaw2d_patch_t *coarse_patch = old_patch;
        int coarse_patchno = old_patchno;

        fclaw2d_patch_t *fine_siblings = new_patch;

        // Loop over four siblings (z-ordering)
        for (int igrid = 0; igrid < NumSiblings; igrid++)
        {
            fclaw2d_patch_t *fine_patch = &fine_siblings[igrid];
            int fine_patchno = new_patchno + igrid;

            // Create new ClawPatch and assign patch pointer to it.
            set_clawpatch(new_domain, fine_patch, blockno, fine_patchno);

            // Do one-time setup on new patch
            fclaw2d_solver_functions_t *sf = get_solver_functions(old_domain);
            (sf->f_patch_setup)(new_domain,fine_patch,blockno,fine_patchno);

            // Initialize new fine patch by either calling an init function or
            // by interpolating from coarser grid.
            fclaw2d_regrid_functions_t *rf = get_regrid_functions(old_domain);
            (rf->f_patch_interpolate2fine)(new_domain,coarse_patch,fine_patch,
                                           blockno,coarse_patchno,fine_patchno,igrid);
        }
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        // Old grids are the finer grids;  new grid is the coarsened grid
        fclaw2d_patch_t *fine_siblings = old_patch;
        int fine_patchno = old_patchno;

        fclaw2d_patch_t *coarse_patch = new_patch;
        int coarse_patchno = new_patchno;

        set_clawpatch(new_domain,coarse_patch,blockno,coarse_patchno);

        fclaw2d_solver_functions_t *sf = get_solver_functions(old_domain);
        (sf->f_patch_setup)(new_domain,coarse_patch,blockno,coarse_patchno);

        fclaw2d_regrid_functions_t *rf = get_regrid_functions(old_domain);
        (rf->f_patch_average2coarse)(new_domain,fine_siblings,coarse_patch,
                                     blockno,coarse_patchno, fine_patchno);
    }
    else
    {
        printf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
}


void regrid(fclaw2d_domain_t **domain)
{

    const amr_options_t *gparms = get_domain_parms(*domain);
    double t = get_domain_time(*domain);

    int minlevel = gparms->minlevel;
    int maxlevel = gparms->maxlevel;

    for(int level = maxlevel; level > minlevel; level--)
    {
        double alpha = 0;
        exchange_with_coarse(*domain,level,t,alpha);
    }

    // First determine which families should be coarsened.
    fclaw2d_domain_iterate_families(*domain, cb_tag4coarsening,
                                    (void*) NULL);

    // Then refine.
    fclaw2d_domain_iterate_patches(*domain, cb_tag4refinement,
                                   (void *) NULL);

    // Rebuild domain if necessary
    // Will return be NULL if no refining was done?
    fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);

    if (new_domain != NULL)
    {
        // This is just for fun; remove when it gets annoying.
        // fclaw2d_domain_list_adapted(*domain, new_domain, SC_LP_STATISTICS);

        // Allocate memory for user data types (but they don't get set)
        init_domain_data(new_domain);
        copy_domain_data(*domain,new_domain);
        set_domain_time(new_domain,t);

        // Average or interpolate to new grids.
        init_block_and_patch_data(new_domain);

        // Physical BCs are needed in boundary level exchange
        // Assume only one block, since we are assuming mthbc
        int num = new_domain->num_blocks;
        for (int i = 0; i < num; i++)
        {
            fclaw2d_block_t *block = &new_domain->blocks[i];
            // This is kind of dumb for now, since block won't in general
            // have the same physical boundary conditions types.
            set_block_data(block,gparms->mthbc);
        }

        fclaw2d_domain_iterate_adapted(*domain, new_domain,cb_domain_adapt,
                                       (void *) NULL);

        // Set some of the user data types.  Some of this is done in
        // 'amr_set_base_level',
        // I should probably come up with a more general way to do this.
        // set_domain_data(new_domain, gparms);

        // Level stuff to make sure all
        for (int level = minlevel; level <= maxlevel; level++)
        {
            // Only do a level exchange;  Coarse and fine grid exchanges will happen
            // during time stepping as needed.
            level_exchange(new_domain,level);
            set_phys_bc(new_domain,level,t);
        }

        // free all memory associated with old domain
        amrreset(domain);
        *domain = new_domain;

        fclaw2d_domain_t *domain_partitioned =
            fclaw2d_domain_partition (*domain);

        if (domain_partitioned != NULL)
        {
            // TODO: allocate patch and block etc. memory for domain_partitioned

            // TODO: write a function to transfer values in parallel */

            /* then the old domain is no longer necessary */
            amrreset(domain);
            *domain = domain_partitioned;

            /* internal clean up */
            fclaw2d_domain_complete(*domain);
        }
    }
}
