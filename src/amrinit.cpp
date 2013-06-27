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
#include "clawpack_fort.H"
#include "fclaw2d_solvers.H"

// Put this here so that I don't have to include "ClawPatch.H"
void set_clawpatch(fclaw2d_domain_t* domain, fclaw2d_patch_t *this_patch,
                   int blockno, int patchno);

// This is essentially the same function that is in amr_regrid.cpp
static
void cb_tag4refinement_init(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx,
                            void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    int maxlevel = gparms->maxlevel;
    int level = this_patch->level;
    int initflag = 1;

    if (level < maxlevel)
    {
        fclaw2d_regrid_functions_t* rf = get_regrid_functions(domain);
        fclaw_bool refine_patch =
            (rf->f_patch_tag4refinement)(domain,this_patch,this_block_idx,
                                         this_patch_idx,initflag);
        if (refine_patch)
        {
            fclaw2d_patch_mark_refine(domain, this_block_idx, this_patch_idx);
        }
    }
}

/* -----------------------------------------------------------------
   Initial grid
   ----------------------------------------------------------------- */
static
void cb_initialize (fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    void *user)
{
    set_clawpatch(domain,this_patch,this_block_idx, this_patch_idx);

    fclaw2d_solver_functions_t *sf = get_solver_functions(domain);

    // One-time setup of patch
    (sf->f_patch_setup)(domain,this_patch,this_block_idx,this_patch_idx);

    // Set initial values on patch
    (sf->f_patch_initialize)(domain,this_patch,this_block_idx,this_patch_idx);
}

static
void cb_domain_adapt_init (fclaw2d_domain_t * old_domain,
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
        // Patch has been refined.
        fclaw2d_patch_t *fine_siblings = new_patch;

        // Loop over four siblings (z-ordering) and interpolate
        // data from coarser parent grid to each finer sibling grid.
        for (int igrid = 0; igrid < NumSiblings; igrid++)
        {
            fclaw2d_patch_t *fine_patch = &fine_siblings[igrid];
            int fine_patchno = new_patchno + igrid;

            // Create new ClawPatch and assign patch pointer to it.
            set_clawpatch(new_domain, fine_patch, blockno, fine_patchno);

            // Do one-time setup on new patch
            fclaw2d_solver_functions_t *sf = get_solver_functions(old_domain);
            (sf->f_patch_setup)(new_domain,fine_patch,blockno,fine_patchno);

            // This is only used here, since only in the initial grid layout do we
            // create fine grids from coarser grids.
            (sf->f_patch_initialize)(new_domain,fine_patch,blockno,fine_patchno);
        }
    }
    else if (newsize == FCLAW2D_PATCH_DOUBLESIZE)
    {
        // We don't coarsen for the initial time step
    }
    else
    {
        printf("cb_adapt_domain : newsize not recognized\n");
        exit(1);
    }
}


// Initialize a base level of grids
void amrinit (fclaw2d_domain_t **domain)
{
    const amr_options_t *gparms = get_domain_parms(*domain);
    fclaw2d_domain_data_t* ddata = get_domain_data(*domain);

    // Set problem dependent parameters for Riemann solvers, etc.
    (ddata->f_problem_setup)(*domain);

    double t = 0;

    set_domain_time(*domain,t);

    int minlevel = gparms->minlevel;
    int maxlevel = gparms->maxlevel;

    init_block_and_patch_data(*domain);  /* Allocate block and patch data */

    // Initialize base level grid - combine with 'amr_set_base_level' above?
    fclaw2d_domain_iterate_level(*domain, minlevel, cb_initialize,
                                 (void *) NULL);

    int num = (*domain)->num_blocks;
    for (int i = 0; i < num; i++)
    {
        fclaw2d_block_t *block = &(*domain)->blocks[i];
        set_block_data(block,gparms->mthbc);
    }

    set_phys_bc(*domain,minlevel,t);

    // Refine as needed.

    for (int level = minlevel; level < maxlevel; level++)
    {
        // TODO: don't use level_refined since it is not agreed upon in parallel
        // the fclaw2d_domain_adapt and _partition calls work fine in parallel

        fclaw2d_domain_iterate_level(*domain, level,
                                     cb_tag4refinement_init,
                                     (void *) NULL);

        // Rebuild domain if necessary
        fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);

        if (new_domain != NULL)
        {
            // This is just for fun; remove when it gets annoying.
            // fclaw2d_domain_list_adapted(*domain, new_domain, SC_LP_STATISTICS);

            // Allocate memory for user data types (but they don't get set)
            init_domain_data(new_domain);
            copy_domain_data(*domain,new_domain);

            // Initialize new grids.  Assume that all ghost cells are filled
            //in by qinit.
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

            fclaw2d_domain_iterate_adapted(*domain, new_domain,
                                           cb_domain_adapt_init,
                                           (void *) NULL);

            int new_level = level+1;
            // Upon initialization, we don't do any ghost cell exchanges, because we assume
            // that the initial conditions have set all internal ghost cells.
            set_phys_bc(new_domain,new_level,t);

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
        else
        {
            // exit loop;  we are done refining
            break;
        }
    }
}
