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
#include "amr_solver_typedefs.H"

#include "ClawPatch.H"

void cb_tag4refinement(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       void *user);

void cb_domain_adapt(fclaw2d_domain_t * old_domain,
                     fclaw2d_patch_t * old_patch,
                     fclaw2d_domain_t * new_domain,
                     fclaw2d_patch_t * new_patch,
                     fclaw2d_patch_relation_t newsize,
                     int blockno,
                     int old_patchno, int new_patchno,
                     void *user);


static
void cb_init_base(fclaw2d_domain_t *domain,
                  fclaw2d_patch_t *this_patch,
                  int this_block_idx,
                  int this_patch_idx,
                  void *user)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    ClawPatch *cp = new ClawPatch();

    cp->define(this_patch->xlower,
               this_patch->ylower,
               this_patch->xupper,
               this_patch->yupper,
               this_block_idx,
               gparms);
    set_patch_data(this_patch,cp);

    /* The user can now retrieve the ClawPatch from 'this_patch' and set
       up whatever they need to set up. */
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    (ddata->f_patch_setup_ptr)(domain,this_patch,this_block_idx,this_patch_idx);
}

static
void set_base_level(fclaw2d_domain_t *domain, const int& level)
{
    fclaw2d_domain_iterate_level(domain, level, cb_init_base,(void *) NULL);
}

/* -----------------------------------------------------------------
   Initial grid
   ----------------------------------------------------------------- */
static
void cb_amrinit(fclaw2d_domain_t *domain,
                fclaw2d_patch_t *this_patch,
                int this_block_idx,
                int this_patch_idx,
                void *user)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    (ddata->f_patch_initialize_ptr)(domain,this_patch,this_block_idx,this_patch_idx);
}

// Initialize a base level of grids
void amrinit(fclaw2d_domain_t **domain)
{
    const amr_options_t *gparms = get_domain_parms(*domain);
    double t = 0;

    set_domain_time(*domain,t);

    int minlevel = gparms->minlevel;
    int maxlevel = gparms->maxlevel;

    // Set problem dependent parameters for Riemann solvers, etc.
    // Values are typically stored in Fortran common blocks, and are not
    // available outside of Fortran.
    set_problem_parameters();

    // Set up storage for base level grids so we can initialize them
    // Allocates per-block and per-patch user data

    // This function is redundant, and should be made more general.
    cout << "Setting base level " << endl;
    set_base_level(*domain,minlevel);

    cout << "Done with amr_set_base_level " << endl;


    // Initialize base level grid - combine with 'amr_set_base_level' above?
    fclaw2d_domain_iterate_level(*domain, minlevel, cb_amrinit,
                                 (void *) NULL);

    cout << "Done with domain adaptation " << endl;

    int num = (*domain)->num_blocks;
    for (int i = 0; i < num; i++)
    {
        fclaw2d_block_t *block = &(*domain)->blocks[i];
        set_block_data(block,gparms->mthbc);
    }

    set_phys_bc(*domain,minlevel,t);

    // Refine as needed.

    bool init_flag = true;
    for (int level = minlevel; level < maxlevel; level++)
    {
        cout << "amrinit : level = " << level << endl << endl;

        // TODO: don't use level_refined since it is not agreed upon in parallel
        // the fclaw2d_domain_adapt and _partition calls work fine in parallel

        fclaw2d_domain_iterate_level(*domain, level,
                                     cb_tag4refinement,
                                     (void *) &init_flag);

        // Rebuild domain if necessary
        cout << "amrinit : Building new domain " << endl;
        fclaw2d_domain_t *new_domain = fclaw2d_domain_adapt(*domain);
        cout << "amrinit : Done building new domain " << endl << endl;

        if (new_domain != NULL)
        {
            // This is just for fun; remove when it gets annoying.
            // fclaw2d_domain_list_adapted(*domain, new_domain, SC_LP_STATISTICS);

            // Allocate memory for user data types (but they don't get set)
            allocate_user_data(new_domain);
            copy_domain_data(*domain,new_domain);

            // Initialize new grids.  Assume that all ghost cells are filled
            //in by qinit.
            fclaw2d_domain_iterate_adapted(*domain, new_domain,
                                           cb_domain_adapt,
                                           (void *) &init_flag);

            // Set some of the user data types.  Some of this is done
            // in 'amr_set_base_level',
            // I should probably come up with a more general way to do this.

            // Not needed, because of copy above.
            // set_domain_data(new_domain, gparms);
            // set_domain_time(new_domain,t);

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
        cout << "amrinit : done with level " << endl << endl;;
    }
    cout << "Done with building initial grid structure " << endl;
}
