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
#include "amr_mol.H"


/* ---------------------------------------------
   Store data for use by the function that evaluates
   the right hand side needed by MOL solver.
   ---------------------------------------------*/
static fclaw2d_f_rhs_data_fort_t *s_f_rhs_data = NULL;


// ----------------------------------------------------------
// Vectorize all patch data at a given level
// ----------------------------------------------------------
static
void cb_vectorize_patch_data(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             void *user)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    double *curr_data = cp->current_data_ptr();

    fclaw_level_mol_data_t *mol_data = (fclaw_level_mol_data_t*) user;
    int size = mol_data->patch_size;
    int count = mol_data->count;
    mol_data->patches[count] = this_patch;
    mol_data->block_indices[count] = this_block_idx;
    mol_data->patch_indices[count] = this_patch_idx;

    // Copy patch data to  vector.
    for (int i = 0; i < size; i++)
    {
        mol_data->patch_data[count*size + i] = curr_data[i];
    }
    mol_data->count++;
}

static
fclaw_level_mol_data_t* vectorize_patch_data(fclaw2d_domain_t *domain,
                                               int a_level)
{
    int num_patches_at_level = num_patches(domain,a_level);


    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    fclaw_level_mol_data_t *mol_data  = new fclaw_level_mol_data_t;
    mol_data->patches = new fclaw2d_patch_t*[num_patches_at_level];
    int neqn = (mx+2*mbc)*(my+2*mbc)*meqn;
    mol_data->patch_size = neqn;
    mol_data->patch_data = new double[num_patches_at_level*neqn];
    mol_data->patch_indices = new int[num_patches_at_level];
    mol_data->block_indices = new int[num_patches_at_level];
    mol_data->level = a_level;

    // Store patches in newly created vector mol_data->patch_data
    mol_data->count = 0;
    fclaw2d_domain_iterate_level(domain, a_level,
                                 cb_vectorize_patch_data,
                                 (void *) mol_data);
    return mol_data;
}

// ----------------------------------------------------------
// Restore vectorized data to patches in tree.
// ----------------------------------------------------------
static
    void restore_patch_data(fclaw2d_domain_t *domain,
                            int a_level,
                            fclaw_level_mol_data_t *mol_data)
{
    int count = mol_data->count; // total number of patches.
    int patch_size = mol_data->patch_size;
    for (int i = 0; i < count; i++)
    {
        fclaw2d_patch_t *this_patch = mol_data->patches[i];
        ClawPatch *cp = get_clawpatch(this_patch);
        double *cp_data = cp->current_data_ptr();
        double *patch_data = &mol_data->patch_data_restore[i*patch_size];
        for (int j = 0; j < patch_size; j++)
        {
            cp_data[j] = patch_data[j];
        }
    }
}

// Evaluate the right hand side, after first restoring
// q to the tree, and getting boundary conditions.
void fclaw_mol_rhs(const double& t_inner,
                   double *q, double *rhs)
{
    // Copy mol data back to domain level.  Use statically defined variables.
    fclaw2d_domain_t *domain = s_f_rhs_data->domain;
    fclaw2d_level_time_data_t *time_data = s_f_rhs_data->time_data;
    fclaw_level_mol_data_t *mol_data = s_f_rhs_data->mol_data;


    // The q pointer may not be the same one that is currently in
    // mol_data->patch_data, so we restore the q that is passed in.
    mol_data->patch_data_restore = q;

    // Things that are needed from the above
    const amr_options_t *gparms = get_domain_parms(domain);

    // From time_data
    double t_level = time_data->t_level;
    double t_initial = time_data->t_initial;
    double dt_coarse = time_data->dt_coarse;
    int level = mol_data->level;

    int count = mol_data->count; // Total number of patches at this level.
    int size = mol_data->patch_size;

    restore_patch_data(domain,level,mol_data);

    level_exchange(domain,level);
    set_phys_bc(domain,level, t_inner);

    // Set ghost cell values at newly updated data on level 'level'.
    if (t_inner > t_level)
    {
        // level_exchange(domain,level);

        if (!time_data->is_coarsest && gparms->subcycle)
        {
            double alpha = (t_inner - t_initial)/dt_coarse;
            printf("%d %16.8e %16.8e %16.8e\n",level,t_initial,t_inner,alpha);
            exchange_with_coarse(domain, level, t_inner, alpha);
        }
    }

    for(int i = 0; i < count; i++)
    {
        fclaw2d_patch_t *this_patch = mol_data->patches[i];
        int this_patch_idx = mol_data->patch_indices[i];
        int this_block_idx = mol_data->block_indices[i];

        // This is the user defined routine.  The value of
        // this function pointer is set in main.
        fclaw2d_domain_data_t *ddata = get_domain_data(domain);
        fclaw_mol_rhs_patch_t f_mol_rhs_patch_ptr = ddata->f_mol_rhs_patch;
        f_mol_rhs_patch_ptr(domain,this_patch,this_block_idx,
                            this_patch_idx,t_inner,&rhs[i*size]);
    }
}

/* ---------------------------------------------------------------
   Take a parabolic step of 'dt'.
--------------------------------------------------------------- */
void fclaw_mol_step(fclaw2d_domain_t *domain,
                    int level,
                    fclaw2d_level_time_data_t *time_data)
{
    double t_level = time_data->t_level;
    double dt = time_data->dt;

    fclaw_level_mol_data_t *mol_data = vectorize_patch_data(domain,level);

    int count = mol_data->count;
    int size = mol_data->patch_size;
    int neqn = size*count;

    // Set static data that is needed by 'f_exp' below.
    s_f_rhs_data = new fclaw2d_f_rhs_data_fort_t;
    s_f_rhs_data->domain = domain;
    s_f_rhs_data->mol_data = mol_data;
    s_f_rhs_data->time_data = time_data;

    // Call time step method and solve to time t.
    // The function below is an interface to the
    // actual solver.
    double *q = mol_data->patch_data;


    fclaw2d_domain_data *ddata = get_domain_data(domain);
    fclaw_mol_solver_t f_mol_solver_ptr = ddata->f_mol_solver;
    f_mol_solver_ptr(neqn,q,t_level,dt);

    // Return patch data to tree
    restore_patch_data(domain,level,mol_data);

    delete [] mol_data->patches;
    delete [] mol_data->patch_data;

    delete s_f_rhs_data;
}
