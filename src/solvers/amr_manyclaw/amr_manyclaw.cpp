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
#include "amr_manyclaw.H"

#include <manyclaw/manyclaw.h>

void set_manyclaw_parms(fclaw2d_domain_t* domain,amr_manyclaw_parms_t* manyclaw_parms)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    ddata->manyclaw_parms = (void*) manyclaw_parms;
}

amr_manyclaw_parms_t* get_manyclaw_parms(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    amr_manyclaw_parms_t *manyclaw_parms = (amr_manyclaw_parms_t*) ddata->manyclaw_parms;
    return manyclaw_parms;
}

static
    void reorder2new(fclaw2d_domain_t *domain, double *qin, double *qout)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    manyclaw_reorder2new_(mx,my,mbc, meqn,qin,qout);
}

static
    void reorder2old(fclaw2d_domain_t *domain, double *qin, double *qout)
{
    const amr_options_t *gparms = get_domain_parms(domain);

    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    manyclaw_reorder2old_(mx,my,mbc, meqn,qin,qout);
}



static
amr_manyclaw_patch_data_t* get_manyclaw_patch_data(ClawPatch *cp)
{
    /* We need the this cast here, because ClawPatch::manyclaw_patch_data() only returns a
       void* object */
    amr_manyclaw_patch_data_t *wp = (amr_manyclaw_patch_data_t*) cp->manyclaw_patch_data();
    return wp;
}

void amr_manyclaw_setprob(fclaw2d_domain_t* domain)
{
    setprob_();
}

void manyclaw_set_riemann_solvers(fclaw2d_patch_t *this_patch,
                                  rp_grid_eval_t rp_grid_eval, 
				  void* aux_global,
				  updater_t update)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    amr_manyclaw_patch_data_t *mc_data = get_manyclaw_patch_data(cp);
    mc_data->rp_grid_eval = rp_grid_eval;
    mc_data->aux_global = aux_global;
    mc_data->update = update;
}

void manyclaw_set_solver(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    amr_manyclaw_parms_t *manyclaw_parms = get_manyclaw_parms(domain);
    ClawPatch *cp = get_clawpatch(this_patch);
    amr_manyclaw_patch_data_t *mc_data = get_manyclaw_patch_data(cp);
    Solver &solver = mc_data->solver;

    /* Define the 'Solver' object */
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int mwaves = manyclaw_parms->mwaves;
    int numcells[2] = {mx, my};
    solver.define(numcells,meqn,mbc,mwaves);
}



/* This should only be called when a new ClawPatch is created. */
void amr_manyclaw_setaux(fclaw2d_domain_t *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    amr_manyclaw_parms_t *manyclaw_parms = get_manyclaw_parms(domain);

    ClawPatch *cp = get_clawpatch(this_patch);
    amr_manyclaw_patch_data_t *manyclaw_patch_data = get_manyclaw_patch_data(cp);

    /* For multiblock problems; otherwise block_idx == 0*/
    set_block_(&this_block_idx);

    int mx = gparms->mx;
    int my = gparms->my;
    int maxmx = mx;
    int maxmy = my;
    int mbc = gparms->mbc;

    /* Allocate aux array for each new patch */
    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-mbc;
    }
    ur[0] = mx + mbc;
    ur[1] = my + mbc;
    Box box(ll,ur);

    FArrayBox &auxarray = manyclaw_patch_data->auxarray; /* No space yet allocated */
    int maux = manyclaw_parms->maux;

    auxarray.define(box,maux);
    manyclaw_patch_data->maux = maux;

    double *aux = auxarray.dataPtr();

    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* Standard clawpack aux routine */
    setaux_(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,maux, aux);
}

void amr_manyclaw_qinit(fclaw2d_domain_t *domain,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx)
{
    const amr_options_t *gparms              = get_domain_parms(domain);
    ClawPatch *cp                            = get_clawpatch(this_patch);
    amr_manyclaw_parms_t *manyclaw_parms     = get_manyclaw_parms(domain);

    amr_manyclaw_patch_data_t *manyclaw_patch_data = get_manyclaw_patch_data(cp);

    set_block_(&this_block_idx);

    double* q = cp->q();

    double* aux = manyclaw_patch_data->auxarray.dataPtr();

    /* maux is also stored in manyclaw_patch_data */
    int maux = manyclaw_parms->maux;

    int mx = gparms->mx;
    int my = gparms->my;
    int maxmx = mx;
    int maxmy = my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    qinit_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux);
}

void amr_manyclaw_b4step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    const amr_options_t *gparms                    = get_domain_parms(domain);
    ClawPatch *cp                            = get_clawpatch(this_patch);
    amr_manyclaw_patch_data_t *manyclaw_patch_data = get_manyclaw_patch_data(cp);
    amr_manyclaw_parms_t *manyclaw_parms     = get_manyclaw_parms(domain);

    set_block_(&this_block_idx);

    double* q = cp->q();

    double* aux = manyclaw_patch_data->auxarray.dataPtr();
    int maux = manyclaw_parms->maux;

    int mx = gparms->mx;
    int my = gparms->my;
    int maxmx = mx;
    int maxmy = my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    b4step2_(maxmx,maxmy,mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux);
}

void amr_manyclaw_src2(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    amr_manyclaw_parms_t *manyclaw_parms = get_manyclaw_parms(domain);
    ClawPatch *cp = get_clawpatch(this_patch);

    /* Solver defined data;  see amr_manyclaw.H */
    amr_manyclaw_patch_data_t *manyclaw_patch_data = get_manyclaw_patch_data(cp);

    set_block_(&this_block_idx);

    double* q = cp->q();

    double* aux = manyclaw_patch_data->auxarray.dataPtr();
    int maux = manyclaw_parms->maux;

    int mx = gparms->mx;
    int my = gparms->my;
    int maxmx = mx;
    int maxmy = my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    src2_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt);
}


/* Use this to return only the right hand side of the manyclaw algorithm */
double amr_manyclaw_step2_rhs(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double *rhs)
{
    /* This should evaluate the right hand side, but not actually do the update.
       This will be useful in cases where we want to use something other than
       a single step method.  For example, in a RK scheme, one might want to
       call the right hand side to evaluate stages. */
    return 0;
}


void amr_manyclaw_bc2(fclaw2d_domain *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      double t,
                      double dt,
                      fclaw_bool intersects_phys_bdry[],
                      fclaw_bool time_interp)
{
    const amr_options_t* gparms              = get_domain_parms(domain);
    ClawPatch *cp                            = get_clawpatch(this_patch);
    amr_manyclaw_parms_t *manyclaw_parms     = get_manyclaw_parms(domain);
    amr_manyclaw_patch_data_t *manyclaw_patch_data = get_manyclaw_patch_data(cp);

    set_block_(&this_block_idx);

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    fclaw2d_block_data_t *bdata = get_block_data(this_block);
    int *block_mthbc = bdata->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[NumFaces];
    for(int i = 0; i < NumFaces; i++)
    {
        if (intersects_phys_bdry[i])
        {
            mthbc[i] = block_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }

    double* q = cp->q_time_sync(time_interp);

    int maux = manyclaw_parms->maux;
    double *aux = manyclaw_patch_data->auxarray.dataPtr();

    /* Global to all patches */
    int mx = gparms->mx;
    int my = gparms->my;
    int maxmx = mx;
    int maxmy = my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* Specific to the patch */
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* Modify for ManyClaw ordering */
    bc2_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc);
}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double amr_manyclaw_step2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    const amr_options_t* gparms              = get_domain_parms(domain);
    ClawPatch *cp                            = get_clawpatch(this_patch);
    amr_manyclaw_patch_data_t *mc_data = get_manyclaw_patch_data(cp);
    amr_manyclaw_parms_t * manyclaw_parms   = get_manyclaw_parms(domain);

    set_block_(&this_block_idx);

    // int level = this_patch->level;

    cp->save_current_step();  /* Save in case we need to retake the step */

    double* aux = mc_data->auxarray.dataPtr();
    void* aux_global  = mc_data->aux_global;
    // int maux = manyclaw_parms->maux;

    // Global to all patches
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    // Specific to the patch
    // double xlower = cp->xlower();
    // double ylower = cp->ylower();
    double dx = cp->dx();
    // double dy = cp->dy();

    // Specific to solver
    int mwaves = manyclaw_parms->mwaves;

    /* ManyClaw update :
       This takes a step on a single grid and updates the grid with the
       & new solution
    */

    Solver &solver = mc_data->solver;
    double dtdx = dt/dx;

    double* qold_mlast = cp->q();
    FArrayBox qold_mfirst_array = cp->newGrid(); /* Create a scratch grid */
    double *qold_mfirst = qold_mfirst_array.dataPtr();
    reorder2new(domain,qold_mlast,qold_mfirst);

    mc_data->rp_grid_eval(&qold_mfirst[0], &aux[0], aux_global,
                          mx, my,
                          &solver.amdq[0],
                          &solver.apdq[0],
                          &solver.wave[0],
                          &solver.wave_speed[0]);

    mc_data->update(&qold_mfirst[0], &aux[0],
                    mx, my,
                    &solver.amdq[0],
                    &solver.apdq[0],
                    &solver.wave[0],
                    &solver.wave_speed[0],
                    mbc, meqn, dtdx);

    reorder2old(domain,qold_mfirst,qold_mlast);

    double cflgrid = solver.cfl(mx, my, mbc, meqn, mwaves, dtdx);


    return cflgrid;
}

/* This routine combines a few of the steps that one might do to update a patch */
double amr_manyclaw_update(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           double t,
                           double dt)
{
    const amr_manyclaw_parms_t* manyclaw_parms = get_manyclaw_parms(domain);

    amr_manyclaw_b4step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    double maxcfl = amr_manyclaw_step2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);

    if (manyclaw_parms->src_term > 0)
    {
        amr_manyclaw_src2(domain,this_patch,this_block_idx,this_patch_idx,t,dt);
    }
    return maxcfl;
}


amr_manyclaw_parms_t*  amr_manyclaw_parms_new(sc_options_t *opt)
{
    /* This is a good example of a place where a [Section] would be nice.
       Something like [manyclaw].  See fclaw_defaults.ini */

    amr_manyclaw_parms_t *manyclaw_parms;

    manyclaw_parms = FCLAW2D_ALLOC_ZERO (amr_manyclaw_parms_t, 1);

    /* Array of SpaceDim many values, with no defaults is set to all 0's */
    amr_options_add_int_array (opt, 0, "order", &manyclaw_parms->order_string, NULL,
                               &manyclaw_parms->order, SpaceDim,
                               "Normal and transverse orders");

    sc_options_add_int (opt, 'M', "mcapa", &manyclaw_parms->mcapa, -1,
                        "Location of capacity function in aux array [-1]");

    sc_options_add_int (opt, 0, "maux", &manyclaw_parms->maux, 0,
                        "Number of auxiliary variables [0]");

    sc_options_add_int (opt, 0, "src_term", &manyclaw_parms->src_term, 0,
                        "Source term option [0]");

    sc_options_add_int (opt, 0, "mwaves", &manyclaw_parms->mwaves, 1,
                        "Number of waves [1]");

    /* Array of mwaves many values */
    amr_options_add_int_array (opt, 0, "mthlim", &manyclaw_parms->mthlim_string, NULL,
                               &manyclaw_parms->mthlim, manyclaw_parms->mwaves,
                               "Waves limiters (one for each wave)");
    /*
     * At this point amropt->mthlim is allocated, but for precisely 1 wave,
     * since that is the default value provided to manyclaw_parms->waves above.
     * If the number of waves (manyclaw_parms->mwaves) is changed by option
     * parsing, then this array needs to be converted anew. This can be done
     * by calling amr_manyclaw_postprocess_parms.
     */


    /* -----------------------------------------------------------------------
       Read in options from file
       ----------------------------------------------------------------------- */
    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-U' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.

       WARNING: The long option name must be unique within the whole program.
                Just 'inifile' is already used in amr_options.c.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'U', "manyclaw_inifile",
                            "Read ManyClaw options from this file [default : fclaw2d_manyclaw.ini]");

    /* -----------------------------------------------------------------------
       This is the default file that will be read if no command line options are
       given.  This file is tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_load (sc_package_id, SC_LP_ALWAYS, opt, "fclaw2d_manyclaw.ini");

    amr_manyclaw_postprocess_parms(manyclaw_parms);

    return manyclaw_parms;

}

void amr_manyclaw_checkparms(amr_manyclaw_parms_t* manyclaw_parms, amr_options_t* gparms)
{
    /* -----------------------------------------------------------------------
       Set up 'method' vector used by Clawpack.
       ------------------------------------------------------------------------ */
    manyclaw_parms->method[0] = gparms->use_fixed_dt;

    manyclaw_parms->method[1] = manyclaw_parms->order[0];
    if (SpaceDim == 2)
    {
        manyclaw_parms->method[2] = manyclaw_parms->order[1];
    }
    else
    {
        manyclaw_parms->method[2] = 10*manyclaw_parms->order[1] + manyclaw_parms->order[2];
    }
    manyclaw_parms->method[3] = gparms->verbosity;
    manyclaw_parms->method[4] = manyclaw_parms->src_term;
    manyclaw_parms->method[5] = manyclaw_parms->mcapa;
    manyclaw_parms->method[6] = manyclaw_parms->maux;

    /* Should also check mthbc, mthlim, etc. */

}

void amr_manyclaw_postprocess_parms(amr_manyclaw_parms_t* manyclaw_parms)
{
    /* -----------------------------------------------------------------------
       Some post-processing of arrays
       ------------------------------------------------------------------------ */

    amr_options_convert_int_array (manyclaw_parms->mthlim_string, &manyclaw_parms->mthlim,
                                   manyclaw_parms->mwaves);

    amr_options_convert_int_array (manyclaw_parms->order_string, &manyclaw_parms->order,
                                   SpaceDim);
}


void amr_manyclaw_parms_delete(amr_manyclaw_parms_t* manyclaw_parms)
{
    SC_FREE(manyclaw_parms->order);
    SC_FREE(manyclaw_parms->mthlim);
    FCLAW2D_FREE(manyclaw_parms);
}

static
void amr_manyclaw_patch_data_new(void** wp)
{
    amr_manyclaw_patch_data_t* manyclaw_patch_data;
    manyclaw_patch_data = new amr_manyclaw_patch_data_t;
    *wp = (void*) manyclaw_patch_data;

    // or?
    // *wp = (void*) FCLAW2D_ALLOC_ZERO (amr_manyclaw_patch_data_t, 1);
}

static
void amr_manyclaw_patch_data_delete(void **wp)
{
    amr_manyclaw_patch_data_t *manyclaw_patch_data = (amr_manyclaw_patch_data_t*) *wp;
    delete manyclaw_patch_data->aux_global;
    delete manyclaw_patch_data;

    // or?
    // FCLAW2D_FREE(manyclaw_patch_data);

    *wp = (void*) NULL;
}


void amr_manyclaw_link_to_clawpatch()
{
    /* These are called whenever a new ClawPatch is created. */
    ClawPatch::f_manyclaw_patch_data_new = &amr_manyclaw_patch_data_new;
    ClawPatch::f_manyclaw_patch_data_delete = &amr_manyclaw_patch_data_delete;
}


/* This the default "linker".   This can be used if the user only wants default
   routines listed below.  */
void  amr_manyclaw_link_solvers(fclaw2d_domain_t* domain)
{
    const amr_manyclaw_parms_t* manyclaw_parms = get_manyclaw_parms(domain);

    fclaw2d_solver_functions_t* sf = get_solver_functions(domain);

    sf->use_single_step_update = fclaw_true;
    sf->use_mol_update = fclaw_false;

    if (manyclaw_parms->maux > 0)
    {
        sf->f_patch_setup          = &amr_manyclaw_setaux;
    }
    else
    {
        sf->f_patch_setup          = &amr_dummy_patch_setup;
    }

    sf->f_patch_initialize         = &amr_manyclaw_qinit;
    sf->f_patch_physical_bc        = &amr_manyclaw_bc2;

    /* Calls b4step2, step2 and src2 */
    sf->f_patch_single_step_update = &amr_manyclaw_update;

    /* This is needed so that a ClawPatch knows how to create data for a manyclaw solver */
    amr_manyclaw_link_to_clawpatch();
}
