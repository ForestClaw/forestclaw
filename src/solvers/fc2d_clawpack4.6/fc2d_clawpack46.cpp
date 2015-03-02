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

#include "amr_includes.H"

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw_options.h"
#include "clawpack_fort.H"
#include "fclaw2d_clawpatch.h"

#include "fc2d_clawpack46_options.h"
#include "fc2d_clawpack46.H"

#include <fclaw_package.h>

static int s_clawpack46_package_id = -1;

static fc2d_clawpack46_vtable_t classic_vt;



void fc2d_clawpack46_set_vtable(const fc2d_clawpack46_vtable_t* user_vt)
{
    classic_vt = *user_vt;
}

void fc2d_clawpack46_init_vtable(fc2d_clawpack46_vtable_t* vt)
{
    /* Required functions  - error if NULL*/
    vt->bc2 = CLAWPACK46_BC2;
    vt->qinit = NULL;
    vt->rpn2 = NULL;
    vt->rpt2 = NULL;

    /* Optional functions - call only if non-NULL */
    vt->setprob = NULL;
    vt->setaux = NULL;
    vt->b4step2 = NULL;
    vt->src2 = NULL;
}





/* Patch data is stored in each ClawPatch */
struct patch_aux_data
{
    FArrayBox auxarray;
    int maux;
};

fc2d_clawpack46_options_t* fc2d_clawpack46_get_options(fclaw2d_domain_t *domain)
{
    int id;
    fclaw_app_t* app;
    app = fclaw2d_domain_get_app(domain);
    id = fc2d_clawpack46_get_package_id();
    return (fc2d_clawpack46_options_t*) fclaw_package_get_options(app,id);
}

static patch_aux_data*
get_patch_data(ClawPatch *cp)
{
    patch_aux_data *wp =
        (patch_aux_data*) cp->clawpack_patch_data(s_clawpack46_package_id);
    return wp;
}

static void*
patch_data_new()
{
    patch_aux_data* data;
    data = new patch_aux_data;
    return (void*) data;
}

static void
patch_data_delete(void *data)
{
    patch_aux_data *pd = (patch_aux_data*) data;
    FCLAW_ASSERT(pd != NULL);
    delete pd;
}

static const fclaw_package_vtable_t clawpack46_patch_vtable = {
    patch_data_new,
    patch_data_delete
};


/* -----------------------------------------------------------
   Public interface to routines in this file
   ----------------------------------------------------------- */
#if 0
void fc2d_clawpack46_package_register(fclaw_app_t* app,
                                      fc2d_clawpack46_options_t *clawopt)
{
    int id;

    /* Don't register a package more than once */
    FCLAW_ASSERT(s_clawpack46_package_id == -1);

    /* Register packages */
    id = fclaw_package_container_add_pkg(app,clawopt,
                                         &clawpack46_patch_vtable);
    s_clawpack46_package_id = id;
}
#endif


void fc2d_clawpack46_register(fclaw_app_t* app, const char *configfile)
{
    fc2d_clawpack46_options_t* clawopt;
    int id;

    /* Register the options */
    clawopt = fc2d_clawpack46_options_register(app,configfile);

    /* And the package */

    /* Don't register a package more than once */
    FCLAW_ASSERT(s_clawpack46_package_id == -1);

    id = fclaw_package_container_add_pkg(app,clawopt,
                                         &clawpack46_patch_vtable);

    s_clawpack46_package_id = id;
}

int fc2d_clawpack46_get_package_id()
{
    return s_clawpack46_package_id;
}


void fc2d_clawpack46_define_auxarray2(fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    fc2d_clawpack46_define_auxarray(domain,cp);
}

void fc2d_clawpack46_aux_data(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              double **aux, int* maux)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    fc2d_clawpack46_get_auxarray(domain,cp,aux,maux);
}

void fc2d_clawpack46_maux(fclaw2d_domain_t* domain,int* maux)
{
    *maux = fc2d_clawpack46_get_maux(domain);
}

int fc2d_clawpack46_get_maux(fclaw2d_domain_t* domain)
{
    fc2d_clawpack46_options_t *clawpack_options;
    clawpack_options = fc2d_clawpack46_get_options(domain);
    return clawpack_options->maux;
}

void fc2d_clawpack46_define_auxarray(fclaw2d_domain_t* domain, ClawPatch *cp)
{
    int maux;
    fc2d_clawpack46_options_t *clawpack_options;
    const amr_options_t* gparms;
    int mx,my,mbc;
    int ll[2], ur[2];

    gparms = get_domain_parms(domain);
    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    clawpack_options = fc2d_clawpack46_get_options(domain);
    maux = clawpack_options->maux;

    /* Construct an index box to hold the aux array */
    ll[0] = 1-mbc;
    ll[1] = 1-mbc;
    ur[0] = mx + mbc;
    ur[1] = my + mbc;    Box box(ll,ur);

    /* get solver specific data stored on this patch */
    patch_aux_data *clawpack_patch_data = get_patch_data(cp);
    clawpack_patch_data->auxarray.define(box,maux);
    clawpack_patch_data->maux = maux; // set maux in solver specific patch data (for completeness)
}

void fc2d_clawpack46_get_auxarray(fclaw2d_domain_t* domain,
                                  ClawPatch *cp, double **aux, int* maux)
{
    fc2d_clawpack46_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack46_get_options(domain);
    *maux = clawpack_options->maux;

    patch_aux_data *clawpack_patch_data = get_patch_data(cp);
    *aux = clawpack_patch_data->auxarray.dataPtr();
}

void fc2d_clawpack46_setprob(fclaw2d_domain_t *domain)
{
    /* Assume that if we are here, then the user has a valid setprob */
    FCLAW_ASSERT(classic_vt.setprob != NULL);

    classic_vt.setprob();
}


/* This should only be called when a new ClawPatch is created. */
void fc2d_clawpack46_setaux(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    FCLAW_ASSERT(classic_vt.setaux != NULL);
    int mx,my,mbc,maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_clawpack46_define_auxarray2(domain,this_patch);

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    maxmx = mx;
    maxmy = my;

    CLAWPACK46_SET_BLOCK(&this_block_idx);
    classic_vt.setaux(&maxmx,&maxmy,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                      &maux,aux);
    CLAWPACK46_UNSET_BLOCK();
}

void fc2d_clawpack46_qinit(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    FCLAW_ASSERT(classic_vt.qinit != NULL); /* Must initialized */
    int mx,my,mbc,meqn,maux,maxmx,maxmy;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    maxmx = mx;
    maxmy = my;

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    CLAWPACK46_SET_BLOCK(&this_block_idx);
    classic_vt.qinit(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                     &maux,aux);
    CLAWPACK46_UNSET_BLOCK();
}

void fc2d_clawpack46_b4step2(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t, double dt)

{
    FCLAW_ASSERT(classic_vt.b4step2 != NULL);

    int mx,my,mbc,meqn, maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    maxmx = mx;
    maxmy = my;

    CLAWPACK46_SET_BLOCK(&this_block_idx);
    classic_vt.b4step2(&maxmx,&maxmy,&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                       &dx,&dy,&t,&dt,&maux,aux);
    CLAWPACK46_UNSET_BLOCK();
}

void fc2d_clawpack46_src2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    FCLAW_ASSERT(classic_vt.src2 != NULL);

    int mx,my,mbc,meqn, maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    maxmx = mx;
    maxmy = my;

    CLAWPACK46_SET_BLOCK(&this_block_idx);
    classic_vt.src2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
                    &dx,&dy,q,&maux,aux,&t,&dt);
    CLAWPACK46_UNSET_BLOCK();
}


/* Use this to return only the right hand side of the clawpack algorithm */
double fc2d_clawpack46_step2_rhs(fclaw2d_domain_t *domain,
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


void fc2d_clawpack46_bc2(fclaw2d_domain *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         fclaw_bool intersects_phys_bdry[],
                         fclaw_bool time_interp)
{
    FCLAW_ASSERT(classic_vt.bc2 != NULL);

    int mx,my,mbc,meqn, maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    ClawPatch *cp = get_clawpatch(this_patch);

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    maxmx = mx;
    maxmy = my;

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

    /*
      We may be imposing boundary conditions on time-interpolated data;
      and is being done just to help with fine grid interpolation.
      In this case, this boundary condition won't be used to update
      anything
    */
    q = cp->q_time_sync(time_interp);

    CLAWPACK46_SET_BLOCK(&this_block_idx);
    classic_vt.bc2(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    CLAWPACK46_UNSET_BLOCK();
}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double fc2d_clawpack46_step2(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    const amr_options_t* gparms;
    fc2d_clawpack46_options_t* clawpack_options;
    ClawPatch *cp;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(classic_vt.rpn2 != NULL);
    FCLAW_ASSERT(classic_vt.rpt2 != NULL);

    gparms = fclaw2d_forestclaw_get_options(domain);
    clawpack_options = fc2d_clawpack46_get_options(domain);

    cp = get_clawpatch(this_patch);

    level = this_patch->level;

    fc2d_clawpack46_aux_data(domain,this_patch,&aux,&maux);

    cp->save_current_step();  // Save for time interpolation

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&qold,&meqn);

    int mwaves = clawpack_options->mwaves;

    int maxm = max(mx,my);

    double cflgrid = 0.0;

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    int ierror = 0;
    fc2d_clawpack46_flux2_t flux2 = clawpack_options->use_fwaves ?
                                    CLAWPACK46_FLUX2FW : CLAWPACK46_FLUX2;

    CLAWPACK46_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, clawpack_options->method,
                          clawpack_options->mthlim, &clawpack_options->mcapa,
                          &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                          work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                          classic_vt.rpn2, classic_vt.rpt2,flux2,
                          cp->block_corner_count(), &ierror);

    FCLAW_ASSERT(ierror == 0);

    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}

double fc2d_clawpack46_update(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double dt)
{
    const fc2d_clawpack46_options_t* clawpack_options;
    clawpack_options = fc2d_clawpack46_get_options(domain);

    if (classic_vt.b4step2 != NULL)
    {
        fc2d_clawpack46_b4step2(domain,
                                this_patch,
                                this_block_idx,
                                this_patch_idx,t,dt);
    }

    double maxcfl = fc2d_clawpack46_step2(domain,
                                          this_patch,
                                          this_block_idx,
                                          this_patch_idx,t,dt);

    if (clawpack_options->src_term > 0)
    {
        fc2d_clawpack46_src2(domain,
                             this_patch,
                             this_block_idx,
                             this_patch_idx,t,dt);
    }
    return maxcfl;
}
