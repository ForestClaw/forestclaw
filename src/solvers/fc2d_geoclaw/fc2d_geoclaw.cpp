/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Yu-Hsuan Shih, Scott Aiton
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

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_options.h"
#include "fc2d_geoclaw_fort.h"
#include "fc2d_geoclaw_output_ascii.h"

#include <fclaw_pointer_map.h>

#include <fclaw_gauges.h>
#include "fc2d_geoclaw_gauges_default.h"

#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

/* Basic forestclaw functions */
#include <fclaw_patch.h>
#include <fclaw_convenience.h>  /* Needed to get search function for gauges */
#include "fclaw_options.h"
#include <fclaw_global.h>
#include <fclaw_vtable.h>
#include <fclaw_diagnostics.h>
#include <fclaw2d_defs.h>

/* Some mapping functions */
#include <fclaw2d_map_brick.h>
#include <fclaw_map.h>
#include <fclaw_map_query.h>


/* Needed for debugging */
#include "types.h"

struct region_type region_type_for_debug;

/* ----------------------------- static function defs ------------------------------- */
static
void geoclaw_setaux(fclaw_global_t *glob,
                    fclaw_patch_t *patch,
                    int blockno,
                    int patchno);



/* --------------------------- Creating/deleting patches ---------------------------- */

static
void geoclaw_patch_setup(fclaw_global_t *glob,
                         fclaw_patch_t *patch,
                         int blockno,
                         int patchno)
{
    geoclaw_setaux(glob,patch,blockno,patchno);
}


/* -------------------------- Clawpack solver functions ------------------------------ */


static
void geoclaw_setprob(fclaw_global_t *glob)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);
    if (geoclaw_vt->setprob != NULL)
        geoclaw_vt->setprob();
}

static
void geoclaw_qinit(fclaw_global_t *glob,
                   fclaw_patch_t *patch,
                   int blockno,
                   int patchno)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);
    FCLAW_ASSERT(geoclaw_vt->qinit != NULL); /* Must initialized */

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double* q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    FC2D_GEOCLAW_SET_BLOCK(&blockno);
    geoclaw_vt->qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                     &maux,aux);
    FC2D_GEOCLAW_UNSET_BLOCK();
}

static
void geoclaw_bc2(fclaw_global_t *glob,
                 fclaw_patch_t *patch,
                 int blockno,
                 int patchno,
                 double t,
                 double dt,
                 int intersects_phys_bdry[],
                 int time_interp)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);

    fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);
    FCLAW_ASSERT(geoclaw_vt->bc2 != NULL);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    int *block_mthbc = geo_opt->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[FCLAW2D_NUMFACES];
    for(int i = 0; i < FCLAW2D_NUMFACES; i++)
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
    int meqn;
    double *q;
    fclaw_clawpatch_timesync_data(glob,patch,time_interp,&q,&meqn);

    FC2D_GEOCLAW_SET_BLOCK(&blockno);
    geoclaw_vt->bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    FC2D_GEOCLAW_UNSET_BLOCK();

}

static
void geoclaw_setaux(fclaw_global_t *glob,
                    fclaw_patch_t *patch,
                    int blockno,
                    int patchno)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);
    FCLAW_ASSERT(geoclaw_vt->setaux != NULL);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    /* If this is a ghost patch, we only set aux values in ghost cells */
    int is_ghost = fclaw_patch_is_ghost(patch);
    int mint = 2*mbc;
    int nghost = mbc;


    FC2D_GEOCLAW_SET_BLOCK(&blockno);
    geoclaw_vt->setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                      &maux,aux,&is_ghost,&nghost,&mint);
    FC2D_GEOCLAW_UNSET_BLOCK();
}


static
void geoclaw_b4step2(fclaw_global_t *glob,
                     fclaw_patch_t *patch,
                     int blockno,
                     int patchno,
                     double t, double dt)

{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);

    if (geoclaw_vt->b4step2 != NULL)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;
        fclaw_clawpatch_2d_grid_data(glob,patch, &mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        int meqn;
        double *q;
        fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

        int maux;
        double *aux;
        fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

        FC2D_GEOCLAW_SET_BLOCK(&blockno);
        geoclaw_vt->b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                            &dx,&dy,&t,&dt,&maux,aux);
        FC2D_GEOCLAW_UNSET_BLOCK();
    }
}

static
void geoclaw_src2(fclaw_global_t *glob,
                  fclaw_patch_t *patch,
                  int blockno,
                  int patchno,
                  double t,
                  double dt)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);
    FCLAW_ASSERT(geoclaw_vt->src2 != NULL);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    FC2D_GEOCLAW_SET_BLOCK(&blockno);
    geoclaw_vt->src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                    &dx,&dy,q,&maux,aux,&t,&dt);
    FC2D_GEOCLAW_UNSET_BLOCK();

}

/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double geoclaw_step2(fclaw_global_t *glob,
                     fclaw_patch_t *patch,
                     int blockno,
                     int patchno,
                     double t,
                     double dt)
{
    fc2d_geoclaw_vtable_t *geoclaw_vt = fc2d_geoclaw_vt(glob);
    FCLAW_ASSERT(geoclaw_vt->rpn2 != NULL);
    FCLAW_ASSERT(geoclaw_vt->rpt2 != NULL);

    fc2d_geoclaw_options_t* geoclaw_options = fc2d_geoclaw_get_options(glob);

    int level = patch->level;

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    /* In case we need to re-take this step */
    fclaw_clawpatch_save_current_step(glob, patch);

    int mx, my, mbc; 
    double xlower, ylower, dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *qold;
    fclaw_clawpatch_soln_data(glob,patch,&qold,&meqn);

    /* Allocate work arrays */
    int mwaves = geoclaw_options->mwaves;
    int maxm = fmax(mx,my);
    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    int* block_corner_count = fclaw_patch_block_corner_count(glob,patch);

    /* Still need conservation fix */


    double cflgrid = 0.0;
    FC2D_GEOCLAW_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, geoclaw_options->method,
                       geoclaw_options->mthlim, &geoclaw_options->mcapa,
                       &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                       work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                       geoclaw_vt->rpn2, geoclaw_vt->rpt2,
                       block_corner_count);

    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}


#if 0
void fc2d_geoclaw_dt(fclaw_global_t *glob,double t, double* dt)
{
    /* What was the idea here? */
}
#endif



static
double geoclaw_update(fclaw_global_t *glob,
                      fclaw_patch_t *patch,
                      int blockno,
                      int patchno,
                      double t,
                      double dt,
                      void* user)
{
    FC2D_GEOCLAW_TOPO_UPDATE(&t);

    geoclaw_b4step2(glob,
                    patch,
                    blockno,
                    patchno,t,dt);

    double maxcfl = geoclaw_step2(glob,
                                  patch,
                                  blockno,
                                  patchno,t,dt);

    const fc2d_geoclaw_options_t* geoclaw_opt = fc2d_geoclaw_get_options(glob);
    if (geoclaw_opt->src_term > 0)
    {
        geoclaw_src2(glob,
                     patch,
                     blockno,
                     patchno,t,dt);
    }

    return maxcfl;
}


/* --------------------------------- Output functions ---------------------------- */

static
void geoclaw_output(fclaw_global_t *glob, int iframe)
{
    const fc2d_geoclaw_options_t*geo_opt = fc2d_geoclaw_get_options(glob);
    if (geo_opt->ascii_out != 0)
        fc2d_geoclaw_output_ascii(glob,iframe);        
}


/* ------------------------------ Regridding functions --------------------------- */


static
int geoclaw_patch_tag4refinement(fclaw_global_t *glob,
                                 fclaw_patch_t *patch,
                                 int blockno, 
                                 int patchno,
                                 int initflag)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    int level = patch->level;
    double t = glob->curr_time;

    int tag_patch;
    /* First check to see if we are forced to refine based on regions 
       If patch intersects a region (including time interval), this routine
       returns :  

          -- level >= maximum level allowed by any region 
             this patch intersects with. (tag_patch = 0)

          -- level < minimum level required by any region
             this patch intersects with. (tag_patch = 1)

        Otherwise, tag_patch = -1 and we should refine using usual criteria.
    */
    double xupper = xlower + mx*dx;
    double yupper = ylower + my*dy;
    int refine = 1;  /* We are tagging for refinement */
    FC2D_GEOCLAW_TEST_REGIONS(&level,&xlower,&ylower,&xupper,&yupper,
                              &t,&refine, &tag_patch);

    if (tag_patch < 0)
    {
        /* Need maxlevel to get length speed_tolerance - hackish? */
        const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);
        int maxlevel = fclaw_opt->maxlevel;
        FC2D_GEOCLAW_FORT_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,
                                         &dx,&dy,&t,&blockno,q,aux,&level,&maxlevel,
                                         &initflag,&tag_patch);
    }

    return tag_patch;
}


static
int geoclaw_patch_tag4coarsening(fclaw_global_t *glob,
                                 fclaw_patch_t *fine_patches,
                                 int blockno, 
                                 int patchno,
                                 int initflag)

{
    int mx,my,mbc,meqn,maux;
    double xlower[4],ylower[4],dx,dy;
    double *q[4], *aux[4];
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
        fclaw_clawpatch_aux_data(glob,&fine_patches[igrid],&aux[igrid],&maux);

        fclaw_clawpatch_2d_grid_data(glob,&fine_patches[igrid],&mx,&my,&mbc,
                                    &xlower[igrid],&ylower[igrid],&dx,&dy);
    }

    int level = fine_patches[0].level;
    double t = glob->curr_time;
    int tag_patch;

    /* Test parent quadrant : If any of the four sibling patches are in the 
       region, we consider that an intersection.  Assume Morton ordering
       on the sibling patches (0=ll, 1=lr, 2=ul, 3=ur) */
    double xupper = xlower[1] + mx*dx;
    double yupper = ylower[2] + my*dy;
    int refine = 0;  /* We are tagging for refinement */
    FC2D_GEOCLAW_TEST_REGIONS(&level,&xlower[0],&ylower[0],&xupper,&yupper,
                              &t,&refine, &tag_patch);

    if (tag_patch < 0) 
    {
        /* Region tagging is inconclusive */
        const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);
        int maxlevel = fclaw_opt->maxlevel;

        FC2D_GEOCLAW_FORT_TAG4COARSENING(&blockno,&mx,&my,&mbc,&meqn,&maux,xlower,ylower,
                                         &dx,&dy, &t,q[0],q[1],q[2],q[3],
                                         aux[0],aux[1],aux[2],aux[3],
                                         &level,&maxlevel, &initflag, &tag_patch);
    }
    return tag_patch;

}

static
void geoclaw_interpolate2fine(fclaw_global_t *glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patches,
                              int blockno, 
                              int coarse_patchno,
                              int fine0_patchno)

{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw_patch_t* fine_patch = &fine_patches[igrid];

        double *qfine;
        fclaw_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        double *auxfine;
        fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        FC2D_GEOCLAW_FORT_INTERPOLATE2FINE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                           &maux,auxcoarse,auxfine, &igrid);
    }
}

static
void geoclaw_average2coarse(fclaw_global_t *glob,
                            fclaw_patch_t *fine_patches,
                            fclaw_patch_t *coarse_patch,
                            int blockno, 
                            int fine0_patchno,
                            int coarse_patchno)

{
    /* Only mx, my are needed here */
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw_patch_t *fine_patch = &fine_patches[igrid];

        double *qfine;
        fclaw_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        double *auxfine;
        fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        const fc2d_geoclaw_options_t* geo_opt = fc2d_geoclaw_get_options(glob);
        int mcapa = geo_opt->mcapa;
        FC2D_GEOCLAW_FORT_AVERAGE2COARSE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                         &maux,auxcoarse,auxfine,&mcapa,&igrid);
    }
}


/* ------------------------- Ghost filling - patch specific ------------------------ */

void geoclaw_average_face(fclaw_global_t *glob,
                          fclaw_patch_t *coarse_patch,
                          fclaw_patch_t *fine_patch,
                          int idir,
                          int iface_coarse,
                          int p4est_refineFactor,
                          int refratio,
                          int time_interp,
                          int igrid,
                          fclaw_patch_transform_data_t* transform_data)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double* qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    /* These will be empty for non-manifolds cases */
    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double *auxfine;
    fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    const fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);
    int mcapa = geo_opt->mcapa;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (manifold != 0)
    {
        fclaw_global_essentialf("fc2d_geoclaw : Manifold case is not handled explicit" \
                                  "in Geoclaw.");
        exit(0);
    }

    FC2D_GEOCLAW_FORT_AVERAGE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                   &maux,auxcoarse,auxfine,&mcapa,
                                   &idir,&iface_coarse,
                                   &igrid,&transform_data);
}

void geoclaw_interpolate_face(fclaw_global_t *glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patch,
                              int idir,
                              int iside,
                              int p4est_refineFactor,
                              int refratio,
                              int time_interp,
                              int igrid,
                              fclaw_patch_transform_data_t* transform_data)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double *auxfine;
    fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,&maux,
                                       auxcoarse,auxfine, &idir, &iside,
                                       &igrid, &transform_data);
}

void geoclaw_average_corner(fclaw_global_t *glob,
                            fclaw_patch_t *coarse_patch,
                            fclaw_patch_t *fine_patch,
                            int coarse_blockno,
                            int fine_blockno,
                            int coarse_corner,
                            int time_interp,
                            fclaw_patch_transform_data_t* transform_data)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double *auxfine;
    fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    const fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);
    int mcapa = geo_opt->mcapa;

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (manifold != 0)
    {
        fclaw_global_essentialf("fc2d_geoclaw : Manifold case is not handled explicit" \
                                  "in Geoclaw.");
        exit(0);
    }

    FC2D_GEOCLAW_FORT_AVERAGE_CORNER(&mx,&my,&mbc,&meqn,
                                     qcoarse,qfine,&maux,auxcoarse,auxfine,
                                     &mcapa,&coarse_corner,
                                     &transform_data);

}

void geoclaw_interpolate_corner(fclaw_global_t* glob,
                                fclaw_patch_t* coarse_patch,
                                fclaw_patch_t* fine_patch,
                                int coarse_blockno,
                                int fine_blockno,
                                int coarse_corner,
                                int time_interp,
                                fclaw_patch_transform_data_t* transform_data)

{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double* qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int maux;
    double *auxcoarse;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double *auxfine;
    fclaw_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER(&mx,&my,&mbc,&meqn,
                                         qcoarse,qfine,&maux,
                                         auxcoarse,auxfine,
                                         &coarse_corner,&transform_data);

}

/* --------------------------- Parallel ghost patches -------------------------------- */

void geoclaw_remote_ghost_setup(fclaw_global_t *glob,
                                fclaw_patch_t *patch,
                                int blockno,
                                int patchno)
{
    fclaw_clawpatch_options_t* clawpatch_options;
    clawpatch_options = fclaw_clawpatch_get_options(glob);

    if (!clawpatch_options->ghost_patch_pack_aux)
    {
        geoclaw_setaux(glob,patch,blockno,patchno);
    }
    else
    {
        /* the aux array data has been packed and transferred as MPI messages */
    }
}

static
void geoclaw_local_ghost_pack_aux(fclaw_global_t *glob,
                                  fclaw_patch_t *patch,
                                  int mint,
                                  double *auxpack,
                                  int auxsize, int packmode,
                                  int* ierror)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        int maux;
        double *aux;
        fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);
        FC2D_GEOCLAW_LOCAL_GHOST_PACK_AUX(&mx,&my,&mbc,&maux,
                                          &mint,aux,auxpack,&auxsize,
                                          &packmode,ierror);
}

/* ------------------------------ Misc access functions ----------------------------- */

/* Called from application routines */
void fc2d_geoclaw_module_setup(fclaw_global_t *glob)
{
    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    const fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);

    FC2D_GEOCLAW_SET_MODULES(&geo_opt->mwaves, 
                             &geo_opt->mcapa,
                             &clawpatch_opt->meqn, 
                             &clawpatch_opt->maux,
                             geo_opt->mthlim, 
                             geo_opt->method,
                             &fclaw_opt->ax, 
                             &fclaw_opt->bx, 
                             &fclaw_opt->ay, 
                             &fclaw_opt->by);
}

/* -------------------------- Virtual table  ---------------------------- */

static
fc2d_geoclaw_vtable_t* fc2d_geoclaw_vt_new()
{
    return (fc2d_geoclaw_vtable_t*) FCLAW_ALLOC_ZERO (fc2d_geoclaw_vtable_t, 1);
}

static
void fc2d_geoclaw_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fc2d_geoclaw_vtable_t* fc2d_geoclaw_vt(fclaw_global_t* glob)
{
	fc2d_geoclaw_vtable_t* geoclaw_vt = (fc2d_geoclaw_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc2d_geoclaw");
	FCLAW_ASSERT(geoclaw_vt != NULL);
	FCLAW_ASSERT(geoclaw_vt->is_set != 0);
	return geoclaw_vt;
}

void fc2d_geoclaw_solver_initialize(fclaw_global_t* glob)
{
	if(glob->domain->refine_dim != 2)
    {
		fclaw_abortf("Domain set to 3d. fc2d_geoclaw is only for 2d");
    }

	fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
	fc2d_geoclaw_options_t* geo_opt = fc2d_geoclaw_get_options(glob);

	if(clawpatch_opt->patch_dim != 2)
    {
		fclaw_abortf("Clawpatch dimension set to 3d. fc2d_geoclaw is only for 2d");
    }

    geo_opt->method[6] = clawpatch_opt->maux;

    /* We have to do this so that we know how to size the ghost patches */
    if (clawpatch_opt->ghost_patch_pack_aux)
    {
        fclaw_opt->ghost_patch_pack_extra = 1;  /* Pack the bathymetry */
        fclaw_opt->ghost_patch_pack_numextrafields = clawpatch_opt->maux;
    }

    int claw_version = 5;
    fclaw_clawpatch_vtable_initialize(glob, claw_version);
    
    fclaw_gauges_vtable_t*           gauges_vt = fclaw_gauges_vt(glob);

    fclaw_vtable_t*                fc_vt = fclaw_vt(glob);
    fclaw_patch_vtable_t*          patch_vt = fclaw_patch_vt(glob);
    fclaw_clawpatch_vtable_t*  clawpatch_vt = fclaw_clawpatch_vt(glob);

    fc2d_geoclaw_vtable_t*  geoclaw_vt = fc2d_geoclaw_vt_new();

    /* ForestClaw virtual tables */
    fc_vt->problem_setup               = geoclaw_setprob;  
    // fclaw_vt->after_regrid                = geoclaw_after_regrid;  /* Handle gauges */

    /* Set basic patch operations */
    patch_vt->setup                       = geoclaw_patch_setup;
    patch_vt->initialize                  = geoclaw_qinit;
    patch_vt->physical_bc                 = geoclaw_bc2;
    patch_vt->single_step_update          = geoclaw_update;  /* Includes b4step2 and src2 */
         
    fc_vt->output_frame                = geoclaw_output;

    /* Regridding */
    patch_vt->tag4refinement              = geoclaw_patch_tag4refinement;
    patch_vt->tag4coarsening              = geoclaw_patch_tag4coarsening;
    patch_vt->interpolate2fine            = geoclaw_interpolate2fine;
    patch_vt->average2coarse              = geoclaw_average2coarse;

    /* Ghost filling */
    clawpatch_vt->d2->fort_copy_face          = FC2D_GEOCLAW_FORT_COPY_FACE;
    clawpatch_vt->d2->fort_copy_corner        = FC2D_GEOCLAW_FORT_COPY_CORNER;

    /* Geoclaw needs specialized averaging and interpolation routines */
    patch_vt->average_face                    = geoclaw_average_face;
    patch_vt->interpolate_face                = geoclaw_interpolate_face;      
    patch_vt->average_corner                  = geoclaw_average_corner;
    patch_vt->interpolate_corner              = geoclaw_interpolate_corner;

    patch_vt->remote_ghost_setup              = geoclaw_remote_ghost_setup;
    clawpatch_vt->d2->fort_local_ghost_pack   = FC2D_GEOCLAW_LOCAL_GHOST_PACK;
    clawpatch_vt->local_ghost_pack_aux        = geoclaw_local_ghost_pack_aux;
  
    /* Diagnostic functions partially implemented in clawpatch */
    clawpatch_vt->d2->fort_compute_error_norm = FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->d2->fort_compute_patch_area = FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->d2->fort_conservation_check = FC2D_GEOCLAW_FORT_CONSERVATION_CHECK;
    clawpatch_vt->d2->fort_timeinterp         = FC2D_GEOCLAW_FORT_TIMEINTERP;

    geoclaw_vt->setprob          = NULL;                   
    geoclaw_vt->setaux           = FC2D_GEOCLAW_SETAUX;
    geoclaw_vt->qinit            = FC2D_GEOCLAW_QINIT;
    geoclaw_vt->bc2              = FC2D_GEOCLAW_BC2;
    geoclaw_vt->b4step2          = FC2D_GEOCLAW_B4STEP2;
    geoclaw_vt->src2             = FC2D_GEOCLAW_SRC2;
    geoclaw_vt->rpn2             = FC2D_GEOCLAW_RPN2;
    geoclaw_vt->rpt2             = FC2D_GEOCLAW_RPT2;

    gauges_vt->set_gauge_data     = geoclaw_read_gauges_data_default;
    gauges_vt->create_gauge_files = geoclaw_create_gauge_files_default; 
    gauges_vt->normalize_coordinates = geoclaw_gauge_normalize_coordinates;

    gauges_vt->update_gauge       = geoclaw_gauge_update_default;
    gauges_vt->print_gauge_buffer = geoclaw_print_gauges_default;

    geoclaw_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc2d_geoclaw") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc2d_geoclaw", geoclaw_vt, fc2d_geoclaw_vt_destroy);
}

