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


#include "fclaw2d_clawpatch_conservation.h"

#include "fclaw2d_clawpatch.h"
#include "fclaw2d_clawpatch_options.h"

#include <fclaw2d_patch.h>
#include <fclaw2d_options.h>

#include <fclaw2d_global.h>
#include <fclaw_math.h>

/* -------------------------------------------------------------
    Four routines here : 
    1. fclaw2d_clawpatch_cons_update_new
    2. fclaw2d_clawpatch_cons_update_reset
    3. fclaw2d_clawpatch_cons_update_delete
    4. fclaw2d_clawpatch_time_sync_fine_to_coarse
    5. fclaw2d_clawpatch_time_sync_copy
*/


static
void set_snan(double* f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)&f) = 0x7ff0000000000001LL;
    // *f = -9999.0;
}


void fclaw2d_clawpatch_cons_update_new (fclaw2d_global_t* glob,
                                       fclaw2d_patch_t* this_patch,
                                       int blockno,int patchno,
                                       fclaw2d_clawpatch_cons_update_t **cons_update)
{
    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int k;

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int meqn = clawpatch_opt->meqn;

    fclaw2d_clawpatch_cons_update_t *cu = *cons_update;

    cu = FCLAW_ALLOC(fclaw2d_clawpatch_cons_update_t,1);
    *cons_update = cu;

    /* Iterate over sides 0,1,3,4 */
    for(k = 0; k < 2; k++)
    {        
        /* Accumulators */
        cu->fp[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);
        cu->fm[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);

        cu->gp[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);
        cu->gm[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);

        cu->edgelengths[k]   = FCLAW_ALLOC(double,my);
        cu->edgelengths[k+2] = FCLAW_ALLOC(double,mx);
        cu->area[k]          = FCLAW_ALLOC(double,my);
        cu->area[k+2]        = FCLAW_ALLOC(double,mx);

        cu->edge_fluxes[k]   = FCLAW_ALLOC(double,2*my*meqn);
        cu->edge_fluxes[k+2] = FCLAW_ALLOC(double,2*mx*meqn);
#if 0         
        cu->edge_fluxes_save[k]   = FCLAW_ALLOC(double,2*my*meqn);
        cu->edge_fluxes_save[k+2] = FCLAW_ALLOC(double,2*mx*meqn);
#endif        
    }

#if 0
    /* Set qc, auxc to signaling NANS, to make sure we only use valid values */
    for(k = 0; k < 2; k++)
    {
        for(j = 0; j < meqn*my; j++)
        {
            set_snan(&cu->qc[k][j]);
            set_snan(&cu->auxc[k][j]);
            set_snan(&cu->qc[k+2][j]);
            set_snan(&cu->auxc[k+2][j]);
        }

        for(i = 0; i < meqn*mx; i++)
        {
            set_snan(&cu->qc[k][i]);
            set_snan(&cu->auxc[k][i]);
            set_snan(&cu->qc[k+2][i]);
            set_snan(&cu->auxc[k+2][i]);
        }
    }
#endif     
}


static
void cb_cons_update_reset (fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    int mx,my,meqn,level,minlevel;
    int i,j,k,idir;

    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;

    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(g->glob);

    fclaw2d_clawpatch_cons_update_t* cu = fclaw2d_clawpatch_get_cons_update(g->glob,
                                                                            this_patch);

    minlevel = *((int*) g->user);
    level = this_patch->level;

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    meqn = clawpatch_opt->meqn;

    fclaw2d_patch_data_t* pdata = fclaw2d_patch_get_user_data(this_patch);

    for(k = 0; k < 4; k++)
    {
        idir = k/2;
        if (((pdata->face_neighbors[k] == FCLAW2D_PATCH_DOUBLESIZE) && level > minlevel)
            || ((pdata->face_neighbors[k] == FCLAW2D_PATCH_HALFSIZE)))
        {
            if (idir == 0)
            {
                for(j = 0; j < meqn*my; j++)
                {                 
                    cu->fm[k][j] = 0;
                    cu->fp[k][j] = 0;
                }
            }
            else
            {
                for(i = 0; i < meqn*mx; i++)
                {
                    cu->gm[k-2][i] = 0;
                    cu->gp[k-2][i] = 0;
                }
            }
        }     /* coarse grid neighbor */
    }
}

void fclaw2d_clawpatch_cons_update_reset (fclaw2d_global_t* glob,int minlevel,
                                         int maxlevel)
{
    int level;

    for(level = minlevel; level <= maxlevel; level++)
    {
        fclaw2d_global_iterate_level(glob, level, cb_cons_update_reset, &minlevel);
    }
}



void fclaw2d_clawpatch_cons_update_delete (fclaw2d_clawpatch_cons_update_t **cons_update)
{
    int k;

    fclaw2d_clawpatch_cons_update_t *cu = *cons_update;

    for(k = 0; k < 2; k++)
    {
        /* Accumulators */
        FCLAW_FREE(cu->fp[k]);
        FCLAW_FREE(cu->fm[k]);    

        FCLAW_FREE(cu->gp[k]);    
        FCLAW_FREE(cu->gm[k]);  

        /* COARSE GRID information */
        FCLAW_FREE(cu->edge_fluxes[k]); 
        FCLAW_FREE(cu->edge_fluxes[k+2]); 
#if 0        
        FCLAW_FREE(cu->edge_fluxes_save[k]); 
        FCLAW_FREE(cu->edge_fluxes_save[k+2]); 
#endif        

        FCLAW_FREE(cu->edgelengths[k]);
        FCLAW_FREE(cu->area[k]);
        FCLAW_FREE(cu->edgelengths[k+2]);
        FCLAW_FREE(cu->area[k+2]);
     }
     FCLAW_FREE(*cons_update);
     *cons_update = NULL;
}



/* This is a patch call-back */
void fclaw2d_clawpatch_time_sync_fine_to_coarse (fclaw2d_global_t* glob,
                                                 fclaw2d_patch_t* coarse_patch,
                                                 fclaw2d_patch_t* fine_patch,
                                                 int idir,
                                                 int igrid,
                                                 int iface_coarse,
                                                 int time_interp,
                                                 fclaw2d_transform_data_t* transform_data)
{

    int meqn,mx,my,mbc;
    double *qcoarse, *qfine;
    int coarse_blockno, coarse_patchno,globnum,level;
    int fine_blockno, fine_patchno;


    if (time_interp != NULL)
    {
        return;
    }

    /* We don't correct time interpolated grids, since we assume that the time 
    interpolated average will already have correction information from times 
    n and n+1 patches.  But we do correct ghost patches, since corrections will be 
    needed for copying and interpolation.  */

    qcoarse = fclaw2d_clawpatch_get_q(glob,coarse_patch);
    qfine   = fclaw2d_clawpatch_get_q(glob,fine_patch);

    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    meqn = clawpatch_opt->meqn;


    fclaw2d_clawpatch_cons_update_t* cucoarse = 
    fclaw2d_clawpatch_get_cons_update(glob,coarse_patch);

    fclaw2d_clawpatch_cons_update_t* cufine = 
    fclaw2d_clawpatch_get_cons_update(glob,fine_patch);

    /* create dummy fine grid to handle indexing between blocks */
    double *qneighbor_dummy = FCLAW_ALLOC_ZERO(double,meqn*(mx+2*mbc)*(my+2*mbc));
    int *maskneighbor = FCLAW_ALLOC_ZERO(int,(mx+2*mbc)*(my+2*mbc));

    /* Include this for debugging */
    fclaw2d_patch_get_info2(glob->domain,coarse_patch,&coarse_blockno, &coarse_patchno,
                            &globnum,&level);

    fclaw2d_patch_get_info2(glob->domain,fine_patch,&fine_blockno, 
                            &fine_patchno,&globnum,&level);

    clawpatch_vt->fort_time_sync_fine_to_coarse(&mx,&my,&mbc,&meqn,&idir,&iface,
                                                cucoarse->area[0], cucoarse->area[1], 
                                                cucoarse->area[2], cucoarse->area[3],
                                                qcoarse,
                                                cucoarse->fp[0],cucoarse->fm[1],
                                                cucoarse->gp[0],cucoarse->gm[1],
                                                cufine->fm[0],cufine->fp[1],
                                                cufine->gm[0],cufine->gp[1],
                                                cucoarse->edge_fluxes[0],
                                                cucoarse->edge_fluxes[1],
                                                cucoarse->edge_fluxes[2],
                                                cucoarse->edge_fluxes[3],
                                                cufine->edge_fluxes[0],
                                                cufine->edge_fluxes[1],
                                                cufine->edge_fluxes[2],
                                                cufine->edge_fluxes[3],
                                                maskneighbor,qneighbor_dummy,
                                                &transform_data);

    FCLAW_FREE(qneighbor_dummy);
    FCLAW_FREE(maskneighbor);       
}


void fclaw2d_clawpatch_time_sync_copy (struct fclaw2d_global* glob,
                                       struct fclaw2d_patch* this_patch,
                                       struct fclaw2d_patch* neighbor_patch,
                                       int this_iface,
                                       fclaw2d_transform_data_t* transform_data)
{
    /* We don't correct time interpolated grids, since we assume that the time 
    interpolated average will already have correction information from times 
    n and n+1 patches.  But we do correct ghost patches, since corrections will be 
    needed for copying and interpolation.  */

    fclaw2d_clawpatch_cons_update_t* cuthis = 
    fclaw2d_clawpatch_get_cons_update(glob,this_patch);

    fclaw2d_clawpatch_cons_update_t* cuneighbor = 
    fclaw2d_clawpatch_get_cons_update(glob,neighbor_patch);

    /* create dummy fine grid to handle indexing between blocks */
    double *qneighbor_dummy = FCLAW_ALLOC_ZERO(double,meqn*(mx+2*mbc)*(my+2*mbc));
    int *maskneighbor = FCLAW_ALLOC_ZERO(int,(mx+2*mbc)*(my+2*mbc));

    /* Include this for debugging */
    int this_blockno, this_patchno,globnum,level;
    int neighbor_blockno, neighbor_patchno;
    fclaw2d_patch_get_info2(glob->domain,this_patch,&neighbor_blockno, &this_patchno,
                            &globnum,&level);

    fclaw2d_patch_get_info2(glob->domain,neighbor_patch,&neighbor_blockno, 
                            &neighbor_patchno,&globnum,&level);

    clawpatch_vt->fort_time_sync_fine_to_coarse(&mx,&my,&mbc,&meqn,&idir,&iface,
                                                cuthis->area[0], cuthis->area[1], 
                                                cuthis->area[2], cuthis->area[3],
                                                qthis,
                                                cuthis->fp[0],cuthis->fm[1],
                                                cuthis->gp[0],cuthis->gm[1],
                                                cuneighbor->fm[0],cuneighbor->fp[1],
                                                cuneighbor->gm[0],cuneighbor->gp[1],
                                                cuthis->edge_fluxes[0],
                                                cuthis->edge_fluxes[1],
                                                cuthis->edge_fluxes[2],
                                                cuthis->edge_fluxes[3],
                                                cuneighbor->edge_fluxes[0],
                                                cuneighbor->edge_fluxes[1],
                                                cuneighbor->edge_fluxes[2],
                                                cuneighbor->edge_fluxes[3],
                                                maskneighbor,qneighbor_dummy,
                                                &transform_data);

    FCLAW_FREE(qneighbor_dummy);
    FCLAW_FREE(maskneighbor);       

}





