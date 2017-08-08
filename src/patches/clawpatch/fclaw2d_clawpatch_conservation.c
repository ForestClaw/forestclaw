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

#include "fclaw2d_patch.h"

#include <fclaw2d_global.h>
#include <fclaw_math.h>


typedef struct cons_update_info
{
    double dt;
    int minlevel;
    int maxlevel;
} cons_update_info_t;

static
void set_snan(double* f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)&f) = 0x7ff0000000000001LL;
}


void fclaw2d_clawpatch_cons_update_new(fclaw2d_global_t* glob,
                                       fclaw2d_patch_t* this_patch,
                                       int blockno,int patchno,
                                       fclaw2d_clawpatch_cons_update_t **cons_update)
{
    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    int i,j,k;

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int maux = clawpatch_opt->maux;
    int meqn = clawpatch_opt->meqn;

    fclaw2d_clawpatch_cons_update_t *cu = *cons_update;

    cu = FCLAW_ALLOC(fclaw2d_clawpatch_cons_update_t,1);
    *cons_update = cu;

    cu->dt = -1;      /* We don't have a time step yet */

    /* Iterate over sides 0,1,3,4 */
    for(k = 0; k < 2; k++)
    {        
        /* Accumulators */
        cu->fp[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);
        cu->fm[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);
        cu->rp[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);       /* apdq + amdq = f(qghost) - f(qcourse) */

        cu->gp[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);
        cu->gm[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);
        cu->rp[k+2]   = FCLAW_ALLOC_ZERO(double,mx*meqn);     

        /* Coarse grid information */
        cu->qc[k]     = FCLAW_ALLOC(double,my*meqn);      /* left, right side */
        cu->auxc[k]   = FCLAW_ALLOC(double,my*maux);      /* left, right side */

        cu->qc[k+2]   = FCLAW_ALLOC(double,mx*meqn);      /* bottom, top side */
        cu->auxc[k+2] = FCLAW_ALLOC(double,mx*maux);      /* bottom, top side */ 
    }
#ifdef FCLAW_ENABLE_DEBUG
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
void cb_cons_update_reset(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          void *user)
{
    fclaw2d_global_iterate_t* g = (fclaw2d_global_iterate_t*) user;
 
    cons_update_info_t cons_info = *((cons_update_info_t*) g->user);

    fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(g->glob);

    fclaw2d_clawpatch_cons_update_t* cu = fclaw2d_clawpatch_get_cons_update(g->glob,
                                                                            this_patch);

    int level = this_patch->level;
    int minlevel = cons_info.minlevel;
    double dt = cons_info.dt;

    cu->dt = dt/pow_int(2,level-minlevel);

    int i, j, k, idir;

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int meqn = clawpatch_opt->meqn;

    fclaw2d_patch_data_t* pdata = fclaw2d_patch_get_user_data(this_patch);

    for(k = 0; k < 4; k++)
    {
        if (pdata->face_neighbors[k] == FCLAW2D_PATCH_DOUBLESIZE)
        {
            idir = k % 2;
            if (idir == 0)
            {
                for(j = 0; j < meqn*my; j++)
                {                    
                    if (k == 0)
                    {
                        cu->fm[0][j] = 0;
                    }
                    else if (k == 1)
                    {
                        cu->fp[1][j] = 0;
                    }
                }
            }
            else
            {
                for(i = 0; i < meqn*mx; i++)
                {
                    if (k == 2)
                    {
                        cu->gm[0][i] = 0;
                    }
                    else if (k == 3)
                    {
                        cu->gp[1][i] = 0;
                    }
                }
            }
        }     /* coarse grid neighbor */
    }

#if 0
    for(k = 0; k < 2; k++)
    {
        for(j = 0; j < meqn*my; j++)
        {
            cu->fp[k][j] = 0;
            cu->fm[k][j] = 0;
            cu->rp[k][j] = 0;
        }

        for(i = 0; i < meqn*mx; i++)
        {
            cu->gp[k][i] = 0;
            cu->gm[k][i] = 0;
            cu->rp[k+2][i] = 0;  
        }
    }
#endif    
}

void fclaw2d_clawpatch_cons_update_reset(fclaw2d_global_t* glob,int minlevel,
                                         int maxlevel, double dt)
{
    int level;

    cons_update_info_t cons_info;

    cons_info.dt = dt;
    cons_info.minlevel = minlevel;
    cons_info.maxlevel = maxlevel;

    /* Reset fine grids only;  coarse grids do not accumulate fluxes; they are 
    set only once */
    for(level = minlevel+1; level <= maxlevel; level++)
    {
        fclaw2d_global_iterate_level(glob, level, cb_cons_update_reset, &cons_info);
    }
}



void fclaw2d_clawpatch_cons_update_delete(fclaw2d_clawpatch_cons_update_t **cons_update)
{
    int k;

    fclaw2d_clawpatch_cons_update_t *cu = *cons_update;

    for(k = 0; k < 2; k++)
    {
        /* Accumulators */
        FCLAW_FREE(cu->fp[k]);
        FCLAW_FREE(cu->fm[k]);    
        FCLAW_FREE(cu->rp[k]);    

        FCLAW_FREE(cu->gp[k]);    
        FCLAW_FREE(cu->gm[k]);    
        FCLAW_FREE(cu->rp[k+2]);

        /* COARSE GRID information */
        FCLAW_FREE(cu->qc[k]); 
        FCLAW_FREE(cu->qc[k+2]); 
        FCLAW_FREE(cu->auxc[k]); 
        FCLAW_FREE(cu->auxc[k+2]); 
     }
     FCLAW_FREE(*cons_update);
     *cons_update = NULL;
}



