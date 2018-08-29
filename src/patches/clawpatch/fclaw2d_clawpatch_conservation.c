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
#include "fclaw2d_clawpatch_conservation_fort.h"

#include <fclaw2d_time_sync.h>

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

		cu->edge_fluxes[k]   = FCLAW_ALLOC_ZERO(double,2*my*meqn);
		cu->edge_fluxes[k+2] = FCLAW_ALLOC_ZERO(double,2*mx*meqn);

		cu->edgelengths[k]   = FCLAW_ALLOC(double,my);
		cu->edgelengths[k+2] = FCLAW_ALLOC(double,mx);
		cu->area[k]          = FCLAW_ALLOC(double,my);
		cu->area[k+2]        = FCLAW_ALLOC(double,mx);		
	}

}


void fclaw2d_clawpatch_time_sync_reset(fclaw2d_global_t *glob,
                                       fclaw2d_patch_t *this_patch, 
                                       int coarse_level,
                                       int reset_mode)
{
	int mx,my,meqn;
	int i,j,k,idir;

	fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

	fclaw2d_clawpatch_cons_update_t* cu = fclaw2d_clawpatch_get_cons_update(glob,
																			this_patch);
	mx = clawpatch_opt->mx;
	my = clawpatch_opt->my;
	meqn = clawpatch_opt->meqn;

	fclaw2d_patch_data_t* pdata = fclaw2d_patch_get_patch_data(this_patch);

	int fine_level = coarse_level+1;
	int reset_flux;

	for(k = 0; k < 4; k++)
	{
		idir = k/2;
		reset_flux = 0;
		if (reset_mode == FCLAW2D_TIME_SYNC_RESET_F2C)
		{
			int is_coarse = (pdata->face_neighbors[k] == FCLAW2D_PATCH_HALFSIZE)
          			&& (this_patch->level == coarse_level);
		    int is_fine = (pdata->face_neighbors[k] == FCLAW2D_PATCH_DOUBLESIZE) &&
		          (this_patch->level == fine_level);
		    reset_flux = is_coarse || is_fine;
		}
		else if (reset_mode == FCLAW2D_TIME_SYNC_RESET_SAMESIZE)
		{
			reset_flux = (pdata->face_neighbors[k] == FCLAW2D_PATCH_SAMESIZE) ||   
			      (this_patch->level == coarse_level);
		}
		else if (reset_mode == FCLAW2D_TIME_SYNC_RESET_LEVEL)
		{
			reset_flux = this_patch->level == coarse_level;
		}

		if (reset_flux)
		{
			if (idir == 0)
			{
				for(j = 0; j < meqn*my; j++)
				{
					cu->fm[k][j] = 0;
					cu->fp[k][j] = 0;
					cu->edge_fluxes[k][j] = 0;
					cu->edge_fluxes[k][j+meqn*my] = 0;  /* Two fluxes stored at each edge point */
				}
			}
			else
			{
				for(i = 0; i < meqn*mx; i++)
				{
					cu->gm[k-2][i] = 0;
					cu->gp[k-2][i] = 0;
					cu->edge_fluxes[k][i] = 0;
					cu->edge_fluxes[k][i + meqn*mx] = 0;
				}
			}
		}     /* coarse grid neighbor */
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

		FCLAW_FREE(cu->edgelengths[k]);
		FCLAW_FREE(cu->area[k]);
		FCLAW_FREE(cu->edgelengths[k+2]);
		FCLAW_FREE(cu->area[k+2]);
	 }
	 FCLAW_FREE(*cons_update);
	 *cons_update = NULL;
}


void fclaw2d_clawpatch_update_cons_metric(fclaw2d_global_t* glob,
										  fclaw2d_patch_t* this_patch,
										  int blockno,int patchno)
{
	int mx,my,mbc;
	double dx,dy,xlower,ylower;
	double *area, *edgelengths, *curvature;

	fclaw2d_clawpatch_cons_update_t *cu = 
	                       fclaw2d_clawpatch_get_cons_update(glob,this_patch);

	const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_metric_scalar(glob,this_patch,
									&area, &edgelengths,&curvature);

	CLAWPATCH_UPDATE_CONS_METRIC(&mx,&my,&mbc,&dx,&dy,area,edgelengths,
	                             cu->area[0],cu->area[1],
	                             cu->area[2],cu->area[3],
	                             cu->edgelengths[0],cu->edgelengths[1],
	                             cu->edgelengths[2],cu->edgelengths[3],
	                             &fclaw_opt->manifold);
}


/* This is a patch call-back */
void fclaw2d_clawpatch_time_sync_f2c(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* coarse_patch,
                                     fclaw2d_patch_t* fine_patch,
                                     int idir,
                                     int igrid,
                                     int iface_coarse,
                                     int time_interp,
                                     fclaw2d_patch_transform_data_t
                                     *transform_data)
{

	int meqn,mx,my,mbc;
	double dx,dy,xlower,ylower;
	double *qcoarse, *qfine;
	int coarse_blockno, coarse_patchno,globnum,level;
	int fine_blockno, fine_patchno;

	fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

	if (time_interp)
	{
		return;
	}

	/* We don't correct time interpolated grids, since we assume that the time 
	interpolated average will already have correction information from times 
	n and n+1 patches.  But we do correct ghost patches, since corrections will be 
	needed for copying and interpolation.  */

	fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);
	fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

	fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_cons_update_t* cucoarse = 
	           fclaw2d_clawpatch_get_cons_update(glob,coarse_patch);

	fclaw2d_clawpatch_cons_update_t* cufine = 
	           fclaw2d_clawpatch_get_cons_update(glob,fine_patch);

	/* create dummy fine grid to handle indexing between blocks */
	double *qneighbor_dummy = FCLAW_ALLOC_ZERO(double,meqn*(mx+2*mbc)*(my+2*mbc));
	int *maskneighbor = FCLAW_ALLOC_ZERO(int,(mx+2*mbc)*(my+2*mbc));

	/* Include this for debugging */
#if 1	
	fclaw2d_patch_get_info2(glob->domain,coarse_patch,&coarse_blockno, &coarse_patchno,
							&globnum,&level);

	fclaw2d_patch_get_info2(glob->domain,fine_patch,&fine_blockno, 
							&fine_patchno,&globnum,&level);
#endif							

	clawpatch_vt->fort_time_sync_fine_to_coarse(&mx,&my,&mbc,&meqn,&idir,&iface_coarse,
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


void fclaw2d_clawpatch_time_sync_samesize (struct fclaw2d_global* glob,
                                           struct fclaw2d_patch* this_patch,
                                           struct fclaw2d_patch* neighbor_patch,
                                           int this_iface,int idir,
                                           fclaw2d_patch_transform_data_t
                                           *transform_data)
{
	/* We don't correct time interpolated grids, since we assume that the time 
	interpolated average will already have correction information from times 
	n and n+1 patches.  But we do correct ghost patches, since corrections will be 
	needed for copying and interpolation.  */
	int mx,my,meqn,mbc;
	double dx,dy,xlower,ylower;
	double *qthis;

	fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();

	fclaw2d_clawpatch_soln_data(glob,this_patch,&qthis,&meqn);

	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

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
	fclaw2d_patch_get_info2(glob->domain,this_patch,&this_blockno, &this_patchno,
							&globnum,&level);

	fclaw2d_patch_get_info2(glob->domain,neighbor_patch,&neighbor_blockno, 
							&neighbor_patchno,&globnum,&level);

    if (1)
    {
    	/* Distribute 0.5 each correction from each side */
    	clawpatch_vt->fort_time_sync_samesize(&mx,&my,&mbc,&meqn,&idir,&this_iface,
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
    }
	FCLAW_FREE(qneighbor_dummy);
	FCLAW_FREE(maskneighbor);       

}





