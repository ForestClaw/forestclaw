/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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


#include <fclaw2d_time_sync.h>

#include <fclaw2d_patch.h>
#include <fclaw_options.h>

#include <fclaw_global.h>
#include <fclaw_math.h>

#include "fclaw2d_clawpatch_conservation.h"
#include "fclaw2d_clawpatch_conservation_fort.h"

#include "fclaw_clawpatch.h"
#include "fclaw_clawpatch_options.h"

/* -------------------------------------------------------------
	Four routines here : 
	1. fclaw2d_clawpatch_cons_update_new
	2. fclaw2d_clawpatch_cons_update_reset
	3. fclaw2d_clawpatch_cons_update_delete
	4. fclaw2d_clawpatch_time_sync_fine_to_coarse
	5. fclaw2d_clawpatch_time_sync_copy
*/



void fclaw2d_clawpatch_time_sync_new (fclaw_global_t* glob,
                                      fclaw_patch_t* this_patch,
                                      int blockno,int patchno,
                                      fclaw2d_clawpatch_registers_t **registers)
{
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);

	int k;

	int mx;
	int my;

	if (clawpatch_opt->dim == 2)
	{
		mx = clawpatch_opt->d2->mx;
		my = clawpatch_opt->d2->my;
	}
	else
	{
		mx = clawpatch_opt->d3->mx;
		my = clawpatch_opt->d3->my;
	}

	int meqn = clawpatch_opt->meqn;

	fclaw2d_clawpatch_registers_t *cr = *registers;  /* cr = clawpatch registers */

	cr = FCLAW_ALLOC(fclaw2d_clawpatch_registers_t,1);
	*registers = cr;

	/* Iterate over sides 0,1,3,4 */
	for(k = 0; k < 2; k++)
	{
		/* Accumulators */
		cr->fp[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);
		cr->fm[k]     = FCLAW_ALLOC_ZERO(double,my*meqn);

		cr->gp[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);
		cr->gm[k]     = FCLAW_ALLOC_ZERO(double,mx*meqn);

		cr->edge_fluxes[k]   = FCLAW_ALLOC_ZERO(double,2*my*meqn);
		cr->edge_fluxes[k+2] = FCLAW_ALLOC_ZERO(double,2*mx*meqn);

		cr->edgelengths[k]   = FCLAW_ALLOC(double,my);
		cr->edgelengths[k+2] = FCLAW_ALLOC(double,mx);
		cr->area[k]          = FCLAW_ALLOC(double,my);
		cr->area[k+2]        = FCLAW_ALLOC(double,mx);		
	}
}

void fclaw2d_clawpatch_time_sync_pack_registers(fclaw_global_t *glob,
                                                fclaw_patch_t *this_patch,
                                                double *qpack,
                                                int frsize, 
                                                fclaw_clawpatch_packmode_t packmode, 
                                                int *ierror)
{
	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);

	int mx;
	int my;

	if (clawpatch_opt->dim == 2)
	{
		mx = clawpatch_opt->d2->mx;
		my = clawpatch_opt->d2->my;
	}
	else
	{
		mx = clawpatch_opt->d3->mx;
		my = clawpatch_opt->d3->my;
	}

	int meqn = clawpatch_opt->meqn;

	fclaw2d_clawpatch_registers_t* cr = fclaw2d_clawpatch_get_registers(glob,this_patch);

	int cnt = 0;
    /* Cycle over four edges */
	for(int k = 0; k < 4; k++)
	{
		int idir = k/2;
		if (idir == 0)
		{
  		    /* k = 0,1 : Pack registers at faces 0,1 */
			for(int j = 0; j < meqn*my; j++)
			{
				if(packmode == CLAWPATCH_REGISTER_PACK)
				{
					qpack[cnt++] = cr->fm[k][j];
					qpack[cnt++] = cr->fp[k][j];
					qpack[cnt++] = cr->edge_fluxes[k][j];
					qpack[cnt++] = cr->edge_fluxes[k][j+meqn*my];
				}
				else
				{
					cr->fm[k][j] = qpack[cnt++];
					cr->fp[k][j] = qpack[cnt++];
					cr->edge_fluxes[k][j] = qpack[cnt++];
					cr->edge_fluxes[k][j+meqn*my] = qpack[cnt++]; 
				}
			}
			for(int j = 0; j < my; j++)
			{
				if (packmode == CLAWPATCH_REGISTER_PACK)
				{
					qpack[cnt++] = cr->edgelengths[k][j];
					qpack[cnt++] = cr->area[k][j];					
				}
				else
				{
					cr->edgelengths[k][j] = qpack[cnt++];
					cr->area[k][j] = qpack[cnt++];									
				}
			}
		}
		else if (idir == 1)
		{
            /* k = 2,3 : Pack registers at faces 2,3 */
			for(int i = 0; i < meqn*mx; i++)
			{
				if(packmode == CLAWPATCH_REGISTER_PACK)
				{
					qpack[cnt++] = cr->gm[k-2][i];
					qpack[cnt++] = cr->gp[k-2][i];
					qpack[cnt++] = cr->edge_fluxes[k][i];
					qpack[cnt++] = cr->edge_fluxes[k][i + meqn*mx];
				}
				else
				{
					cr->gm[k-2][i] = qpack[cnt++];
					cr->gp[k-2][i] = qpack[cnt++];
					cr->edge_fluxes[k][i] = qpack[cnt++];
					cr->edge_fluxes[k][i + meqn*mx] = qpack[cnt++];		
				}
			}
			for(int i = 0; i < mx; i++)
			{
				if (packmode == CLAWPATCH_REGISTER_PACK)
				{
					qpack[cnt++] = cr->edgelengths[k][i];
					qpack[cnt++] = cr->area[k][i];					
				}
				else
				{
					cr->edgelengths[k][i] = qpack[cnt++];
					cr->area[k][i] = qpack[cnt++];										
				}
			}
		}
	}
	/* Check that we packed exactly the right number of elements */
	*ierror = (cnt == frsize) ? 0 : 1;
}



void fclaw2d_clawpatch_time_sync_reset(fclaw_global_t *glob,
                                       fclaw_patch_t *this_patch, 
                                       int coarse_level,
                                       int reset_mode)
{
	int meqn;
	int i,j,k,idir;

	fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);

	fclaw2d_clawpatch_registers_t* cr = fclaw2d_clawpatch_get_registers(glob,this_patch);

	int mx;
	int my;

	if (clawpatch_opt->dim == 2)
	{
		mx = clawpatch_opt->d2->mx;
		my = clawpatch_opt->d2->my;
	}
	else
	{
		mx = clawpatch_opt->d3->mx;
		my = clawpatch_opt->d3->my;
	}
	
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
			/* Reset registers at interface between levels 
			   'coarse_level' and 'fine_level' */
			int is_coarse = (pdata->face_neighbors[k] == FCLAW_PATCH_HALFSIZE)
          			&& (this_patch->level == coarse_level);
		    int is_fine = (pdata->face_neighbors[k] == FCLAW_PATCH_DOUBLESIZE) &&
		          (this_patch->level == fine_level);
		    reset_flux = is_coarse || is_fine;
		}
		else if (reset_mode == FCLAW2D_TIME_SYNC_RESET_SAMESIZE)
		{
			/* Reset registers at interfaces between same size grids on coarse level */
			reset_flux = (pdata->face_neighbors[k] == FCLAW_PATCH_SAMESIZE) &&   
			      (this_patch->level == coarse_level);
		}
		else if (reset_mode == FCLAW2D_TIME_SYNC_RESET_PHYS)
		{
			/* Reset flux registers at physical boundaries (not actually used, 
			   but they are accumulated, so should be cleared out) */
			reset_flux = pdata->face_neighbors[k] == FCLAW_PATCH_BOUNDARY;
		}

		if (reset_flux)
		{
			if (idir == 0)
			{
				for(j = 0; j < meqn*my; j++)
				{
					/* k = 0,1 */
					cr->fm[k][j] = 0;
					cr->fp[k][j] = 0;
					cr->edge_fluxes[k][j] = 0;
					cr->edge_fluxes[k][j+meqn*my] = 0;  /* Two fluxes stored at each edge point */
				}
			}
			else
			{
				for(i = 0; i < meqn*mx; i++)
				{
					/* k = 2,3 */
					cr->gm[k-2][i] = 0;
					cr->gp[k-2][i] = 0;
					cr->edge_fluxes[k][i] = 0;
					cr->edge_fluxes[k][i + meqn*mx] = 0;
				}
			}
		}     
	}
}


void fclaw2d_clawpatch_time_sync_delete (fclaw2d_clawpatch_registers_t **registers)
{
	int k;

	fclaw2d_clawpatch_registers_t *cr = *registers;

	for(k = 0; k < 2; k++)
	{
		/* Accumulators */
		FCLAW_FREE(cr->fp[k]);
		FCLAW_FREE(cr->fm[k]);    

		FCLAW_FREE(cr->gp[k]);    
		FCLAW_FREE(cr->gm[k]);  

		/* COARSE GRID information */
		FCLAW_FREE(cr->edge_fluxes[k]); 
		FCLAW_FREE(cr->edge_fluxes[k+2]); 

		FCLAW_FREE(cr->edgelengths[k]);
		FCLAW_FREE(cr->area[k]);
		FCLAW_FREE(cr->edgelengths[k+2]);
		FCLAW_FREE(cr->area[k+2]);
	 }
	 FCLAW_FREE(*registers);
	 *registers = NULL;
}


void fclaw2d_clawpatch_time_sync_setup(fclaw_global_t* glob,
                                       fclaw_patch_t* this_patch,
                                       int blockno,int patchno)
{
	double *area, *edgelengths, *curvature;

	fclaw2d_clawpatch_registers_t *cr = 
	                       fclaw2d_clawpatch_get_registers(glob,this_patch);

	const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);


    int mx,my,mbc;
    double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_metric_scalar(glob,this_patch,
									&area, &edgelengths,&curvature);

    if (fclaw_opt->manifold)
    {
        /* only need area != NULL if manifold == 1 */
        FCLAW_ASSERT(area != NULL); 
    }
	CLAWPATCH_TIME_SYNC_SETUP(&mx,&my,&mbc,&dx,&dy,area,edgelengths,
	                          cr->area[0],cr->area[1],
	                          cr->area[2],cr->area[3],
	                          cr->edgelengths[0],cr->edgelengths[1],
	                          cr->edgelengths[2],cr->edgelengths[3],
	                          &fclaw_opt->manifold);
}


/* This is a patch call-back */
void fclaw2d_clawpatch_time_sync_f2c(fclaw_global_t* glob,
                                     fclaw_patch_t* coarse_patch,
                                     fclaw_patch_t* fine_patch,
                                     int coarse_blockno, int fine_blockno,
                                     int coarse_patchno, 
                                     int idir,
                                     int igrid,
                                     int iface_coarse,
                                     int time_interp,
                                     fclaw2d_patch_transform_data_t
                                     *transform_data)
{
	/* We don't correct time interpolated grids, since we assume that the time 
	interpolated average will already have correction information from times 
	n and n+1 patches.  But we do correct ghost patches, since corrections will be 
	needed for copying and interpolation.  */


    int meqn;
	double *qcoarse;
	fclaw_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

	double *qfine;
	fclaw_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

	int mx,my,mbc;
	double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_registers_t* crcoarse = 
	           fclaw2d_clawpatch_get_registers(glob,coarse_patch);

	fclaw2d_clawpatch_registers_t* crfine = 
	           fclaw2d_clawpatch_get_registers(glob,fine_patch);


  	/* create dummy fine grid to handle indexing between blocks */
	double *qneighbor_dummy = FCLAW_ALLOC_ZERO(double,meqn*(mx+4*mbc)*(my+4*mbc));

	int normal_match = fclaw2d_patch_normal_match(glob->domain, coarse_blockno, 
	                                              coarse_patchno, iface_coarse);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

	/* This function is defined in fc2d_clawpack4.6 and fc2d_clawpack5 */
	clawpatch_vt->d2->fort_time_sync_f2c(&mx,&my,&mbc,&meqn,&idir,&iface_coarse,
	                                     &coarse_blockno, &fine_blockno,
	                                     &normal_match,
										 crcoarse->area[0], crcoarse->area[1], 
										 crcoarse->area[2], crcoarse->area[3],
										 qcoarse,
										 crcoarse->fp[0],crcoarse->fm[1],
										 crcoarse->gp[0],crcoarse->gm[1],
										 crfine->fm[0],crfine->fp[1],
										 crfine->gm[0],crfine->gp[1],
										 crcoarse->edge_fluxes[0],
										 crcoarse->edge_fluxes[1],
										 crcoarse->edge_fluxes[2],
										 crcoarse->edge_fluxes[3],
										 crfine->edge_fluxes[0],
										 crfine->edge_fluxes[1],
										 crfine->edge_fluxes[2],
										 crfine->edge_fluxes[3],
										 qneighbor_dummy,
										 &transform_data);

	FCLAW_FREE(qneighbor_dummy);
}


void fclaw2d_clawpatch_time_sync_samesize (struct fclaw_global* glob,
                                           struct fclaw_patch* this_patch,
                                           struct fclaw_patch* neighbor_patch,
                                           int this_iface,int idir,
                                           fclaw2d_patch_transform_data_t
                                           *transform_data)
{
	/* We don't correct time interpolated grids, since we assume that the time 
	interpolated average will already have correction information from times 
	n and n+1 patches.  But we do correct ghost patches, since corrections will be 
	needed for copying and interpolation.  */

	fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

	double *qthis;
	int meqn;
	fclaw_clawpatch_soln_data(glob,this_patch,&qthis,&meqn);

	int mx,my,mbc;
	double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
	                            &xlower,&ylower,&dx,&dy);

	fclaw2d_clawpatch_registers_t* crthis = 
	fclaw2d_clawpatch_get_registers(glob,this_patch);

	fclaw2d_clawpatch_registers_t* crneighbor = 
	fclaw2d_clawpatch_get_registers(glob,neighbor_patch);

	/* create dummy fine grid to handle indexing between blocks */
	double *qneighbor_dummy = FCLAW_ALLOC_ZERO(double,meqn*(mx+2*mbc)*(my+2*mbc));

	/* Include this for debugging */
	int this_blockno, this_patchno,globnum,level;
	int neighbor_blockno, neighbor_patchno;
	fclaw2d_patch_get_info2(glob->domain,this_patch,&this_blockno, &this_patchno,
							&globnum,&level);

	fclaw2d_patch_get_info2(glob->domain,neighbor_patch,&neighbor_blockno, 
							&neighbor_patchno,&globnum,&level);

    /* This function is defined in fc2d_clawpack4.6 and fc2d_clawpack5  */

    /* Distribute 0.5 of each correction from each side */
	clawpatch_vt->d2->fort_time_sync_samesize(&mx,&my,&mbc,&meqn,&idir,&this_iface,
	                                          &this_blockno, &neighbor_blockno,
	                                          crthis->area[0], crthis->area[1], 
	                                          crthis->area[2], crthis->area[3],
	                                          qthis,
	                                          crthis->fp[0],crthis->fm[1],
	                                          crthis->gp[0],crthis->gm[1],
	                                          crneighbor->fm[0],crneighbor->fp[1],
	                                          crneighbor->gm[0],crneighbor->gp[1],
	                                          crthis->edge_fluxes[0],
	                                          crthis->edge_fluxes[1],
	                                          crthis->edge_fluxes[2],
	                                          crthis->edge_fluxes[3],
	                                          crneighbor->edge_fluxes[0],
	                                          crneighbor->edge_fluxes[1],
	                                          crneighbor->edge_fluxes[2],
	                                          crneighbor->edge_fluxes[3],
	                                          qneighbor_dummy,
	                                          &transform_data);    
	FCLAW_FREE(qneighbor_dummy);

}





