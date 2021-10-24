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

#include <fclaw2d_patch.h>

#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>

struct fclaw2d_patch_transform_data;

static fclaw2d_patch_vtable_t s_patch_vt;

/* ------------------------------- static access functions ---------------------------- */
static
fclaw2d_patch_data_t* get_patch_data(fclaw2d_patch_t* patch)
{
	fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) patch->user;
	FCLAW_ASSERT(pdata != NULL);
	return pdata;
}

static
void* get_user_patch(fclaw2d_patch_t* patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	return pdata->user_patch;
}



/* ----------------------------- Creating/deleting patches ---------------------------- */

static
void patch_data_new(fclaw2d_global_t* glob,
					fclaw2d_patch_t* this_patch,
					int this_block_idx, int this_patch_idx)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

	/* Initialize user data */
	fclaw2d_patch_data_t *pdata = FCLAW2D_ALLOC(fclaw2d_patch_data_t, 1);
	this_patch->user = (void *) pdata;

#if 0   
    /* This check is dubious, since glob->domain is the old domain */
#ifdef FCLAW_ENABLE_DEBUG
    fclaw2d_block_t *block = &glob->domain->blocks[this_block_idx];
#endif
    if (!(0 <= this_patch_idx && this_patch_idx < block->num_patches))
    {
    	// fclaw_global_essentialf("patch index is incorrect\n");
    	// FCLAW_ASSERT (0 <= this_patch_idx && this_patch_idx < block->num_patches);    
    }
#endif    

	pdata->patch_idx = this_patch_idx;
	pdata->block_idx = this_block_idx;

	/* create new user data */
	FCLAW_ASSERT(patch_vt->patch_new != NULL);
	pdata->user_patch = patch_vt->patch_new();

	fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(glob->domain);
	++ddata->count_set_patch; //this is now in cb_fclaw2d_regrid_repopulate 
	pdata->neighbors_set = 0;
}

void fclaw2d_patch_reset_data(fclaw2d_global_t* glob,
							  fclaw2d_patch_t* old_patch,
							  fclaw2d_patch_t* new_patch,
							  int blockno,int old_patchno, int new_patchno)
{
	fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) new_patch->user;
	pdata->block_idx = blockno;
#if 0 
    /* This check may be bogus, since glob->domain is the old domain */
#ifdef FCLAW_ENABLE_DEBUG
    fclaw2d_block_t *block = glob->domain->blocks + blockno;
#endif
    FCLAW_ASSERT (0 <= new_patchno && new_patchno < block->num_patches);    
#endif    

	pdata->patch_idx = new_patchno;

	/* Everything else will hopefully sort itself out ... */
}

void fclaw2d_patch_data_delete(fclaw2d_global_t *glob,
							   fclaw2d_patch_t *this_patch)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->patch_delete != NULL);
	fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) this_patch->user;

    if (pdata != NULL)
    {
        if (patch_vt->destroy_user_data)
        {
            patch_vt->destroy_user_data(glob,this_patch);
        }

        fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(glob->domain);        
        patch_vt->patch_delete(pdata->user_patch);
        ++ddata->count_delete_patch;

		FCLAW2D_FREE(pdata);
		this_patch->user = NULL;
	}
}

void fclaw2d_patch_build(fclaw2d_global_t *glob,
						 fclaw2d_patch_t *this_patch,
						 int blockno,
						 int patchno,
						 void *user)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

	/* This is where we store 'user' (from the point of view of a p4est user) data */
	patch_data_new(glob,this_patch,blockno, patchno);

	FCLAW_ASSERT(patch_vt->build != NULL);
	patch_vt->build(glob,
					this_patch,
					blockno,
					patchno,
					user);

#if 0
	fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
	pdata->blockno = blockno;
#endif    


    if (patch_vt->setup != NULL)
    {
        patch_vt->setup(glob,this_patch,blockno,patchno);
    }
    if (patch_vt->create_user_data)
    {
        patch_vt->create_user_data(glob,this_patch);
    }
}



void fclaw2d_patch_build_from_fine(fclaw2d_global_t *glob,
                                   fclaw2d_patch_t *fine_patches,
                                   fclaw2d_patch_t *coarse_patch,
                                   int blockno,
                                   int coarse_patchno,
                                   int fine0_patchno,
                                   fclaw2d_build_mode_t build_mode)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    FCLAW_ASSERT(patch_vt->build_from_fine != NULL);

    patch_data_new(glob,coarse_patch,blockno, coarse_patchno);

    fclaw2d_patch_data_t *pdata = get_patch_data(coarse_patch);
    pdata->block_idx = blockno;

    patch_vt->build_from_fine(glob,
                              fine_patches,
                              coarse_patch,
                              blockno,
                              coarse_patchno,
                              fine0_patchno,
                              build_mode);

    if (patch_vt->setup != NULL && build_mode == FCLAW2D_BUILD_FOR_UPDATE)
    {
        patch_vt->setup(glob,coarse_patch,blockno,coarse_patchno);
    }
    if (patch_vt->create_user_data)
    {
        patch_vt->create_user_data(glob,coarse_patch);
    }    
}

#if 0
void fclaw2d_patch_create_user_data(fclaw2d_global_t* glob,
                                    fclaw2d_patch_t* patch)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    FCLAW_ASSERT(patch_vt->create_user_data != NULL);
    patch_vt->create_user_data(glob,patch);
}

void fclaw2d_patch_destroy_user_data(fclaw2d_global_t* glob,
                                     fclaw2d_patch_t* patch)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    FCLAW_ASSERT(patch_vt->destroy_user_data != NULL);
    patch_vt->destroy_user_data(glob,patch);
}
#endif

#if 0
void* fclaw2d_patch_get_user_data(fclaw2d_global_t* glob,
                              fclaw2d_patch_t* this_patch)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    return pdata->user_data;
}

void fclaw2d_patch_set_user_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch, 
                                 void* user)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->user_data = user;
}
#endif


/* --------------------------- Solver specific functions ------------------------------ */

void fclaw2d_patch_initialize(fclaw2d_global_t *glob,
							  fclaw2d_patch_t *this_patch,
							  int this_block_idx,
							  int this_patch_idx)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    if (patch_vt->initialize != NULL)
    {
        patch_vt->initialize(glob,this_patch,this_block_idx,this_patch_idx);
    }
}



void fclaw2d_patch_physical_bc(fclaw2d_global_t *glob,
							   fclaw2d_patch_t *this_patch,
							   int this_block_idx,
							   int this_patch_idx,
							   double t,
							   double dt,
							   int *intersects_bc,
							   int time_interp)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->physical_bc != NULL);
	patch_vt->physical_bc(glob,this_patch,this_block_idx,this_patch_idx,
							t,dt,intersects_bc,time_interp);
}

double fclaw2d_patch_single_step_update(fclaw2d_global_t *glob,
                                        fclaw2d_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        double t,
                                        double dt, 
                                        void* user)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->single_step_update != NULL);

    double maxcfl = patch_vt->single_step_update(glob,this_patch,this_block_idx,
                                                   this_patch_idx,t,dt, user);
    return maxcfl;
}


void fclaw2d_patch_set_rhs(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *patch,
                           int blockno,
                           int patchno)
{
    fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
    
    FCLAW_ASSERT(patch_vt->rhs != NULL);

    patch_vt->rhs(glob,patch, blockno, patchno);
}


/* ------------------------------------ time stepping --------------------------------- */

void fclaw2d_patch_restore_step(fclaw2d_global_t* glob,
								fclaw2d_patch_t* this_patch)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->restore_step != NULL);

	patch_vt->restore_step(glob, this_patch);
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw2d_patch_save_step(fclaw2d_global_t* glob,
							 fclaw2d_patch_t* this_patch)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->save_step != NULL);
	patch_vt->save_step(glob, this_patch);
}

void fclaw2d_patch_setup_timeinterp(fclaw2d_global_t *glob,
									fclaw2d_patch_t *this_patch,
									double alpha)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->setup_timeinterp != NULL);

	patch_vt->setup_timeinterp(glob,this_patch,alpha);
}
	
/* ---------------------------------- Ghost filling  ---------------------------------- */

void fclaw2d_patch_copy_face(fclaw2d_global_t* glob,
							 fclaw2d_patch_t *this_patch,
							 fclaw2d_patch_t *neighbor_patch,
							 int iface,
							 int time_interp,
							 struct fclaw2d_patch_transform_data *transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->copy_face != NULL);

	patch_vt->copy_face(glob,this_patch,neighbor_patch,iface,
					   time_interp,transform_data);
}

void fclaw2d_patch_average_face(fclaw2d_global_t* glob,
								fclaw2d_patch_t *coarse_patch,
								fclaw2d_patch_t *fine_patch,
								int idir,
								int iface_coarse,
								int RefineFactor,
								int refratio,
								int time_interp,
								int igrid,
								struct fclaw2d_patch_transform_data* transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->average_face != NULL);

	patch_vt->average_face(glob,coarse_patch,fine_patch,idir,
						   iface_coarse,RefineFactor,
						   refratio,time_interp,igrid,
						   transform_data);
}

void fclaw2d_patch_interpolate_face(fclaw2d_global_t* glob,
									fclaw2d_patch_t *coarse_patch,
									fclaw2d_patch_t *fine_patch,
									int idir,
									int iside,
									int RefineFactor,
									int refratio,
									int time_interp,
									int igrid,
									struct fclaw2d_patch_transform_data* transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->interpolate_face != NULL);
	patch_vt->interpolate_face(glob,coarse_patch,fine_patch,idir,
								 iside,RefineFactor,refratio,
								 time_interp,igrid,transform_data);
}


void fclaw2d_patch_copy_corner(fclaw2d_global_t* glob,
							   fclaw2d_patch_t *this_patch,
							   fclaw2d_patch_t *corner_patch,
							   int coarse_blockno,
							   int fine_blockno,
							   int is_block_corner,
							   int icorner,
							   int time_interp,
							   struct fclaw2d_patch_transform_data *transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	if (!is_block_corner)
	{
		FCLAW_ASSERT(patch_vt->copy_corner != NULL);        
		patch_vt->copy_corner(glob,this_patch,corner_patch,
							  coarse_blockno,fine_blockno,
							  icorner,time_interp,transform_data);
	}
	else
	{
		FCLAW_ASSERT(patch_vt->copy_block_corner != NULL);        
		patch_vt->copy_block_corner(glob,this_patch,corner_patch,
									coarse_blockno,fine_blockno,
									icorner,time_interp,transform_data);        
	}
}

void fclaw2d_patch_average_corner(fclaw2d_global_t* glob,
								  fclaw2d_patch_t *coarse_patch,
								  fclaw2d_patch_t *fine_patch,
								  int coarse_blockno,
								  int fine_blockno,
								  int is_block_corner,
								  int coarse_corner,
								  int time_interp,
								  struct fclaw2d_patch_transform_data* transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	if (!is_block_corner)
	{
		FCLAW_ASSERT(patch_vt->average_corner != NULL);
		patch_vt->average_corner(glob,coarse_patch,fine_patch,
								 coarse_blockno,fine_blockno,
								 coarse_corner,
								 time_interp,transform_data);
	}
	else
	{
		FCLAW_ASSERT(patch_vt->average_block_corner != NULL);
		patch_vt->average_block_corner(glob,coarse_patch,fine_patch,
									   coarse_blockno,fine_blockno,
									   coarse_corner,
									   time_interp,transform_data);        
	}
}

void fclaw2d_patch_interpolate_corner(fclaw2d_global_t* glob,
									  fclaw2d_patch_t* coarse_patch,
									  fclaw2d_patch_t* fine_patch,
									  int coarse_blockno,
									  int fine_blockno,
									  int is_block_corner,
									  int coarse_corner,
									  int time_interp,
									  struct fclaw2d_patch_transform_data
									  *transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	if (!is_block_corner)
	{
		FCLAW_ASSERT(patch_vt->interpolate_corner != NULL);
		patch_vt->interpolate_corner(glob,coarse_patch,fine_patch,
									 coarse_blockno,fine_blockno,
									 coarse_corner,time_interp,
									 transform_data);        
	}
	else
	{
		FCLAW_ASSERT(patch_vt->interpolate_block_corner != NULL);
		patch_vt->interpolate_block_corner(glob,coarse_patch,fine_patch,
										   coarse_blockno,fine_blockno,
										   coarse_corner,time_interp,
										   transform_data);

	}
}



/* -------------------------- Transform functions ------------------------------------- */

void fclaw2d_patch_transform_init_data(struct fclaw2d_global* glob,
									   struct fclaw2d_patch* patch,
									   int blockno, int patchno,
									   struct fclaw2d_patch_transform_data *ftransform)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->transform_init_data != NULL);
	patch_vt->transform_init_data(glob,patch,blockno,patchno,ftransform);    
}

/* This is the transform at block faces.  */
void fclaw2d_patch_transform_blockface(int faceno, int rfaceno,
									   int ftransform[])
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->transform_face != NULL);
	patch_vt->transform_face(faceno, rfaceno, ftransform);
}

/* Transform within a block, where we only have the identity transform */
void fclaw2d_patch_transform_blockface_intra(int ftransform[])
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->transform_face_intra != NULL);
	patch_vt->transform_face_intra(ftransform);
}


/* ------------------------------- Regridding functions ------------------------------- */

int fclaw2d_patch_tag4refinement(fclaw2d_global_t *glob,
								 fclaw2d_patch_t *this_patch,
								 int blockno, int patchno,
								 int initflag)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->tag4refinement != NULL);

	return patch_vt->tag4refinement(glob,this_patch,blockno,
									patchno, initflag);
}

int fclaw2d_patch_tag4coarsening(fclaw2d_global_t *glob,
									  fclaw2d_patch_t *fine_patches,
									  int blockno,
									  int patchno,
                                      int initflag)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->tag4coarsening != NULL);

	return patch_vt->tag4coarsening(glob,fine_patches,
									blockno, patchno,initflag);
}

void fclaw2d_patch_average2coarse(fclaw2d_global_t *glob,
								  fclaw2d_patch_t *fine_patches,
								  fclaw2d_patch_t *coarse_patch,
								  int blockno, int fine0_patchno,
								  int coarse_patchno)

{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->average2coarse != NULL);

	patch_vt->average2coarse(glob,fine_patches,coarse_patch,
							 blockno,fine0_patchno,coarse_patchno);
}

void fclaw2d_patch_interpolate2fine(fclaw2d_global_t* glob,
									fclaw2d_patch_t* coarse_patch,
									fclaw2d_patch_t* fine_patches,
									int this_blockno, int coarse_patchno,
									int fine0_patchno)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->interpolate2fine != NULL);

	patch_vt->interpolate2fine(glob,coarse_patch,fine_patches,
							   this_blockno,coarse_patchno,
							   fine0_patchno);
}

/* ---------------------------- Ghost patches (local and remote) ---------------------- */

size_t fclaw2d_patch_ghost_packsize(fclaw2d_global_t* glob)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->ghost_packsize != NULL);
	return patch_vt->ghost_packsize(glob);
}

void fclaw2d_patch_local_ghost_alloc(fclaw2d_global_t* glob,
									 void** q)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	if (patch_vt->local_ghost_alloc != NULL)
	{
		patch_vt->local_ghost_alloc(glob, q);
	}
	else
	{
		/* This assumes that sizeof(char) = 1 */
		size_t psize = fclaw2d_patch_ghost_packsize(glob);
		*q = (void*) FCLAW_ALLOC(char,psize);  /* sizeof(datatype) included in psize */
		FCLAW_ASSERT(*q != NULL);
	}
}

void fclaw2d_patch_local_ghost_free(fclaw2d_global_t* glob,
									void **q)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	if (patch_vt->local_ghost_free != NULL)
	{
		patch_vt->local_ghost_free(glob,q);
	}
	else
	{
		FCLAW_FREE(*q);
		*q = NULL;
	}
}

void fclaw2d_patch_local_ghost_pack(fclaw2d_global_t *glob,
									fclaw2d_patch_t *this_patch,
									void *patch_data,
									int time_interp)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->local_ghost_pack != NULL);
	patch_vt->local_ghost_pack(glob,
							   this_patch,
							   patch_data,
							   time_interp);
}

void fclaw2d_patch_remote_ghost_build(fclaw2d_global_t *glob,
									  fclaw2d_patch_t *this_patch,
									  int blockno,
									  int patchno,
									  fclaw2d_build_mode_t build_mode)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

	FCLAW_ASSERT(patch_vt->remote_ghost_build != NULL);

	patch_data_new(glob,this_patch,blockno,patchno);

	patch_vt->remote_ghost_build(glob,this_patch,blockno,
							patchno,build_mode);
	if (patch_vt->remote_ghost_setup != NULL)
	{
		patch_vt->remote_ghost_setup(glob,this_patch,blockno,patchno);
	}
}


void fclaw2d_patch_remote_ghost_delete(fclaw2d_global_t *glob,
									   fclaw2d_patch_t *this_patch)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	
	fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);


	if (pdata != NULL)
	{
		FCLAW_ASSERT(patch_vt->patch_delete != NULL);
		patch_vt->remote_ghost_delete(pdata->user_patch);

		FCLAW2D_FREE(pdata);
		this_patch->user = NULL;

		fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data(glob->domain);
		++ddata->count_delete_patch;
	}
}

void fclaw2d_patch_remote_ghost_unpack(fclaw2d_global_t* glob,
									   fclaw2d_patch_t* this_patch,
									   int this_block_idx,
									   int this_patch_idx,
									   void *qdata, int time_interp)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->remote_ghost_unpack != NULL);
	patch_vt->remote_ghost_unpack(glob, this_patch, this_block_idx,
								  this_patch_idx, qdata, time_interp);
}


/* ----------------------------------- Partitioning ----------------------------------- */

size_t fclaw2d_patch_partition_packsize(fclaw2d_global_t* glob)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->partition_packsize != NULL);

	return patch_vt->partition_packsize(glob);
}

void fclaw2d_patch_partition_pack(fclaw2d_global_t *glob,
								  fclaw2d_patch_t *this_patch,
								  int this_block_idx,
								  int this_patch_idx,
								  void* pack_data_here)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->partition_pack != NULL);

	patch_vt->partition_pack(glob,
							 this_patch,
							 this_block_idx,
							 this_patch_idx,
							 pack_data_here);
}


void fclaw2d_patch_partition_unpack(fclaw2d_global_t *glob,  
									fclaw2d_domain_t *new_domain,
									fclaw2d_patch_t *this_patch,
									int this_block_idx,
									int this_patch_idx,
									void *unpack_data_from_here)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();

	fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;

	fclaw2d_patch_build(glob,this_patch,this_block_idx,
						this_patch_idx,(void*) &build_mode);

	/* This copied q data from memory */
	FCLAW_ASSERT(patch_vt->partition_unpack != NULL);

	patch_vt->partition_unpack(glob,  /* contains old domain */
							   new_domain,
							   this_patch,
							   this_block_idx,
							   this_patch_idx,
							   unpack_data_from_here);
}

/* ----------------------------- Conservative updates --------------------------------- */

/* We need to virtualize this because we call it from fclaw2d_face_neighbors */

void fclaw2d_patch_time_sync_f2c(fclaw2d_global_t* glob,
								 fclaw2d_patch_t *coarse_patch,
								 fclaw2d_patch_t *fine_patch,
								 int coarse_blockno, int fine_blockno,
								 int coarse_patchno, 
								 int idir,
								 int igrid,
								 int iface_coarse,
								 int time_interp,
								 fclaw2d_patch_transform_data_t* transform_data)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->time_sync_f2c != NULL);

	patch_vt->time_sync_f2c(glob,coarse_patch,fine_patch,
	                        coarse_blockno, fine_blockno,
	                        coarse_patchno, 
	                        idir, igrid,iface_coarse,time_interp,
	                        transform_data);    
}

/* Correct for metric discontinuities at block boundaries */
void fclaw2d_patch_time_sync_samesize(fclaw2d_global_t* glob,
                                      fclaw2d_patch_t *this_patch,
                                      fclaw2d_patch_t *neighbor_patch,
                                      int iface, int idir,
                                      struct fclaw2d_patch_transform_data *transform_data)

{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->time_sync_samesize != NULL);

	patch_vt->time_sync_samesize(glob,this_patch,neighbor_patch,iface,idir,
	                             transform_data);    
}

void fclaw2d_patch_time_sync_reset(fclaw2d_global_t* glob,
                                   fclaw2d_patch_t* this_patch,
                                   int coarse_level,
                                   int reset_mode)
{
	fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt();
	FCLAW_ASSERT(patch_vt->time_sync_reset != NULL);

	patch_vt->time_sync_reset(glob,this_patch,coarse_level, reset_mode);

}

/* ----------------------------------- Virtual table ---------------------------------- */

static
fclaw2d_patch_vtable_t* patch_vt_init()
{
	s_patch_vt.is_set = 0;
	return &s_patch_vt;
}

void fclaw2d_patch_vtable_initialize()
{
	fclaw2d_patch_vtable_t *patch_vt = patch_vt_init();

#if 0
	/* Function pointers are set to NULL by default so are not set here */

	/* These must be redefined by the solver and user */
	patch_vt->initialize          = NULL;
	patch_vt->physical_bc         = NULL;
	patch_vt->single_step_update  = NULL;

	patch_vt->time_sync_f2c       = NULL;
	patch_vt->time_sync_samesize  = NULL;

	/* These are optional */
	patch_vt->setup               = NULL;
	patch_vt->remote_ghost_setup  = NULL;
#endif    

	patch_vt->is_set = 1;

}

/* ------------------------------ User access functions ------------------------------- */

fclaw2d_patch_vtable_t* fclaw2d_patch_vt()
{
	FCLAW_ASSERT(s_patch_vt.is_set != 0);
	return &s_patch_vt;
}



/* Use this one if you have the patch or block number */
void fclaw2d_patch_get_info(fclaw2d_domain_t * domain,
							fclaw2d_patch_t * patch,
							int blockno, int patchno,
							int *global_num, int *local_num, 
                            int *level)

{
	fclaw2d_block_t *block = &domain->blocks[blockno];

    /* For blockno = 0, we have 
               patchno == local_num.
       For blockno > 0, the patchno resets on each new block and we have
               patchno + block->num_patches_before = local_num.  */

    *local_num = block->num_patches_before + patchno;

	*global_num = domain->global_num_patches_before + *local_num;

	*level = patch->level;

}

/* Use this one if you don't have the patch or block number */
void fclaw2d_patch_get_info2(fclaw2d_domain_t * domain,
							fclaw2d_patch_t * this_patch,
							int *this_block_idx, int *this_patch_idx,
							int *global_num, int *level)

{

    /* I don't completely trust this that blockno and patchno are consistent
       with p4est numbering. */
	fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t*) this_patch->user;
	*this_patch_idx = pdata->patch_idx;
	*this_block_idx = pdata->block_idx;

	fclaw2d_block_t *this_block = &domain->blocks[*this_block_idx];

	*global_num = domain->global_num_patches_before +
		(this_block->num_patches_before + *this_patch_idx);

	*level = this_patch->level;

}


void*
fclaw2d_patch_get_user_patch(fclaw2d_patch_t* patch)

{
	return get_user_patch(patch);
}

fclaw2d_patch_data_t*
fclaw2d_patch_get_patch_data(fclaw2d_patch_t* patch)
{
	return get_patch_data(patch);
}

void* fclaw2d_patch_get_user_data(fclaw2d_global_t* glob,
                                  fclaw2d_patch_t* this_patch)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    return pdata->user_data;
}

void fclaw2d_patch_set_user_data(fclaw2d_global_t* glob,
                                 fclaw2d_patch_t* this_patch, 
                                 void* user)
{
    fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->user_data = user;
}


void* fclaw2d_patch_metric_patch(struct fclaw2d_patch *patch)
{
	fclaw2d_patch_vtable_t* patch_vt = fclaw2d_patch_vt();    
	FCLAW_ASSERT(patch_vt->metric_patch != NULL);
	return patch_vt->metric_patch(patch);
}

int
fclaw2d_patch_get_blockno(fclaw2d_patch_t* this_patch)
{
	/* This is is really only here if we fail to pass in the blockno to a patch 
	routine.  */
	fclaw2d_patch_data_t* pdata = get_patch_data(this_patch);
	FCLAW_ASSERT(pdata->block_idx >= 0);
	return pdata->block_idx;
}

int
fclaw2d_patch_get_patchno(fclaw2d_patch_t* this_patch)
{
	/* This is is really only here if we fail to pass in the patchno to a patch 
	routine.  */
	fclaw2d_patch_data_t* pdata = get_patch_data(this_patch);
	FCLAW_ASSERT(pdata->patch_idx >= 0);
	return pdata->patch_idx;
}


/* -------------------------- Internal ForestClaw functions --------------------------- */

void fclaw2d_patch_set_face_type(fclaw2d_patch_t *patch,int iface,
								 fclaw2d_patch_relation_t face_type)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	pdata->face_neighbors[iface] = face_type;
}

void fclaw2d_patch_set_corner_type(fclaw2d_patch_t *patch,int icorner,
								   fclaw2d_patch_relation_t corner_type)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	pdata->corner_neighbors[icorner] = corner_type;
	pdata->corners[icorner] = 1;
}

void fclaw2d_patch_set_missing_corner(fclaw2d_patch_t *patch,int icorner)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	pdata->corners[icorner] = 0;
}

fclaw2d_patch_relation_t fclaw2d_patch_get_face_type(fclaw2d_patch_t* patch,
													 int iface)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	FCLAW_ASSERT(pdata->neighbors_set != 0);
	FCLAW_ASSERT(0 <= iface && iface < 4);
	return pdata->face_neighbors[iface];
}

fclaw2d_patch_relation_t fclaw2d_patch_get_corner_type(fclaw2d_patch_t* patch,
													   int icorner)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	FCLAW_ASSERT(pdata->corners[icorner] != 0);
	FCLAW_ASSERT(pdata->neighbors_set != 0);
	return pdata->corner_neighbors[icorner];
}

int fclaw2d_patch_corner_is_missing(fclaw2d_patch_t* patch,
									int icorner)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	return !pdata->corners[icorner];
}

void fclaw2d_patch_neighbors_set(fclaw2d_patch_t* patch)
{
	int iface, icorner;
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	FCLAW_ASSERT(pdata->neighbors_set == 0);

	pdata->has_finegrid_neighbors = 0;
	pdata->on_coarsefine_interface = 0;
	for (iface = 0; iface < 4; iface++)
	{
		fclaw2d_patch_relation_t nt;
		nt = pdata->face_neighbors[iface];
		if (nt == FCLAW2D_PATCH_HALFSIZE || (nt == FCLAW2D_PATCH_DOUBLESIZE))
		{
			pdata->on_coarsefine_interface = 1;
			if (nt == FCLAW2D_PATCH_HALFSIZE)
			{
				pdata->has_finegrid_neighbors = 1;
			}
		}
	}

	for (icorner = 0; icorner < 4; icorner++)
	{
		fclaw2d_patch_relation_t nt;
		int has_corner = pdata->corners[icorner];
		if (has_corner)
		{
			nt = pdata->corner_neighbors[icorner];
			if ((nt == FCLAW2D_PATCH_HALFSIZE) || (nt == FCLAW2D_PATCH_DOUBLESIZE))
			{
				pdata->on_coarsefine_interface = 1;
				if (nt == FCLAW2D_PATCH_HALFSIZE)
				{
					pdata->has_finegrid_neighbors = 1;
				}
			}
		}
	}
	pdata->neighbors_set = 1;
}

void fclaw2d_patch_neighbors_reset(fclaw2d_patch_t* patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	pdata->neighbors_set = 0;
}

int fclaw2d_patch_neighbor_type_set(fclaw2d_patch_t* patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	return pdata->neighbors_set;
}


int fclaw2d_patch_has_finegrid_neighbors(fclaw2d_patch_t *patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	return pdata->has_finegrid_neighbors;
}

int fclaw2d_patch_on_coarsefine_interface(fclaw2d_patch_t *patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	return pdata->on_coarsefine_interface;
}

int
fclaw2d_patch_on_parallel_boundary (const fclaw2d_patch_t * patch)
{
	return patch->flags & FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY ? 1 : 0;
}

int* fclaw2d_patch_block_corner_count(fclaw2d_global_t* glob,
									  fclaw2d_patch_t* this_patch)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
	return pdata->block_corner_count;
}

void fclaw2d_patch_set_block_corner_count(fclaw2d_global_t* glob,
										  fclaw2d_patch_t* this_patch,
										  int icorner, int block_corner_count)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(this_patch);
	pdata->block_corner_count[icorner] = block_corner_count;
}




