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

#include "forestclaw.h"
#include <fclaw_pointer_map.h>
#include <fclaw_global.h>
#include <fclaw_domain.h>

#include <fclaw_patch.h>
#include <fclaw2d_defs.h>
#include <fclaw3d_defs.h>

struct fclaw_patch_transform_data;

fclaw_patch_transform_data_t* fclaw_patch_transform_data_new(int dim)
{
    fclaw_patch_transform_data_t* tdata = FCLAW_ALLOC_ZERO(fclaw_patch_transform_data_t,1);
    tdata->dim = dim;
    return tdata;
}

void fclaw_patch_transform_data_destroy(fclaw_patch_transform_data_t *tdata)
{
    FCLAW_FREE(tdata);
}

/* ------------------------------- static access functions ---------------------------- */
static
fclaw_patch_data_t* get_patch_data(fclaw_patch_t* patch)
{
    fclaw_patch_data_t *pdata = (fclaw_patch_data_t*) patch->user;
    FCLAW_ASSERT(pdata != NULL);
    return pdata;
}

static
void* get_user_patch(fclaw_patch_t* patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    return pdata->user_patch;
}



/* ----------------------------- Creating/deleting patches ---------------------------- */

static
void patch_data_new(fclaw_global_t* glob,
                    fclaw_patch_t* this_patch,
                    int this_block_idx, int this_patch_idx)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    /* Initialize user data */
    fclaw_patch_data_t *pdata = FCLAW_ALLOC_ZERO(fclaw_patch_data_t, 1);
    if(glob->domain->dim == 2)
    {
        pdata->dim = 2;
        pdata->d2 = FCLAW_ALLOC_ZERO(fclaw_patch_data_d2_t, 1);
    }
    else 
    {
        pdata->dim = 3;
        pdata->d3 = FCLAW_ALLOC_ZERO(fclaw_patch_data_d3_t, 1);
    }
    this_patch->user = (void *) pdata;

    pdata->patch_idx = this_patch_idx;
    pdata->block_idx = this_block_idx;

    /* create new user data */
    FCLAW_ASSERT(patch_vt->patch_new != NULL);
    pdata->user_patch = patch_vt->patch_new();

    fclaw_domain_data_t *ddata = glob->domain->domain_data;
    ++ddata->count_set_patch; //this is now in cb_fclaw_regrid_repopulate 
    pdata->neighbors_set = 0;
}

void fclaw_patch_reset_data(fclaw_global_t* glob,
                              fclaw_patch_t* old_patch,
                              fclaw_patch_t* new_patch,
                              int blockno,int old_patchno, int new_patchno)
{
    fclaw_patch_data_t *pdata = (fclaw_patch_data_t*) new_patch->user;
    pdata->block_idx = blockno;

    pdata->patch_idx = new_patchno;

    /* Everything else will hopefully sort itself out ... */
}

void fclaw_patch_data_delete(fclaw_global_t *glob,
                               fclaw_patch_t *this_patch)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->patch_delete != NULL);
    fclaw_patch_data_t *pdata = (fclaw_patch_data_t*) this_patch->user;

    if (pdata != NULL)
    {
        if (patch_vt->destroy_user_data)
        {
            patch_vt->destroy_user_data(glob,this_patch);
        }

        patch_vt->patch_delete(pdata->user_patch);
        fclaw_domain_data_t *ddata = glob->domain->domain_data;
        ++ddata->count_delete_patch;

        if(pdata->dim == 2)
        {
            FCLAW_FREE(pdata->d2);
        }
        else 
        {
            FCLAW_FREE(pdata->d3);
        }
        FCLAW_FREE(pdata);
        this_patch->user = NULL;
    }
}

void fclaw_patch_build(fclaw_global_t *glob,
                         fclaw_patch_t *this_patch,
                         int blockno,
                         int patchno,
                         void *user)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    /* This is where we store 'user' (from the point of view of a p4est user) data */
    patch_data_new(glob,this_patch,blockno, patchno);

    FCLAW_ASSERT(patch_vt->build != NULL);
    patch_vt->build(glob,
                    this_patch,
                    blockno,
                    patchno,
                    user);

#if 0
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->blockno = blockno;
#endif    

    /* This should go first, in case the setup function relies on this 
       user data (or just needs to set it up)
    */
    if (patch_vt->create_user_data != NULL)
    {
        patch_vt->create_user_data(glob,this_patch);
    }

    if (patch_vt->setup != NULL)
    {
        patch_vt->setup(glob,this_patch,blockno,patchno);
    }
}



void fclaw_patch_build_from_fine(fclaw_global_t *glob,
                                   fclaw_patch_t *fine_patches,
                                   fclaw_patch_t *coarse_patch,
                                   int blockno,
                                   int coarse_patchno,
                                   int fine0_patchno,
                                   fclaw_build_mode_t build_mode)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->build_from_fine != NULL);

    patch_data_new(glob,coarse_patch,blockno, coarse_patchno);

    fclaw_patch_data_t *pdata = get_patch_data(coarse_patch);
    pdata->block_idx = blockno;

    patch_vt->build_from_fine(glob,
                              fine_patches,
                              coarse_patch,
                              blockno,
                              coarse_patchno,
                              fine0_patchno,
                              build_mode);

    if (patch_vt->create_user_data != NULL)
    {
        patch_vt->create_user_data(glob,coarse_patch);
    }    

    if (patch_vt->setup != NULL && build_mode == FCLAW_BUILD_FOR_UPDATE)
    {
        patch_vt->setup(glob,coarse_patch,blockno,coarse_patchno);
    }
}

/* --------------------------- Solver specific functions ------------------------------ */

void fclaw_patch_initialize(fclaw_global_t *glob,
                              fclaw_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if (patch_vt->initialize != NULL)
    {
        patch_vt->initialize(glob,this_patch,this_block_idx,this_patch_idx);
    }
}



void fclaw_patch_physical_bc(fclaw_global_t *glob,
                               fclaw_patch_t *this_patch,
                               int this_block_idx,
                               int this_patch_idx,
                               double t,
                               double dt,
                               int *intersects_bc,
                               int time_interp)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->physical_bc != NULL);
    patch_vt->physical_bc(glob,this_patch,this_block_idx,this_patch_idx,
                            t,dt,intersects_bc,time_interp);
}

double fclaw_patch_single_step_update(fclaw_global_t *glob,
                                        fclaw_patch_t *this_patch,
                                        int this_block_idx,
                                        int this_patch_idx,
                                        double t,
                                        double dt, 
                                        void* user)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->single_step_update != NULL);

    double maxcfl = patch_vt->single_step_update(glob,this_patch,this_block_idx,
                                                   this_patch_idx,t,dt, user);
    return maxcfl;
}


void fclaw_patch_set_rhs(fclaw_global_t *glob,
                           fclaw_patch_t *patch,
                           int blockno,
                           int patchno)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    
    FCLAW_ASSERT(patch_vt->rhs != NULL);

    patch_vt->rhs(glob,patch, blockno, patchno);
}


/* ------------------------------------ time stepping --------------------------------- */

void fclaw_patch_restore_step(fclaw_global_t* glob,
                                fclaw_patch_t* this_patch)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->restore_step != NULL);

    patch_vt->restore_step(glob, this_patch);
}

/* This is called from libraries routines (clawpack4.6, clawpack5, etc) */
void fclaw_patch_save_step(fclaw_global_t* glob,
                             fclaw_patch_t* this_patch)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->save_step != NULL);
    patch_vt->save_step(glob, this_patch);
}

void fclaw_patch_setup_timeinterp(fclaw_global_t *glob,
                                    fclaw_patch_t *this_patch,
                                    double alpha)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->setup_timeinterp != NULL);

    patch_vt->setup_timeinterp(glob,this_patch,alpha);
}
    
/* ---------------------------------- Ghost filling  ---------------------------------- */

void fclaw_patch_copy_face(fclaw_global_t* glob,
                             fclaw_patch_t *this_patch,
                             fclaw_patch_t *neighbor_patch,
                             int iface,
                             int time_interp,
                             fclaw_patch_transform_data_t *transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    if(patch_vt->dim == 2)
    {
        FCLAW_ASSERT(patch_vt->d2->copy_face != NULL);
        patch_vt->d2->copy_face(glob,this_patch,neighbor_patch,iface,
                                time_interp,transform_data);
    }
    else 
    {
        FCLAW_ASSERT(patch_vt->d3->copy_face != NULL);
        patch_vt->d3->copy_face(glob,this_patch,neighbor_patch,iface,
                                time_interp,transform_data);
    }
}

void fclaw_patch_average_face(fclaw_global_t* glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patch,
                              int idir,
                              int iface_coarse,
                              int RefineFactor,
                              int refratio,
                              int time_interp,
                              int igrid,
                              fclaw_patch_transform_data_t* transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    if(patch_vt->dim == 2)
    {
        FCLAW_ASSERT(patch_vt->d2->average_face != NULL);
        patch_vt->d2->average_face(glob,coarse_patch,fine_patch,idir,
                                   iface_coarse,RefineFactor,
                                   refratio,time_interp,igrid,
                                   transform_data);
    }
    else 
    {
        FCLAW_ASSERT(patch_vt->d3->average_face != NULL);
        patch_vt->d3->average_face(glob,coarse_patch,fine_patch,idir,
                                   iface_coarse,RefineFactor,
                                   refratio,time_interp,igrid,
                                   transform_data);
    }
}

void fclaw_patch_interpolate_face(fclaw_global_t* glob,
                                    fclaw_patch_t *coarse_patch,
                                    fclaw_patch_t *fine_patch,
                                    int idir,
                                    int iside,
                                    int RefineFactor,
                                    int refratio,
                                    int time_interp,
                                    int igrid,
                                    fclaw_patch_transform_data_t* transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if(patch_vt->dim == 2)
    {
        FCLAW_ASSERT(patch_vt->d2->interpolate_face != NULL);
        patch_vt->d2->interpolate_face(glob,coarse_patch,fine_patch,idir,
                                       iside,RefineFactor,refratio,
                                       time_interp,igrid,transform_data);
    }
    else 
    {
        FCLAW_ASSERT(patch_vt->d3->interpolate_face != NULL);
        patch_vt->d3->interpolate_face(glob,coarse_patch,fine_patch,idir,
                                       iside,RefineFactor,refratio,
                                       time_interp,igrid,transform_data);
    }
}

void fclaw_patch_copy_edge(fclaw_global_t* glob,
                           fclaw_patch_t *this_patch,
                           fclaw_patch_t *edge_patch,
                           int coarse_blockno,
                           int fine_blockno,
                           int icorner,
                           int time_interp,
                           fclaw_patch_transform_data_t *transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
	FCLAW_ASSERT(patch_vt->dim == 3);
    FCLAW_ASSERT(patch_vt->d3->copy_edge != NULL);        
    patch_vt->d3->copy_edge(glob,this_patch,edge_patch,
                            coarse_blockno,fine_blockno,
                            icorner,time_interp,transform_data);
}

void fclaw_patch_average_edge(fclaw_global_t* glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patch,
                              int iedge_coarse,
                              int time_interp,
                              fclaw_patch_transform_data_t* transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
	FCLAW_ASSERT(patch_vt->dim == 3);
    FCLAW_ASSERT(patch_vt->d3->copy_edge != NULL);        
    patch_vt->d3->average_edge(glob,coarse_patch,fine_patch,
                               iedge_coarse,
                               time_interp,transform_data);

}

void fclaw_patch_copy_corner(fclaw_global_t* glob,
                               fclaw_patch_t *this_patch,
                               fclaw_patch_t *corner_patch,
                               int coarse_blockno,
                               int fine_blockno,
                               int is_block_corner,
                               int icorner,
                               int time_interp,
                               fclaw_patch_transform_data_t *transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if(patch_vt->dim == 2)
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d2->copy_corner != NULL);        
            patch_vt->d2->copy_corner(glob,this_patch,corner_patch,
                                      coarse_blockno,fine_blockno,
                                      icorner,time_interp,transform_data);
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d2->copy_block_corner != NULL);        
            patch_vt->d2->copy_block_corner(glob,this_patch,corner_patch,
                                            coarse_blockno,fine_blockno,
                                            icorner,time_interp,transform_data);        
        }
    }
    else 
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d3->copy_corner != NULL);        
            patch_vt->d3->copy_corner(glob,this_patch,corner_patch,
                                      coarse_blockno,fine_blockno,
                                      icorner,time_interp,transform_data);
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d3->copy_block_corner != NULL);        
            patch_vt->d3->copy_block_corner(glob,this_patch,corner_patch,
                                            coarse_blockno,fine_blockno,
                                            icorner,time_interp,transform_data);        
        }
    }
}

void fclaw_patch_average_corner(fclaw_global_t* glob,
                                  fclaw_patch_t *coarse_patch,
                                  fclaw_patch_t *fine_patch,
                                  int coarse_blockno,
                                  int fine_blockno,
                                  int is_block_corner,
                                  int coarse_corner,
                                  int time_interp,
                                  fclaw_patch_transform_data_t* transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if(patch_vt->dim == 2)
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d2->copy_corner != NULL);        
            patch_vt->d2->average_corner(glob,coarse_patch,fine_patch,
                                      coarse_blockno,fine_blockno,
                                      coarse_corner,time_interp,transform_data);
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d2->copy_block_corner != NULL);        
            patch_vt->d2->average_block_corner(glob,coarse_patch,fine_patch,
                                               coarse_blockno,fine_blockno,
                                               coarse_corner,time_interp,transform_data);        
        }
    }
    else 
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d3->copy_corner != NULL);        
            patch_vt->d3->average_corner(glob,coarse_patch,fine_patch,
                                         coarse_blockno,fine_blockno,
                                         coarse_corner,time_interp,transform_data);
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d3->copy_block_corner != NULL);        
            patch_vt->d3->average_block_corner(glob,coarse_patch,fine_patch,
                                               coarse_blockno,fine_blockno,
                                               coarse_corner,time_interp,transform_data);        
        }
    }
}

void fclaw_patch_interpolate_corner(fclaw_global_t* glob,
                                      fclaw_patch_t* coarse_patch,
                                      fclaw_patch_t* fine_patch,
                                      int coarse_blockno,
                                      int fine_blockno,
                                      int is_block_corner,
                                      int coarse_corner,
                                      int time_interp,
                                      fclaw_patch_transform_data_t
                                      *transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if(patch_vt->dim == 2)
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d2->interpolate_corner != NULL);
            patch_vt->d2->interpolate_corner(glob,coarse_patch,fine_patch,
                                             coarse_blockno,fine_blockno,
                                             coarse_corner,time_interp,
                                             transform_data);        
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d2->interpolate_block_corner != NULL);
            patch_vt->d2->interpolate_block_corner(glob,coarse_patch,fine_patch,
                                                   coarse_blockno,fine_blockno,
                                                   coarse_corner,time_interp,
                                                   transform_data);
        }
    }
    else
    {
        if (!is_block_corner)
        {
            FCLAW_ASSERT(patch_vt->d3->interpolate_corner != NULL);
            patch_vt->d3->interpolate_corner(glob,coarse_patch,fine_patch,
                                             coarse_blockno,fine_blockno,
                                             coarse_corner,time_interp,
                                             transform_data);        
        }
        else
        {
            FCLAW_ASSERT(patch_vt->d3->interpolate_block_corner != NULL);
            patch_vt->d3->interpolate_block_corner(glob,coarse_patch,fine_patch,
                                                   coarse_blockno,fine_blockno,
                                                   coarse_corner,time_interp,
                                                   transform_data);
        }
    }
}



/* -------------------------- Transform functions ------------------------------------- */

void fclaw_patch_transform_init_data(fclaw_global_t* glob,
                                       fclaw_patch_t* patch,
                                       int blockno, int patchno,
                                       fclaw_patch_transform_data_t *ftransform)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->transform_init_data != NULL);
    patch_vt->transform_init_data(glob,patch,blockno,patchno,ftransform);    
}

/* This is the transform at block faces.  */
void fclaw_patch_transform_blockface(fclaw_global_t* glob,
                                       int faceno, int rfaceno,
                                       int ftransform[])
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->transform_face != NULL);
    patch_vt->transform_face(faceno, rfaceno, ftransform);
}

/* Transform within a block, where we only have the identity transform */
void fclaw_patch_transform_blockface_intra(fclaw_global_t* glob,
                                             int ftransform[])
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->transform_face_intra != NULL);
    patch_vt->transform_face_intra(ftransform);
}


/* ------------------------------- Regridding functions ------------------------------- */

int fclaw_patch_tag4refinement(fclaw_global_t *glob,
                                 fclaw_patch_t *this_patch,
                                 int blockno, int patchno,
                                 int initflag)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->tag4refinement != NULL);

    return patch_vt->tag4refinement(glob,this_patch,blockno,
                                    patchno, initflag);
}

int fclaw_patch_tag4coarsening(fclaw_global_t *glob,
                                      fclaw_patch_t *fine_patches,
                                      int blockno,
                                      int patchno,
                                      int initflag)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->tag4coarsening != NULL);

    return patch_vt->tag4coarsening(glob,fine_patches,
                                    blockno, patchno,initflag);
}

void fclaw_patch_average2coarse(fclaw_global_t *glob,
                                  fclaw_patch_t *fine_patches,
                                  fclaw_patch_t *coarse_patch,
                                  int blockno, int fine0_patchno,
                                  int coarse_patchno)

{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->average2coarse != NULL);

    patch_vt->average2coarse(glob,fine_patches,coarse_patch,
                             blockno,fine0_patchno,coarse_patchno);
}

void fclaw_patch_interpolate2fine(fclaw_global_t* glob,
                                    fclaw_patch_t* coarse_patch,
                                    fclaw_patch_t* fine_patches,
                                    int this_blockno, int coarse_patchno,
                                    int fine0_patchno)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->interpolate2fine != NULL);

    patch_vt->interpolate2fine(glob,coarse_patch,fine_patches,
                               this_blockno,coarse_patchno,
                               fine0_patchno);
}

/* ---------------------------- Ghost patches (local and remote) ---------------------- */

size_t fclaw_patch_ghost_packsize(fclaw_global_t* glob)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->ghost_packsize != NULL);
    return patch_vt->ghost_packsize(glob);
}

void fclaw_patch_local_ghost_alloc(fclaw_global_t* glob,
                                     void** q)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    if (patch_vt->local_ghost_alloc != NULL)
    {
        patch_vt->local_ghost_alloc(glob, q);
    }
    else
    {
        /* This assumes that sizeof(char) = 1 */
        size_t psize = fclaw_patch_ghost_packsize(glob);
        *q = (void*) FCLAW_ALLOC(char,psize);  /* sizeof(datatype) included in psize */
        FCLAW_ASSERT(*q != NULL);
    }
}

void fclaw_patch_local_ghost_free(fclaw_global_t* glob,
                                    void **q)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
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

void fclaw_patch_local_ghost_pack(fclaw_global_t *glob,
                                    fclaw_patch_t *this_patch,
                                    void *patch_data,
                                    int time_interp)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->local_ghost_pack != NULL);
    patch_vt->local_ghost_pack(glob,
                               this_patch,
                               patch_data,
                               time_interp);
}

void fclaw_patch_remote_ghost_build(fclaw_global_t *glob,
                                      fclaw_patch_t *this_patch,
                                      int blockno,
                                      int patchno,
                                      fclaw_build_mode_t build_mode)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    FCLAW_ASSERT(patch_vt->remote_ghost_build != NULL);

    patch_data_new(glob,this_patch,blockno,patchno);

    patch_vt->remote_ghost_build(glob,this_patch,blockno,
                            patchno,build_mode);
    if (patch_vt->remote_ghost_setup != NULL)
    {
        patch_vt->remote_ghost_setup(glob,this_patch,blockno,patchno);
    }
}


void fclaw_patch_remote_ghost_delete(fclaw_global_t *glob,
                                       fclaw_patch_t *this_patch)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);


    if (pdata != NULL)
    {
        FCLAW_ASSERT(patch_vt->patch_delete != NULL);
        patch_vt->remote_ghost_delete(pdata->user_patch);

        if(pdata->dim == 2)
        {
            FCLAW_FREE(pdata->d2);
        }
        else 
        {
            FCLAW_FREE(pdata->d3);
        }
        FCLAW_FREE(pdata);
        this_patch->user = NULL;

        fclaw_domain_data_t *ddata = glob->domain->domain_data;
        ++ddata->count_delete_patch;
    }
}

void fclaw_patch_remote_ghost_unpack(fclaw_global_t* glob,
                                       fclaw_patch_t* this_patch,
                                       int this_block_idx,
                                       int this_patch_idx,
                                       void *qdata, int time_interp)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->remote_ghost_unpack != NULL);
    patch_vt->remote_ghost_unpack(glob, this_patch, this_block_idx,
                                  this_patch_idx, qdata, time_interp);
}


/* ----------------------------------- Partitioning ----------------------------------- */

size_t fclaw_patch_partition_packsize(fclaw_global_t* glob)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->partition_packsize != NULL);

    return patch_vt->partition_packsize(glob);
}

void fclaw_patch_partition_pack(fclaw_global_t *glob,
                                  fclaw_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx,
                                  void* pack_data_here)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->partition_pack != NULL);

    patch_vt->partition_pack(glob,
                             this_patch,
                             this_block_idx,
                             this_patch_idx,
                             pack_data_here);
}


void fclaw_patch_partition_unpack(fclaw_global_t *glob,  
                                    fclaw_domain_t *new_domain,
                                    fclaw_patch_t *this_patch,
                                    int this_block_idx,
                                    int this_patch_idx,
                                    void *unpack_data_from_here)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);

    fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;

    fclaw_patch_build(glob,this_patch,this_block_idx,
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

/* We need to virtualize this because we call it from fclaw_face_neighbors */

void fclaw_patch_time_sync_f2c(fclaw_global_t* glob,
                                 fclaw_patch_t *coarse_patch,
                                 fclaw_patch_t *fine_patch,
                                 int coarse_blockno, int fine_blockno,
                                 int coarse_patchno, 
                                 int idir,
                                 int igrid,
                                 int iface_coarse,
                                 int time_interp,
                                 fclaw_patch_transform_data_t* transform_data)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->time_sync_f2c != NULL);

    patch_vt->time_sync_f2c(glob,coarse_patch,fine_patch,
                            coarse_blockno, fine_blockno,
                            coarse_patchno, 
                            idir, igrid,iface_coarse,time_interp,
                            transform_data);    
}

/* Correct for metric discontinuities at block boundaries */
void fclaw_patch_time_sync_samesize(fclaw_global_t* glob,
                                      fclaw_patch_t *this_patch,
                                      fclaw_patch_t *neighbor_patch,
                                      int iface, int idir,
                                      fclaw_patch_transform_data_t *transform_data)

{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->time_sync_samesize != NULL);

    patch_vt->time_sync_samesize(glob,this_patch,neighbor_patch,iface,idir,
                                 transform_data);    
}

void fclaw_patch_time_sync_reset(fclaw_global_t* glob,
                                   fclaw_patch_t* this_patch,
                                   int coarse_level,
                                   int reset_mode)
{
    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    FCLAW_ASSERT(patch_vt->time_sync_reset != NULL);

    patch_vt->time_sync_reset(glob,this_patch,coarse_level, reset_mode);

}

/* ----------------------------------- Virtual table ---------------------------------- */

static
fclaw_patch_vtable_t* patch_vt_new(int dim)
{
    fclaw_patch_vtable_t* vt = FCLAW_ALLOC_ZERO (fclaw_patch_vtable_t, 1);
    vt->dim = dim;
    if(dim == 2)
    {
        vt->d2 = FCLAW_ALLOC_ZERO (fclaw_patch_vtable_d2_t, 1);
    }
    else 
    {
        vt->d3 = FCLAW_ALLOC_ZERO (fclaw_patch_vtable_d3_t, 1);
    }
    return vt;
}

static
void patch_vt_destroy(void* vt)
{
    fclaw_patch_vtable_t* patch_vt = (fclaw_patch_vtable_t*) vt;
    if(patch_vt->dim == 2)
    {
        FCLAW_FREE (patch_vt->d2);
    }
    else
    {
        FCLAW_FREE (patch_vt->d3);
    }
    FCLAW_FREE (patch_vt);
}

void fclaw_patch_vtable_initialize(fclaw_global_t* glob)
{
    FCLAW_ASSERT_MESSAGE(
        glob->domain != NULL, 
        "domain needs to be stored in glob before initializing vtables"
    );

    fclaw_patch_vtable_t *patch_vt = patch_vt_new(glob->domain->dim);

    patch_vt->is_set = 1;

    FCLAW_ASSERT_MESSAGE(
        fclaw_pointer_map_get(glob->vtables,"fclaw_patch") == NULL, 
        "Vtable has already been initialized"
    );
    fclaw_pointer_map_insert(glob->vtables, "fclaw_patch", patch_vt, patch_vt_destroy);
}

/* ------------------------------ User access functions ------------------------------- */

fclaw_patch_vtable_t* fclaw_patch_vt(fclaw_global_t* glob)
{
    fclaw_patch_vtable_t* patch_vt = (fclaw_patch_vtable_t*) 
                                   fclaw_pointer_map_get(glob->vtables, "fclaw_patch");
    FCLAW_ASSERT(patch_vt != NULL);
    FCLAW_ASSERT(patch_vt->is_set != 0);
    return patch_vt;
}



/* Use this one if you have the patch or block number */
void fclaw_patch_get_info(fclaw_domain_t * domain,
                            fclaw_patch_t * patch,
                            int blockno, int patchno,
                            int *global_num, int *local_num, 
                            int *level)

{
    fclaw_block_t *block = &domain->blocks[blockno];

    /* For blockno = 0, we have 
               patchno == local_num.
       For blockno > 0, the patchno resets on each new block and we have
               patchno + block->num_patches_before = local_num.  */

    *local_num = block->num_patches_before + patchno;

    *global_num = domain->global_num_patches_before + *local_num;

    *level = patch->level;

}

/* Use this one if you don't have the patch or block number */
void fclaw_patch_get_info2(fclaw_domain_t * domain,
                            fclaw_patch_t * this_patch,
                            int *this_block_idx, int *this_patch_idx,
                            int *global_num, int *level)

{

    /* I don't completely trust this that blockno and patchno are consistent
       with p4est numbering. */
    fclaw_patch_data_t *pdata = (fclaw_patch_data_t*) this_patch->user;
    *this_patch_idx = pdata->patch_idx;
    *this_block_idx = pdata->block_idx;

    fclaw_block_t *this_block = &domain->blocks[*this_block_idx];

    *global_num = domain->global_num_patches_before +
        (this_block->num_patches_before + *this_patch_idx);

    *level = this_patch->level;

}


void*
fclaw_patch_get_user_patch(fclaw_patch_t* patch)

{
    return get_user_patch(patch);
}

fclaw_patch_data_t*
fclaw_patch_get_patch_data(fclaw_patch_t* patch)
{
    return get_patch_data(patch);
}

void* fclaw_patch_get_user_data(fclaw_global_t* glob,
                                  fclaw_patch_t* this_patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);
    return pdata->user_data;
}

void fclaw_patch_set_user_data(fclaw_global_t* glob,
                                 fclaw_patch_t* this_patch, 
                                 void* user)
{
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);
    pdata->user_data = user;
}


void* fclaw_patch_metric_patch(fclaw_global_t* glob,
                                 fclaw_patch_t *patch)
{
    fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(glob);    
    FCLAW_ASSERT(patch_vt->metric_patch != NULL);
    return patch_vt->metric_patch(patch);
}

int
fclaw_patch_get_blockno(fclaw_patch_t* this_patch)
{
    /* This is is really only here if we fail to pass in the blockno to a patch 
    routine.  */
    fclaw_patch_data_t* pdata = get_patch_data(this_patch);
    FCLAW_ASSERT(pdata->block_idx >= 0);
    return pdata->block_idx;
}

int
fclaw_patch_get_patchno(fclaw_patch_t* this_patch)
{
    /* This is is really only here if we fail to pass in the patchno to a patch 
    routine.  */
    fclaw_patch_data_t* pdata = get_patch_data(this_patch);
    FCLAW_ASSERT(pdata->patch_idx >= 0);
    return pdata->patch_idx;
}


/* -------------------------- Internal ForestClaw functions --------------------------- */

void fclaw_patch_set_face_type(fclaw_patch_t *patch,int iface,
                                 fclaw_patch_relation_t face_type)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    if(patch->dim == 2){
        pdata->d2->face_neighbors[iface] = face_type;
    } else {
        pdata->d3->face_neighbors[iface] = face_type;
    }
}

void fclaw_patch_set_corner_type(fclaw_patch_t *patch,int icorner,
                                   fclaw_patch_relation_t corner_type)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    if(patch->dim == 2){
        pdata->d2->corner_neighbors[icorner] = corner_type;
        pdata->d2->corners[icorner] = 1;
    } else {
        pdata->d3->corner_neighbors[icorner] = corner_type;
        pdata->d3->corners[icorner] = 1;
    }
}

void fclaw_patch_set_missing_corner(fclaw_patch_t *patch,int icorner)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    if(patch->dim == 2){
        pdata->d2->corners[icorner] = 0;
    } else {
        pdata->d3->corners[icorner] = 0;
    }
}

fclaw_patch_relation_t fclaw_patch_get_face_type(fclaw_patch_t* patch,
                                                     int iface)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    FCLAW_ASSERT(pdata->neighbors_set != 0);
    FCLAW_ASSERT(0 <= iface && iface < (patch->dim == 2 ? FCLAW2D_NUMFACES : FCLAW3D_NUMFACES));
    if(patch->dim == 2){
        return pdata->d2->face_neighbors[iface];
    } else {
        return pdata->d3->face_neighbors[iface];
    }
}

fclaw_patch_relation_t fclaw_patch_get_corner_type(fclaw_patch_t* patch,
                                                       int icorner)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    FCLAW_ASSERT(pdata->d2->corners[icorner] != 0);
    FCLAW_ASSERT(pdata->neighbors_set != 0);
    if(patch->dim == 2){
        return pdata->d2->corner_neighbors[icorner];
    } else {
        return pdata->d3->corner_neighbors[icorner];
    }
}

int fclaw_patch_corner_is_missing(fclaw_patch_t* patch,
                                    int icorner)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    if(patch->dim == 2){
        return !pdata->d2->corners[icorner];
    } else {
        return !pdata->d3->corners[icorner];
    }
}

void fclaw_patch_neighbors_set(fclaw_patch_t* patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    FCLAW_ASSERT(pdata->neighbors_set == 0);

    pdata->has_finegrid_neighbors = 0;
    pdata->on_coarsefine_interface = 0;
    int num_faces                          = patch->dim == 2 ? FCLAW2D_NUMFACES : FCLAW3D_NUMFACES;
    fclaw_patch_relation_t* face_neighbors = patch->dim == 2 ? pdata->d2->face_neighbors : pdata->d3->face_neighbors;
    for (int iface = 0; iface < num_faces; iface++)
    {
        fclaw_patch_relation_t nt;
        nt = face_neighbors[iface];
        if (nt == FCLAW_PATCH_HALFSIZE || (nt == FCLAW_PATCH_DOUBLESIZE))
        {
            pdata->on_coarsefine_interface = 1;
            if (nt == FCLAW_PATCH_HALFSIZE)
            {
                pdata->has_finegrid_neighbors = 1;
            }
        }
    }

    if(patch->dim == 3)
    {
        for (int iedge = 0; iedge < FCLAW3D_NUMEDGES; iedge++)
        {
            fclaw_patch_relation_t nt;
            nt = pdata->d3->edge_neighbors[iedge];
            if (nt == FCLAW_PATCH_HALFSIZE || (nt == FCLAW_PATCH_DOUBLESIZE))
            {
                pdata->on_coarsefine_interface = 1;
                if (nt == FCLAW_PATCH_HALFSIZE)
                {
                    pdata->has_finegrid_neighbors = 1;
                }
            }
        }
    }

    int num_corners;
    fclaw_patch_relation_t* corner_neighbors;
    int* corners;
    num_corners      = patch->dim == 2 ? FCLAW2D_NUMCORNERS : FCLAW3D_NUMCORNERS;
    corner_neighbors = patch->dim == 2 ? pdata->d2->corner_neighbors : pdata->d3->corner_neighbors;
    corners          = patch->dim == 2 ? pdata->d2->corners : pdata->d3->corners;
    for (int icorner = 0; icorner < num_corners; icorner++)
    {
        fclaw_patch_relation_t nt;
        int has_corner = corners[icorner];
        if (has_corner)
        {
            nt = corner_neighbors[icorner];
            if ((nt == FCLAW_PATCH_HALFSIZE) || (nt == FCLAW_PATCH_DOUBLESIZE))
            {
                pdata->on_coarsefine_interface = 1;
                if (nt == FCLAW_PATCH_HALFSIZE)
                {
                    pdata->has_finegrid_neighbors = 1;
                }
            }
        }
    }
    pdata->neighbors_set = 1;
}

void fclaw_patch_neighbors_reset(fclaw_patch_t* patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    pdata->neighbors_set = 0;
}

int fclaw_patch_neighbor_type_set(fclaw_patch_t* patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    return pdata->neighbors_set;
}


int fclaw_patch_has_finegrid_neighbors(fclaw_patch_t *patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    return pdata->has_finegrid_neighbors;
}

int fclaw_patch_on_coarsefine_interface(fclaw_patch_t *patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    return pdata->on_coarsefine_interface;
}

int
fclaw_patch_on_parallel_boundary (const fclaw_patch_t * patch)
{
    return patch->flags & FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY ? 1 : 0;
}

int* fclaw_patch_block_corner_count(fclaw_global_t* glob,
                                      fclaw_patch_t* this_patch)
{
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);
    if(this_patch->dim == 2)
    {
        return pdata->d2->block_corner_count;
    }
    else
    {
        return pdata->d3->block_corner_count;
    }
}

void fclaw_patch_set_block_corner_count(fclaw_global_t* glob,
                                          fclaw_patch_t* this_patch,
                                          int icorner, int block_corner_count)
{
    fclaw_patch_data_t *pdata = get_patch_data(this_patch);
    if(this_patch->dim == 2)
    {
        pdata->d2->block_corner_count[icorner] = block_corner_count;
    }
    else
    {
        pdata->d3->block_corner_count[icorner] = block_corner_count;
    }
}

void fclaw_patch_set_edge_type(fclaw_patch_t *patch,int iedge,
                               fclaw_patch_relation_t edge_type)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    pdata->d3->edge_neighbors[iedge] = edge_type;
}

fclaw_patch_relation_t fclaw_patch_get_edge_type(fclaw_patch_t* patch,
                                                       int iedge)
{
    fclaw_patch_data_t *pdata = get_patch_data(patch);
    FCLAW_ASSERT(pdata->neighbors_set != 0);
    FCLAW_ASSERT(0 <= iedge && iedge < FCLAW3D_NUMEDGES);
    return pdata->d3->edge_neighbors[iedge];
}
