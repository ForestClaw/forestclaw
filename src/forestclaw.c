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

#include <forestclaw.h>
#include <forestclaw2d.h>
#include <forestclaw3d.h>
#include <p4est.h>
#include <p8est.h>

fclaw_domain_t* 
fclaw_domain_wrap_2d(fclaw2d_domain_t *domain_2d)
{
    fclaw_domain_t* domain = FCLAW_ALLOC_ZERO(fclaw_domain_t,1);
    domain->dim = 2;
    domain->d2 = FCLAW_ALLOC_ZERO(fclaw_domain_d2_t,1);
    domain->d2->domain = domain_2d;
    return domain;
}
fclaw_domain_t* 
fclaw_domain_wrap_3d(fclaw3d_domain_t *domain_3d)
{
    fclaw_domain_t* domain = FCLAW_ALLOC_ZERO(fclaw_domain_t,1);
    domain->dim = 3;
    domain->d3 = FCLAW_ALLOC_ZERO(fclaw_domain_d3_t,1);
    domain->d3->domain = domain_3d;
    return domain;
}

int fclaw_patch_edge_neighbors (fclaw_domain_t * domain,
                                int blockno, int patchno, int edgeno,
                                int *rproc, int *rblockno, int *rpatchno,
                                int *redge,
                                fclaw_patch_relation_t * neighbor_size)
{
    return fclaw3d_patch_edge_neighbors(domain->d3->domain,blockno,patchno,edgeno,
                                        rproc,rblockno,rpatchno,redge,
                                        (fclaw3d_patch_relation_t*) neighbor_size);
}
void
fclaw_domain_attribute_add (fclaw_domain_t * domain,
                              const char *name, void *attribute)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (!sc_keyvalue_exists (a, name));
    sc_keyvalue_set_pointer (a, name, attribute);
}

void *
fclaw_domain_attribute_access (fclaw_domain_t * domain,
                                 const char *name, void *default_attr)
{
    sc_keyvalue_t *a = domain->attributes;

    FCLAW_ASSERT (a != NULL);
    return sc_keyvalue_get_pointer (a, name, default_attr);
}

void
fclaw_domain_attribute_remove (fclaw_domain_t * domain, const char *name)
{
    sc_keyvalue_t *a = domain->attributes;
#ifdef FCLAW_ENABLE_DEBUG
    sc_keyvalue_entry_type_t et;
#endif

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (sc_keyvalue_exists (a, name));
#ifndef FCLAW_ENABLE_DEBUG
    (void)
#else
    et =
#endif
        sc_keyvalue_unset (a, name);
    FCLAW_ASSERT (et == SC_KEYVALUE_ENTRY_POINTER);
}

int
fclaw_domain_dimension (const fclaw_domain_t * domain)
{
    return domain->dim;
}

int
fclaw_domain_refine_factor (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_CHILDREN/2 : P8EST_CHILDREN/2;
}

int
fclaw_domain_num_children (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_faces (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_FACES : P8EST_FACES;
}

int
fclaw_domain_num_edges (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? 0 : P8EST_EDGES;
}

int
fclaw_domain_num_corners (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_face_corners (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? P4EST_HALF : P8EST_HALF;
}

int
fclaw_domain_num_orientations (const fclaw_domain_t * domain)
{
    return (domain->dim == 2) ? (P4EST_FACES * P4EST_HALF) : (P8EST_FACES * P8EST_HALF);
}

void
fclaw_domain_edge_faces (const fclaw_domain_t * domain,
                         int iedge, int faces[2])
{
    FCLAW_ASSERT (0 <= iedge && iedge < fclaw_domain_num_edges(domain));
    faces[0] = p8est_edge_faces[iedge][0];
    faces[1] = p8est_edge_faces[iedge][1];
}

void
fclaw_domain_corner_faces (const fclaw_domain_t * domain,
                             int icorner, int faces[3])
{
    if(domain->dim == 2)
    {
        return fclaw2d_domain_corner_faces(domain->d2->domain, icorner, faces);
    }
    else if(domain->dim == 3)
    {
        return fclaw3d_domain_corner_faces(domain->d3->domain, icorner, faces);
    }
    else 
    {
        SC_ABORT_NOT_REACHED();
    }
}

int
fclaw_patch_corner_dimension (const fclaw_patch_t * patch, int cornerno)
{
    if(patch->dim == 2)
    {
        return fclaw2d_patch_corner_dimension(patch->d2, cornerno);
    }
    else if(patch->dim == 3)
    {
        return fclaw3d_patch_corner_dimension(patch->d3, cornerno);
    }
    else 
    {
        SC_ABORT_NOT_REACHED();
    }
}

int
fclaw_patch_childid (const fclaw_patch_t * patch)
{
    if(patch->dim == 2)
    {
        return fclaw2d_patch_childid(patch->d2);
    } 
    else if(patch->dim == 3)
    {
        return fclaw3d_patch_childid(patch->d3);
    }
    else 
    {
        SC_ABORT_NOT_REACHED();
    }
}

int
fclaw_patch_is_first_sibling (const fclaw_patch_t * patch)
{
    if(patch->dim == 2)
    {
        return fclaw2d_patch_is_first_sibling(patch->d2);
    }
    else if(patch->dim == 3)
    {
        return fclaw3d_patch_is_first_sibling(patch->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

int
fclaw_patch_is_ghost (const fclaw_patch_t * patch)
{
    if(patch->dim == 2)
    {
        return fclaw2d_patch_is_ghost(patch->d2);
    }
    else if(patch->dim == 3)
    {
        return fclaw3d_patch_is_ghost(patch->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

int 
fclaw_patch_boundary_type(fclaw_domain_t *domain, int blockno, int patchno, int boundaries[6])
{
    if(domain->dim == 2)
    {
        return fclaw2d_patch_boundary_type(domain->d2->domain,blockno,patchno,boundaries);
    }
    else if (domain->dim == 3)
    {
        return fclaw3d_patch_boundary_type(domain->d3->domain,blockno,patchno,boundaries);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

int 
fclaw_patch_normal_match(fclaw_domain_t *domain, int blockno, int patchno, int faceno)
{
    if(domain->dim == 2)
    {
        return fclaw2d_patch_normal_match(domain->d2->domain,blockno,patchno,faceno);
    }
    else if (domain->dim == 3)
    {
        return fclaw3d_patch_normal_match(domain->d3->domain,blockno,patchno,faceno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

fclaw_patch_relation_t 
fclaw_patch_face_neighbors(fclaw_domain_t *domain, int blockno, int patchno, int faceno, int rproc[2], int *rblockno, int rpatchno[2], int *rfaceno)
{
    if(domain->dim == 2)
    {
        return (fclaw_patch_relation_t) fclaw2d_patch_face_neighbors(domain->d2->domain,blockno,patchno,faceno,rproc,rblockno,rpatchno,rfaceno);
    }
    else if (domain->dim == 3)
    {
        return (fclaw_patch_relation_t) fclaw3d_patch_face_neighbors(domain->d3->domain,blockno,patchno,faceno,rproc,rblockno,rpatchno,rfaceno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void 
fclaw_patch_face_swap (int dim, int *faceno, int *rfaceno)
{
    if(dim == 2)
    {
        fclaw2d_patch_face_swap(faceno,rfaceno);
    }
    else if (dim == 3)
    {
        fclaw3d_patch_face_swap(faceno,rfaceno);
    }
    else
    {
        fclaw_abortf("fclaw_patch_face_swap : dim = %d not supported\n",dim);
    }
}

void 
fclaw_patch_face_transformation (int dim, int faceno, int rfaceno,
                                 int ftransform[])
{
    if(dim == 2)
    {
        fclaw2d_patch_face_transformation(faceno,rfaceno,ftransform);
    }
    else if(dim == 3)
    {
        fclaw3d_patch_face_transformation(faceno,rfaceno,ftransform);
    }
    else
    {
        fclaw_abortf("fclaw_patch_face_transformation : dim = %d not supported\n",dim);
    }
}

void 
fclaw_patch_face_transformation_block (int dim, int ftransform[],
                                       int sameblock)
{
    if(dim == 2)
    {
        fclaw2d_patch_face_transformation_block(ftransform, sameblock);
    }
    else if(dim == 3)
    {
        fclaw3d_patch_face_transformation_block(ftransform, sameblock);
    }
    else
    {
        fclaw_abortf("fclaw_patch_face_transformation_block : dim = %d not supported\n",dim);
    }
}

void 
fclaw_patch_face_transformation_intra (int dim, int ftransform[])
{
    if(dim == 2)
    {
        fclaw2d_patch_face_transformation_intra(ftransform);
    }
    else if(dim == 3)
    {
        fclaw3d_patch_face_transformation_intra(ftransform);
    }
    else
    {
        fclaw_abortf("fclaw_patch_face_transformation_intra : dim = %d not supported\n",dim);
    }
}

int 
fclaw_patch_face_transformation_valid (int dim, const int ftransform[])
{
    if(dim == 2)
    {
        return fclaw2d_patch_face_transformation_valid(ftransform);
    }
    else if(dim == 3)
    {
        return fclaw3d_patch_face_transformation_valid(ftransform);
    }
    else
    {
        fclaw_abortf("fclaw_patch_face_transformation_valid : dim = %d not supported\n",dim);
    }
}

void 
fclaw_patch_transform_face_2d (fclaw_patch_t * ipatch,
                               fclaw_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my, int based, int *i, int *j)
{
    fclaw2d_patch_transform_face(ipatch->d2,opatch->d2,ftransform,mx,my,based,i,j);
}

void 
fclaw_patch_transform_face2_2d (fclaw_patch_t * ipatch,
                                fclaw_patch_t * opatch,
                                const int ftransform[],
                                int mx, int my, int based, int i[],
                                int j[])
{
    fclaw2d_patch_transform_face2(ipatch->d2,opatch->d2,ftransform,mx,my,based,i,j);
}

void 
fclaw_patch_transform_face_3d (fclaw_patch_t * ipatch,
                               fclaw_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my, int mz, int based,
                               int *i, int *j, int *k)
{
    fclaw3d_patch_transform_face(ipatch->d3,opatch->d3,ftransform,mx,my,mz,based,i,j,k);
}

void 
fclaw_patch_transform_face2_3d (fclaw_patch_t * ipatch,
                                fclaw_patch_t * opatch,
                                const int ftransform[],
                                int mx, int my, int mz, int based,
                                int i[], int j[], int k[])
{
    fclaw3d_patch_transform_face2(ipatch->d3,opatch->d3,ftransform,mx,my,mz,based,i,j,k);
}


int 
fclaw_patch_corner_neighbors(fclaw_domain_t *domain, int blockno, int patchno, int cornerno, int *rproc, int *rblockno, int *rpatchno, int *rcorner, fclaw_patch_relation_t *neighbor_size)
{
    if(domain->dim == 2)
    {
        fclaw2d_patch_relation_t neighbor_size2d;
        int retval = fclaw2d_patch_corner_neighbors(domain->d2->domain,blockno,patchno,cornerno,rproc,rblockno,rpatchno,rcorner,&neighbor_size2d);
        *neighbor_size = (fclaw_patch_relation_t) neighbor_size2d;
        return retval;
    }
    else if (domain->dim == 3)
    {
        fclaw3d_patch_relation_t neighbor_size3d;
        int retval = fclaw3d_patch_corner_neighbors(domain->d3->domain,blockno,patchno,cornerno,rproc,rblockno,rpatchno,rcorner,&neighbor_size3d);
        *neighbor_size = (fclaw_patch_relation_t) neighbor_size3d;
        return retval;
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_patch_corner_swap (int dim, int *cornerno, int *rcornerno)
{
    if(dim == 2)
    {
        fclaw2d_patch_corner_swap(cornerno, rcornerno);
    }
    else if(dim == 3)
    {
        fclaw3d_patch_corner_swap(cornerno, rcornerno);
    }
    else
    {
        fclaw_abortf("fclaw_patch_corner_swap : dim = %d not supported\n",dim);
    }
}

void fclaw_patch_transform_corner_2d (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my,
                                      int based, int *i, int *j)
{
    fclaw2d_patch_transform_corner(ipatch->d2,opatch->d2,icorner,is_block_boundary,
                                   mx,my,based,i,j);
}

void fclaw_patch_transform_corner2_2d (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int based,
                                       int i[], int j[])
{
    fclaw2d_patch_transform_corner2(ipatch->d2,opatch->d2,icorner,is_block_boundary,
                                    mx,my,based,i,j);
}

/** Transform a patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the two patches are of the SAME size and that the
 * patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] icorner      Corner number of this patch to transform across.
 *                          This function assumes ocorner == icorner ^ 7, so
 *                          ocorner is the opposite corner of icorner.
 * \param [in] is_block_boundary      Set to true for a block corner.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] mz           Number of cells along z direction of patch.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        Integer coordinate along x-axis in \a based .. \a mx.
 * \param [in,out] j        Integer coordinate along y-axis in \a based .. \a my.
 * \param [in,out] k        Integer coordinate along z-axis in \a based .. \a mz.
 */
void fclaw_patch_transform_corner_3d (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my, int mz,
                                      int based, int *i, int *j, int *k)
{
    fclaw3d_patch_transform_corner(ipatch->d3,opatch->d3,icorner,is_block_boundary,
                                   mx,my,mz,based,i,j,k);
}

void fclaw_patch_transform_corner2_3d (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int mz, int based,
                                       int i[], int j[], int k[])
{
    fclaw3d_patch_transform_corner2(ipatch->d3,opatch->d3,icorner,is_block_boundary,
                                    mx,my,mz,based,i,j,k);
}

void
fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                              fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    FCLAW_ASSERT (0 <= level && level <= domain->possible_maxlevel);
    for (i = 0; i < domain->num_blocks; ++i)
    {
        block = domain->blocks + i;
        for (patch = block->patchbylevel[level];
             patch != NULL; patch = patch->u.next)
        {
            j = (int) (patch - block->patches);
            FCLAW_ASSERT (0 <= j && j < block->num_patches);
            FCLAW_ASSERT (patch->level == level);
            pcb (domain, patch, i, j, user);
        }
    }
}

void
fclaw_domain_iterate_patches (fclaw_domain_t * domain,
                                fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = domain->blocks + i;
        for (j = 0; j < block->num_patches; j++)
        {
            patch = block->patches + j;
            pcb (domain, patch, i, j, user);
        }
    }
}

void
fclaw_domain_iterate_families (fclaw_domain_t * domain,
                                 fclaw_patch_callback_t pcb, void *user)
{
    int i, j;
    fclaw_block_t *block;
    fclaw_patch_t *patch;

    for (i = 0; i < domain->num_blocks; i++)
    {
        block = domain->blocks + i;
        for (j = 0; j < block->num_patches; j++)
        {
            patch = block->patches + j;
            if (fclaw_patch_is_first_sibling (patch))
            {
#ifdef FCLAW_ENABLE_DEBUG
                int k;
                for (k = 0; k < fclaw_domain_num_children(domain); ++k)
                {
                    FCLAW_ASSERT (j + k < block->num_patches);
                    FCLAW_ASSERT (fclaw_patch_childid (patch + k) == k);
                }
#endif
                pcb (domain, patch, i, j, user);
                j += fclaw_domain_num_children(domain) - 1;
            }
        }
    }
}





#define FCLAW2D_DOMAIN_TAG_SERIALIZE 4526
#define FCLAW3D_DOMAIN_TAG_SERIALIZE 4527

int get_tag(int dim)
{
    if(dim == 2)
    {
        return FCLAW2D_DOMAIN_TAG_SERIALIZE;
    }
    else
    {
        return FCLAW3D_DOMAIN_TAG_SERIALIZE;
    }
}

double
fclaw_domain_global_maximum (fclaw_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = sc_MPI_Allreduce (&d, &gd, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

double
fclaw_domain_global_sum (fclaw_domain_t * domain, double d)
{
    int mpiret;
    double gd;

    mpiret = sc_MPI_Allreduce (&d, &gd, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return gd;
}

void
fclaw_domain_barrier (fclaw_domain_t * domain)
{
    int mpiret;

    mpiret = sc_MPI_Barrier (domain->mpicomm);
    SC_CHECK_MPI (mpiret);
}

void
fclaw_domain_serialization_enter (fclaw_domain_t * domain)
{
    int mpiret;
    int i;
    sc_MPI_Status status;

    if (domain->mpirank > 0)
    {
        mpiret = sc_MPI_Recv (&i, 1, sc_MPI_INT, domain->mpirank - 1,
                              get_tag(domain->dim), domain->mpicomm,
                              &status);
        SC_CHECK_MPI (mpiret);
        FCLAW_ASSERT (i == 0);
    }
}

void
fclaw_domain_serialization_leave (fclaw_domain_t * domain)
{
    int mpiret;
    int i = 0;

    if (domain->mpirank + 1 < domain->mpisize)
    {
        mpiret = sc_MPI_Send (&i, 1, sc_MPI_INT, domain->mpirank + 1,
                              get_tag(domain->dim), domain->mpicomm);
        SC_CHECK_MPI (mpiret);
    }
}

void
fclaw_domain_set_refinement (fclaw_domain_t * domain,
                             int smooth_refine, int smooth_level,
                             int coarsen_delay)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_set_refinement(domain->d2->domain,smooth_refine,smooth_level,
                                      coarsen_delay);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_set_refinement(domain->d3->domain,smooth_refine,smooth_level,
                                      coarsen_delay);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_patch_mark_refine(fclaw_domain_t *domain, int blockno, int patchno)
{
    if(domain->dim == 2)
    {
        fclaw2d_patch_mark_refine(domain->d2->domain,blockno,patchno);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_patch_mark_refine(domain->d3->domain,blockno,patchno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_patch_mark_coarsen(fclaw_domain_t *domain, int blockno, int patchno)
{
    if(domain->dim == 2)
    {
        fclaw2d_patch_mark_coarsen(domain->d2->domain,blockno,patchno);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_patch_mark_coarsen(domain->d3->domain,blockno,patchno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

typedef struct mcb_wrap_user
{
    fclaw_match_callback_t mcb;
    fclaw_domain_t* old_domain;
    fclaw_domain_t* new_domain;
    void *user;
} mcb_wrap_user_t;

static 
void
mcb_wrap2d(fclaw2d_domain_t * old_domain,
           fclaw2d_patch_t * old_patch_2d,
           fclaw2d_domain_t * new_domain,
           fclaw2d_patch_t * new_patch_2d,
           fclaw2d_patch_relation_t newsize,
           int blockno,
           int old_patchno, int new_patchno,
           void *user)
{
    mcb_wrap_user_t* wrap_user = (mcb_wrap_user_t*) user;

    fclaw_patch_t* old_patch = wrap_user->old_domain->blocks[blockno].patches + old_patchno;
    fclaw_patch_t* new_patch = wrap_user->new_domain->blocks[blockno].patches + new_patchno;

    FCLAW_ASSERT(old_patch->d2 == old_patch_2d);
    FCLAW_ASSERT(new_patch->d2 == new_patch_2d);

    wrap_user->mcb(wrap_user->old_domain, old_patch,
                   wrap_user->new_domain, new_patch,
                   (fclaw_patch_relation_t) newsize,
                   blockno, old_patchno, new_patchno,
                   wrap_user->user);
}

static 
void
mcb_wrap3d(fclaw3d_domain_t * old_domain,
           fclaw3d_patch_t * old_patch_3d,
           fclaw3d_domain_t * new_domain,
           fclaw3d_patch_t * new_patch_3d,
           fclaw3d_patch_relation_t newsize,
           int blockno,
           int old_patchno, int new_patchno,
           void *user)
{
    mcb_wrap_user_t* wrap_user = (mcb_wrap_user_t*) user;

    fclaw_patch_t* old_patch = wrap_user->old_domain->blocks[blockno].patches + old_patchno;
    fclaw_patch_t* new_patch = wrap_user->new_domain->blocks[blockno].patches + new_patchno;

    FCLAW_ASSERT(old_patch->d3 == old_patch_3d);
    FCLAW_ASSERT(new_patch->d3 == new_patch_3d);

    wrap_user->mcb(wrap_user->old_domain, old_patch,
                   wrap_user->new_domain, new_patch,
                   (fclaw_patch_relation_t) newsize,
                   blockno, old_patchno, new_patchno,
                   wrap_user->user);
}

void
fclaw_domain_iterate_adapted(fclaw_domain_t *old_domain, fclaw_domain_t *new_domain, fclaw_match_callback_t mcb, void *user)
{
    mcb_wrap_user_t mcb_wrap_user;
    mcb_wrap_user.mcb = mcb;
    mcb_wrap_user.old_domain = old_domain;
    mcb_wrap_user.new_domain = new_domain;
    mcb_wrap_user.user = user;

    if(old_domain->dim == 2)
    {
        fclaw2d_domain_iterate_adapted(old_domain->d2->domain,
                                       new_domain->d2->domain,
                                       mcb_wrap2d,
                                       &mcb_wrap_user);
    }
    else if (old_domain->dim == 3)
    {
        fclaw3d_domain_iterate_adapted(old_domain->d3->domain,
                                       new_domain->d3->domain,
                                       mcb_wrap3d,
                                       &mcb_wrap_user);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_allocate_before_partition(fclaw_domain_t *domain, size_t data_size, void ***patch_data)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_allocate_before_partition(domain->d2->domain,data_size,patch_data);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_allocate_before_partition(domain->d3->domain,data_size,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_retrieve_after_partition(fclaw_domain_t *domain, void ***patch_data)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_retrieve_after_partition(domain->d2->domain,patch_data);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_retrieve_after_partition(domain->d3->domain,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

typedef struct tcb_wrap_user
{
    fclaw_transfer_callback_t tcb;
    fclaw_domain_t* old_domain;
    fclaw_domain_t* new_domain;
    void *user;
} tcb_wrap_user_t;

static 
void
tcb_wrap2d(fclaw2d_domain_t * old_domain,
           fclaw2d_patch_t * old_patch_2d,
           fclaw2d_domain_t * new_domain,
           fclaw2d_patch_t * new_patch_2d,
           int blockno,
           int old_patchno, int new_patchno,
           void *user)
{
    tcb_wrap_user_t* wrap_user = (tcb_wrap_user_t*) user;

    fclaw_patch_t* old_patch = wrap_user->old_domain->blocks[blockno].patches + old_patchno;
    fclaw_patch_t* new_patch = wrap_user->new_domain->blocks[blockno].patches + new_patchno;

    FCLAW_ASSERT(old_patch->d2 == old_patch_2d);
    FCLAW_ASSERT(new_patch->d2 == new_patch_2d);

    wrap_user->tcb(wrap_user->old_domain, old_patch,
                   wrap_user->new_domain, new_patch,
                   blockno, old_patchno, new_patchno,
                   wrap_user->user);
}

static 
void
tcb_wrap3d(fclaw3d_domain_t * old_domain,
           fclaw3d_patch_t * old_patch_3d,
           fclaw3d_domain_t * new_domain,
           fclaw3d_patch_t * new_patch_3d,
           int blockno,
           int old_patchno, int new_patchno,
           void *user)
{
    tcb_wrap_user_t* wrap_user = (tcb_wrap_user_t*) user;

    fclaw_patch_t* old_patch = wrap_user->old_domain->blocks[blockno].patches + old_patchno;
    fclaw_patch_t* new_patch = wrap_user->new_domain->blocks[blockno].patches + new_patchno;

    FCLAW_ASSERT(old_patch->d3 == old_patch_3d);
    FCLAW_ASSERT(new_patch->d3 == new_patch_3d);

    wrap_user->tcb(wrap_user->old_domain, old_patch,
                   wrap_user->new_domain, new_patch,
                   blockno, old_patchno, new_patchno,
                   wrap_user->user);
}

void fclaw_domain_iterate_partitioned(fclaw_domain_t *old_domain, fclaw_domain_t *new_domain, fclaw_transfer_callback_t tcb, void *user)
{
    tcb_wrap_user_t tcb_wrap_user;
    tcb_wrap_user.tcb = tcb;
    tcb_wrap_user.old_domain = old_domain;
    tcb_wrap_user.new_domain = new_domain;
    tcb_wrap_user.user = user;

    if(old_domain->dim == 2)
    {
        fclaw2d_domain_iterate_partitioned(old_domain->d2->domain,
                                           new_domain->d2->domain,
                                           tcb_wrap2d,&tcb_wrap_user);
    }
    else if (old_domain->dim == 3)
    {
        fclaw3d_domain_iterate_partitioned(old_domain->d3->domain,
                                           new_domain->d3->domain,
                                            tcb_wrap3d,&tcb_wrap_user);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_free_after_partition(fclaw_domain_t *domain, void ***patch_data)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_free_after_partition(domain->d2->domain,patch_data);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_free_after_partition(domain->d3->domain,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_allocate_before_exchange (fclaw_domain_t * domain,
                                       size_t data_size)
{
    if(domain->dim == 2)
    {
        domain->d2->exchange = fclaw2d_domain_allocate_before_exchange(domain->d2->domain,data_size);
    }
    else if (domain->dim == 3)
    {
        domain->d3->exchange = fclaw3d_domain_allocate_before_exchange(domain->d3->domain,data_size);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_ghost_exchange (fclaw_domain_t * domain,
                                  int exchange_minlevel,
                                  int exchange_maxlevel)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_ghost_exchange(domain->d2->domain,
                                      domain->d2->exchange,
                                      exchange_minlevel,
                                      exchange_maxlevel);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_ghost_exchange(domain->d3->domain,
                                      domain->d3->exchange,
                                      exchange_minlevel,
                                      exchange_maxlevel);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_ghost_exchange_begin (fclaw_domain_t * domain,
                                        int exchange_minlevel,
                                        int exchange_maxlevel)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_ghost_exchange_begin(domain->d2->domain,
                                            domain->d2->exchange,
                                            exchange_minlevel,
                                            exchange_maxlevel);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_ghost_exchange_begin(domain->d3->domain,
                                            domain->d3->exchange,
                                            exchange_minlevel,
                                            exchange_maxlevel);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_ghost_exchange_end (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_ghost_exchange_end(domain->d2->domain,
                                          domain->d2->exchange);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_ghost_exchange_end(domain->d3->domain,
                                          domain->d3->exchange);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_free_after_exchange (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_free_after_exchange(domain->d2->domain,
                                           domain->d2->exchange);
        domain->d2->exchange = NULL;
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_free_after_exchange(domain->d3->domain,
                                           domain->d3->exchange);
        domain->d3->exchange = NULL;
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_indirect_begin (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        domain->d2->indirect = fclaw2d_domain_indirect_begin(domain->d2->domain);
    }
    else if (domain->dim == 3)
    {
        domain->d3->indirect = fclaw3d_domain_indirect_begin(domain->d3->domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_indirect_end (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_indirect_end(domain->d2->domain,domain->d2->indirect);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_indirect_end(domain->d3->domain,domain->d3->indirect);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

fclaw_patch_relation_t
fclaw_domain_indirect_neighbors (fclaw_domain_t * domain,
                                 int ghostno, int faceno, int rproc[],
                                 int *rblockno, int rpatchno[],
                                 int *rfaceno)
{
    if(domain->dim == 2)
    {
        return (fclaw_patch_relation_t) fclaw2d_domain_indirect_neighbors(domain->d2->domain,domain->d2->indirect,
                                                                          ghostno,faceno,rproc,rblockno,
                                                                          rpatchno,rfaceno);
    }
    else if (domain->dim == 3)
    {
        return (fclaw_patch_relation_t) fclaw3d_domain_indirect_neighbors(domain->d3->domain,domain->d3->indirect,
                                                                          ghostno,faceno,rproc,rblockno,
                                                                          rpatchno,rfaceno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_indirect_destroy (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_indirect_destroy(domain->d2->domain,domain->d2->indirect);
        domain->d2->indirect = NULL;
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_indirect_destroy(domain->d3->domain,domain->d3->indirect);
        domain->d3->indirect = NULL;
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}


int fclaw_domain_is_meta (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        return fclaw2d_domain_is_meta(domain->d2->domain);
    }
    else if (domain->dim == 3)
    {
        return fclaw3d_domain_is_meta(domain->d3->domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_init_meta (fclaw_domain_t *domain, int mpirank)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_init_meta(domain->d2->domain,mpirank);
    }
    else if (domain->dim == 3)
    {
        fclaw3d_domain_init_meta(domain->d3->domain,mpirank);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}


