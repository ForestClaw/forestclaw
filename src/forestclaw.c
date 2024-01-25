/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
#include <fclaw2d_wrap.h>
#include <fclaw3d_wrap.h>
#include <forestclaw2d.h>
#include <forestclaw3d.h>
#include <p4est.h>
#include <p8est.h>

int fclaw_patch_edge_neighbors (fclaw_domain_t * domain,
                                int blockno, int patchno, int edgeno,
                                int *rproc, int *rblockno, int *rpatchno,
                                int *redge,
                                fclaw_patch_relation_t * neighbor_size)
{
    FCLAW_ASSERT(domain->refine_dim == 3);
    return fclaw3d_patch_edge_neighbors(domain->d3,blockno,patchno,edgeno,
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
    return domain->refine_dim;
}

int
fclaw_domain_refine_factor (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? P4EST_CHILDREN/2 : P8EST_CHILDREN/2;
}

int
fclaw_domain_num_siblings (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_faces (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? P4EST_FACES : P8EST_FACES;
}

int
fclaw_domain_num_edges (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? 0 : P8EST_EDGES;
}

int
fclaw_domain_num_corners (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? P4EST_CHILDREN : P8EST_CHILDREN;
}

int
fclaw_domain_num_face_corners (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? P4EST_HALF : P8EST_HALF;
}

int
fclaw_domain_num_orientations (const fclaw_domain_t * domain)
{
    return (domain->refine_dim == 2) ? (P4EST_FACES * P4EST_HALF) : (P8EST_FACES * P8EST_HALF);
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
    if(domain->refine_dim == 2)
    {
        return fclaw2d_domain_corner_faces(domain->d2, icorner, faces);
    }
    else if(domain->refine_dim == 3)
    {
        return fclaw3d_domain_corner_faces(domain->d3, icorner, faces);
    }
    else 
    {
        SC_ABORT_NOT_REACHED();
    }
}

int
fclaw_patch_corner_dimension (const fclaw_patch_t * patch, int cornerno)
{
    if(patch->refine_dim == 2)
    {
        return fclaw2d_patch_corner_dimension(patch->d2, cornerno);
    }
    else if(patch->refine_dim == 3)
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
    if(patch->refine_dim == 2)
    {
        return fclaw2d_patch_childid(patch->d2);
    } 
    else if(patch->refine_dim == 3)
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
    if(patch->refine_dim == 2)
    {
        return fclaw2d_patch_is_first_sibling(patch->d2);
    } 
    else if(patch->refine_dim == 3)
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
    if(patch->refine_dim == 2)
    {
        return fclaw2d_patch_is_ghost(patch->d2);
    }
    else if(patch->refine_dim == 3)
    {
        return fclaw3d_patch_is_ghost(patch->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

int 
fclaw_patch_boundary_type(fclaw_domain_t *domain, int blockno, int patchno, int boundaries[])
{
    if(domain->refine_dim == 2)
    {
        return fclaw2d_patch_boundary_type(domain->d2,blockno,patchno,boundaries);
    }
    else if (domain->refine_dim == 3)
    {
        return fclaw3d_patch_boundary_type(domain->d3,blockno,patchno,boundaries);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

int 
fclaw_patch_normal_match(fclaw_domain_t *domain, int blockno, int patchno, int faceno)
{
    if(domain->refine_dim == 2)
    {
        return fclaw2d_patch_normal_match(domain->d2,blockno,patchno,faceno);
    }
    else if (domain->refine_dim == 3)
    {
        return fclaw3d_patch_normal_match(domain->d3,blockno,patchno,faceno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

fclaw_patch_relation_t 
fclaw_patch_face_neighbors(fclaw_domain_t *domain, int blockno, int patchno, int faceno, int rproc[], int *rblockno, int rpatchno[], int *rfaceno)
{
    if(domain->refine_dim == 2)
    {
        return (fclaw_patch_relation_t) fclaw2d_patch_face_neighbors(domain->d2,blockno,patchno,faceno,rproc,rblockno,rpatchno,rfaceno);
    }
    else if (domain->refine_dim == 3)
    {
        return (fclaw_patch_relation_t) fclaw3d_patch_face_neighbors(domain->d3,blockno,patchno,faceno,rproc,rblockno,rpatchno,rfaceno);
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
fclaw_patch_2d_transform_face (fclaw_patch_t * ipatch,
                               fclaw_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my, int based, int *i, int *j)
{
    FCLAW_ASSERT(ipatch->refine_dim == 2);
    FCLAW_ASSERT(opatch->refine_dim == 2);
    fclaw2d_patch_transform_face(ipatch->d2,opatch->d2,ftransform,mx,my,based,i,j);
}

void 
fclaw_patch_2d_transform_face2 (fclaw_patch_t * ipatch,
                                fclaw_patch_t * opatch,
                                const int ftransform[],
                                int mx, int my, int based, int i[],
                                int j[])
{
    FCLAW_ASSERT(ipatch->refine_dim == 2);
    FCLAW_ASSERT(opatch->refine_dim == 2);
    fclaw2d_patch_transform_face2(ipatch->d2,opatch->d2,ftransform,mx,my,based,i,j);
}

void 
fclaw_patch_3d_transform_face (fclaw_patch_t * ipatch,
                               fclaw_patch_t * opatch,
                               const int ftransform[],
                               int mx, int my, int mz, int based,
                               int *i, int *j, int *k)
{
    FCLAW_ASSERT(ipatch->refine_dim == 3);
    FCLAW_ASSERT(opatch->refine_dim == 3);
    fclaw3d_patch_transform_face(ipatch->d3,opatch->d3,ftransform,mx,my,mz,based,i,j,k);
}

void 
fclaw_patch_3d_transform_face2 (fclaw_patch_t * ipatch,
                                fclaw_patch_t * opatch,
                                const int ftransform[],
                                int mx, int my, int mz, int based,
                                int i[], int j[], int k[])
{
    FCLAW_ASSERT(ipatch->refine_dim == 3);
    FCLAW_ASSERT(opatch->refine_dim == 3);
    fclaw3d_patch_transform_face2(ipatch->d3,opatch->d3,ftransform,mx,my,mz,based,i,j,k);
}


void fclaw_patch_3d_transform_edge (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int iedge, int is_block_boundary,
                                      int mx, int my, int mz,
                                      int based, int *i, int *j, int *k)
{
    FCLAW_ASSERT(ipatch->refine_dim == 3);
    FCLAW_ASSERT(opatch->refine_dim == 3);
    fclaw3d_patch_transform_edge(ipatch->d3,opatch->d3,iedge,is_block_boundary,
                                   mx,my,mz,based,i,j,k);
}

void fclaw_patch_3d_transform_edge2 (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int iedge, int is_block_boundary,
                                       int mx, int my, int mz, int based,
                                       int i[], int j[], int k[])
{
    fclaw3d_patch_transform_edge2(ipatch->d3,opatch->d3,iedge,is_block_boundary,
                                    mx,my,mz,based,i,j,k);
}

int 
fclaw_patch_corner_neighbors(fclaw_domain_t *domain, int blockno, int patchno, int cornerno, int *rproc, int *rblockno, int *rpatchno, int *rcorner, fclaw_patch_relation_t *neighbor_size)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_patch_relation_t neighbor_size2d;
        int retval = fclaw2d_patch_corner_neighbors(domain->d2,blockno,patchno,cornerno,rproc,rblockno,rpatchno,rcorner,&neighbor_size2d);
        *neighbor_size = (fclaw_patch_relation_t) neighbor_size2d;
        return retval;
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_patch_relation_t neighbor_size3d;
        int retval = fclaw3d_patch_corner_neighbors(domain->d3,blockno,patchno,cornerno,rproc,rblockno,rpatchno,rcorner,&neighbor_size3d);
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

void fclaw_patch_2d_transform_corner (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my,
                                      int based, int *i, int *j)
{
    FCLAW_ASSERT(ipatch->refine_dim == 2);
    FCLAW_ASSERT(opatch->refine_dim == 2);
    fclaw2d_patch_transform_corner(ipatch->d2,opatch->d2,icorner,is_block_boundary,
                                   mx,my,based,i,j);
}

void fclaw_patch_2d_transform_corner2 (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int based,
                                       int i[], int j[])
{
    FCLAW_ASSERT(ipatch->refine_dim == 2);
    FCLAW_ASSERT(opatch->refine_dim == 2);
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
void fclaw_patch_3d_transform_corner (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my, int mz,
                                      int based, int *i, int *j, int *k)
{
    FCLAW_ASSERT(ipatch->refine_dim == 3);
    FCLAW_ASSERT(opatch->refine_dim == 3);
    fclaw3d_patch_transform_corner(ipatch->d3,opatch->d3,icorner,is_block_boundary,
                                   mx,my,mz,based,i,j,k);
}

void fclaw_patch_3d_transform_corner2 (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int mz, int based,
                                       int i[], int j[], int k[])
{
    FCLAW_ASSERT(ipatch->refine_dim == 3);
    FCLAW_ASSERT(opatch->refine_dim == 3);
    fclaw3d_patch_transform_corner2(ipatch->d3,opatch->d3,icorner,is_block_boundary,
                                    mx,my,mz,based,i,j,k);
}

void
fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                              fclaw_patch_callback_t pcb, void *user)
{
    fclaw_patch_wrap_user_t wrap;
    wrap.pcb = pcb;
    wrap.user = user;

    if (domain->refine_dim == 2)
    {
        fclaw2d_domain_iterate_level(domain->d2,level,fclaw2d_patch_callback_wrap,
                                     &wrap);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_iterate_level(domain->d3,level,fclaw3d_patch_callback_wrap,
                                     &wrap);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_iterate_patches (fclaw_domain_t * domain,
                                fclaw_patch_callback_t pcb, void *user)
{
    fclaw_patch_wrap_user_t wrap;
    wrap.pcb = pcb;
    wrap.user = user;

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_iterate_patches(domain->d2,fclaw2d_patch_callback_wrap,
                                       &wrap);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_iterate_patches(domain->d3,fclaw3d_patch_callback_wrap,
                                       &wrap);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_iterate_families (fclaw_domain_t * domain,
                                 fclaw_patch_callback_t pcb, void *user)
{
    fclaw_patch_wrap_user_t wrap;
    wrap.pcb = pcb;
    wrap.user = user;

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_iterate_families(domain->d2,fclaw2d_patch_callback_wrap,
                                        &wrap);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_iterate_families(domain->d3,fclaw3d_patch_callback_wrap,
                                        &wrap);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
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
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_serialization_enter(domain->d2);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_serialization_enter(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_serialization_leave (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_serialization_leave(domain->d2);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_serialization_leave(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_set_refinement (fclaw_domain_t * domain,
                             int smooth_refine, int smooth_level,
                             int coarsen_delay)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_set_refinement(domain->d2,smooth_refine,smooth_level,
                                      coarsen_delay);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_set_refinement(domain->d3,smooth_refine,smooth_level,
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
    if(domain->refine_dim == 2)
    {
        fclaw2d_patch_mark_refine(domain->d2,blockno,patchno);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_patch_mark_refine(domain->d3,blockno,patchno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_patch_mark_coarsen(fclaw_domain_t *domain, int blockno, int patchno)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_patch_mark_coarsen(domain->d2,blockno,patchno);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_patch_mark_coarsen(domain->d3,blockno,patchno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void
fclaw_domain_iterate_adapted(fclaw_domain_t *old_domain, fclaw_domain_t *new_domain, fclaw_match_callback_t mcb, void *user)
{
    FCLAW_ASSERT(old_domain->refine_dim == new_domain->refine_dim);

    fclaw_match_wrap_user_t mcb_wrap_user;
    mcb_wrap_user.mcb = mcb;
    mcb_wrap_user.user = user;

    if(old_domain->refine_dim == 2)
    {
        fclaw2d_domain_iterate_adapted(old_domain->d2,
                                       new_domain->d2,
                                       fclaw2d_match_callback_wrap,
                                       &mcb_wrap_user);
    }
    else if (old_domain->refine_dim == 3)
    {
        fclaw3d_domain_iterate_adapted(old_domain->d3,
                                       new_domain->d3,
                                       fclaw3d_match_callback_wrap,
                                       &mcb_wrap_user);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_allocate_before_partition(fclaw_domain_t *domain, size_t data_size, void ***patch_data)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_allocate_before_partition(domain->d2,data_size,patch_data);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_allocate_before_partition(domain->d3,data_size,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_retrieve_after_partition(fclaw_domain_t *domain, void ***patch_data)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_retrieve_after_partition(domain->d2,patch_data);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_retrieve_after_partition(domain->d3,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_iterate_partitioned(fclaw_domain_t *old_domain, fclaw_domain_t *new_domain, fclaw_transfer_callback_t tcb, void *user)
{
    FCLAW_ASSERT(old_domain->refine_dim == new_domain->refine_dim);

    fclaw_transfer_wrap_user_t wrap;
    wrap.tcb = tcb;
    wrap.user = user;

    if(old_domain->refine_dim == 2)
    {
        fclaw2d_domain_iterate_partitioned(old_domain->d2,
                                           new_domain->d2,
                                           fclaw2d_transfer_callback_wrap,
                                           &wrap);
    }
    else if(old_domain->refine_dim == 3)
    {
        fclaw3d_domain_iterate_partitioned(old_domain->d3,
                                           new_domain->d3,
                                           fclaw3d_transfer_callback_wrap,
                                           &wrap);
    }
    else 
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_free_after_partition(fclaw_domain_t *domain, void ***patch_data)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_free_after_partition(domain->d2,patch_data);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_free_after_partition(domain->d3,patch_data);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

fclaw_domain_exchange_t * 
fclaw_domain_allocate_before_exchange (fclaw_domain_t * domain,
                                       size_t data_size)
{
    FCLAW_ASSERT(domain->refine_dim == 2 || domain->refine_dim == 3);

    fclaw_domain_exchange_t * e = FCLAW_ALLOC(fclaw_domain_exchange_t, 1);
    e->refine_dim = domain->refine_dim;

    if(domain->refine_dim == 2)
    {
        e->d2 = fclaw2d_domain_allocate_before_exchange(domain->d2,data_size);
        e->d3 = NULL;
    }
    else
    {
        e->d2 = NULL;
        e->d3 = fclaw3d_domain_allocate_before_exchange(domain->d3,data_size);
    }

    return e;
}

void fclaw_domain_ghost_exchange (fclaw_domain_t * domain,
                                  fclaw_domain_exchange_t *e,
                                  int exchange_minlevel,
                                  int exchange_maxlevel)
{
    FCLAW_ASSERT(domain->refine_dim == 2 || domain->refine_dim == 3);
    FCLAW_ASSERT(e->refine_dim == domain->refine_dim);

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_ghost_exchange(domain->d2,
                                      e->d2,
                                      exchange_minlevel,
                                      exchange_maxlevel);
    }
    else
    {
        fclaw3d_domain_ghost_exchange(domain->d3,
                                      e->d3,
                                      exchange_minlevel,
                                      exchange_maxlevel);
    }
}

void fclaw_domain_ghost_exchange_begin (fclaw_domain_t * domain,
                                        fclaw_domain_exchange_t * e,
                                        int exchange_minlevel,
                                        int exchange_maxlevel)
{
    FCLAW_ASSERT(domain->refine_dim == 2 || domain->refine_dim == 3);
    FCLAW_ASSERT(e->refine_dim == domain->refine_dim);

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_ghost_exchange_begin(domain->d2,
                                            e->d2,
                                            exchange_minlevel,
                                            exchange_maxlevel);
    }
    else
    {
        fclaw3d_domain_ghost_exchange_begin(domain->d3,
                                            e->d3,
                                            exchange_minlevel,
                                            exchange_maxlevel);
    }
}

void fclaw_domain_ghost_exchange_end (fclaw_domain_t * domain,
                                      fclaw_domain_exchange_t * e)
{
    FCLAW_ASSERT(domain->refine_dim == 2 || domain->refine_dim == 3);
    FCLAW_ASSERT(e->refine_dim == domain->refine_dim);

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_ghost_exchange_end(domain->d2,
                                          e->d2);
    }
    else
    {
        fclaw3d_domain_ghost_exchange_end(domain->d3,
                                          e->d3);
    }
}

void fclaw_domain_free_after_exchange (fclaw_domain_t * domain,
                                       fclaw_domain_exchange_t * e)
{
    FCLAW_ASSERT(domain->refine_dim == 2 || domain->refine_dim == 3);
    FCLAW_ASSERT(e->refine_dim == domain->refine_dim);

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_free_after_exchange(domain->d2,
                                           e->d2);
    }
    else
    {
        fclaw3d_domain_free_after_exchange(domain->d3,
                                           e->d3);
    }

    FCLAW_FREE(e);
}

fclaw_domain_indirect_t
 * fclaw_domain_indirect_begin (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        return (fclaw_domain_indirect_t*) fclaw2d_domain_indirect_begin(domain->d2);
    }
    else if (domain->refine_dim == 3)
    {
        return (fclaw_domain_indirect_t*) fclaw3d_domain_indirect_begin(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_indirect_end (fclaw_domain_t * domain,
                                fclaw_domain_indirect_t* indirect)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_indirect_end(domain->d2,(fclaw2d_domain_indirect_t*)indirect);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_indirect_end(domain->d3,(fclaw3d_domain_indirect_t*)indirect);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

fclaw_patch_relation_t
fclaw_domain_indirect_neighbors (fclaw_domain_t * domain,
                                 fclaw_domain_indirect_t * indirect,
                                 int ghostno, int faceno, int rproc[],
                                 int *rblockno, int rpatchno[],
                                 int *rfaceno)
{
    if(domain->refine_dim == 2)
    {
        return (fclaw_patch_relation_t) fclaw2d_domain_indirect_face_neighbors(domain->d2,
                                                                               (fclaw2d_domain_indirect_t*)indirect,
                                                                               ghostno,faceno,rproc,rblockno,
                                                                               rpatchno,rfaceno);
    }
    else if (domain->refine_dim == 3)
    {
        return (fclaw_patch_relation_t) fclaw3d_domain_indirect_face_neighbors(domain->d3,
                                                                               (fclaw3d_domain_indirect_t*)indirect,
                                                                               ghostno,faceno,rproc,rblockno,
                                                                               rpatchno,rfaceno);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_domain_indirect_destroy (fclaw_domain_t * domain,
                                    fclaw_domain_indirect_t * indirect)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_indirect_destroy(domain->d2,(fclaw2d_domain_indirect_t*)indirect);
    }
    else if (domain->refine_dim == 3)
    {
        fclaw3d_domain_indirect_destroy(domain->d3,(fclaw3d_domain_indirect_t*)indirect);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}


int fclaw_domain_is_meta (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        return fclaw2d_domain_is_meta(domain->d2);
    }
    else if (domain->refine_dim == 3)
    {
        return fclaw3d_domain_is_meta(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}
