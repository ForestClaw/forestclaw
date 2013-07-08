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

#ifndef FORESTCLAW2D_H
#define FORESTCLAW2D_H

#include <sc.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct fclaw2d_domain fclaw2d_domain_t;
typedef struct fclaw2d_block fclaw2d_block_t;
typedef struct fclaw2d_patch fclaw2d_patch_t;

typedef enum
{
    FCLAW2D_PATCH_CHILDID = 0x7,
    FCLAW2D_PATCH_FIRST_SIBLING = 0x8,
    FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY = 0x10
}
fclaw2d_patch_flags_t;

/*
 * The domain structure gives a processor local view of the grid hierarchy.
 * Unless explicitly noted otherwise, all variables are processor local,
 * i.e., they are generally different on each processor.
 * Variables that are synchronized and shared between processors
 * are denoted *global*.
 */

struct fclaw2d_patch
{
    int level;                  /* 0 is root, increases if refined */
    int flags;                  /* flags that encode tree information */
    double xlower, xupper;
    double ylower, yupper;
    fclaw2d_patch_t *next;      /* next patch same level same block */
    void *user;
};

struct fclaw2d_block
{
    double xlower, xupper;
    double ylower, yupper;
    int is_boundary[4];         /* physical boundary flag */
    int num_patches;            /* local patches in this block */
    int num_patches_before;     /* in all previous blocks */
    int minlevel, maxlevel;     /* local over this block */
    fclaw2d_patch_t *patches;   /* allocated storage */
    fclaw2d_patch_t **patchbylevel;     /* array of pointers */
    void *user;
};

struct fclaw2d_domain
{
    MPI_Comm mpicomm;           /* MPI communicator */
    int mpisize, mpirank;       /* MPI variables */
    int num_patches_all;        /* sum over all blocks */
    int minlevel_all, maxlevel_all;     /* proc local */
    int global_minlevel, global_maxlevel;       /* global */
    int possible_maxlevel;      /* theoretical maximum */
    int num_ghost_patches;      /* off-proc patches relevant to this proc */
    int num_blocks;
    fclaw2d_block_t *blocks;    /* allocated storage */
    int *patch_to_block;        /* allocated storage */
    void *pp;                   /* opaque backend data */
    int pp_owned;               /* The pp member is owned by this domain */
    void *user;
};

/** Return the space dimension. */
int fclaw2d_domain_dimension (const fclaw2d_domain_t * domain);

/** Return the number of faces of a cube: 4 in 2D, 6 in 3D. */
int fclaw2d_domain_num_faces (const fclaw2d_domain_t * domain);

/** Return the number of corners of a cube: 4 in 2D, 8 in 3D.
 * This is the same as the number of siblings in a refined tree. */
int fclaw2d_domain_num_corners (const fclaw2d_domain_t * domain);

/** Return the number of corners of a cube face: 2 in 2D, 4 in 3D.
 * This is the same as the number of refined (smaller) face neighbors. */
int fclaw2d_domain_num_face_corners (const fclaw2d_domain_t * domain);

/** Return the number of possible orientations of a cube face.
 * This is mostly used for internal encodings.
 */
int fclaw2d_domain_num_orientations (const fclaw2d_domain_t * domain);

/** Find the numbers of faces adjacent to a cube corner: 2 in 2D, 3 in 3D. */
void fclaw2d_domain_corner_faces (const fclaw2d_domain_t * domain,
                                  int icorner, int faces[2]);

/** Return the dimension of a corner.
 * \param [in] patch    A patch with properly set member variables.
 * \param [in] cornerno A corner number in 0..3.
 * \return              0 if the corner is always at a fourfold intersection,
 *                      1 if the corner would end up in the middle of a face
 *                      when there is a coarser neighbor.
 */
int fclaw2d_patch_corner_dimension (const fclaw2d_patch_t * patch,
                                    int cornerno);

/** Return the number of a patch with respect to its parent in the tree.
 * \param [in] patch    A patch with properly set member variables.
 * \return              The child id is a number in 0..3.
 */
int fclaw2d_patch_childid (const fclaw2d_patch_t * patch);

/** Check if a patch is the first in a family of four siblings.
 * \param [in] patch    A patch with properly set member variables.
 * \return              True if patch is the first sibling.
 */
int fclaw2d_patch_is_first_sibling (const fclaw2d_patch_t * patch);

void *fclaw2d_alloc (size_t size);
void *fclaw2d_calloc (size_t nmemb, size_t size);
void *fclaw2d_realloc (void *ptr, size_t size);
void fclaw2d_free (void *ptr);
#define FCLAW2D_ALLOC(t,n)      (t *) fclaw2d_alloc ((n) * sizeof (t))
#define FCLAW2D_ALLOC_ZERO(t,n) (t *) fclaw2d_calloc ((n), sizeof (t))
#define FCLAW2D_REALLOC(p,t,n)  (t *) fclaw2d_realloc ((p), (n) * sizeof (t))
#define FCLAW2D_FREE(p)         fclaw2d_free (p)

/** Callback prototype for the patch iterators.
 * \param [in] domain	General domain structure.
 * \param [in] patch	The patch currently processed by the iterator.
 * \param [in] blockno  Block number of processed patch.
 * \param [in] patchno  Patch number within block of processed patch.
 * \param [in,out] user	Data that was passed into the iterator functions.
 */
typedef void (*fclaw2d_patch_callback_t)
    (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
     int blockno, int patchno, void *user);

/** Iterate over all patches on a given level.
 * \param [in] domain	General domain structure.
 * \param [in] level	Level to iterate.  Ignore patches of other levels.
 * \param [in] pcb	Function called for each patch of matching level.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw2d_domain_iterate_level (fclaw2d_domain_t * domain, int level,
                                   fclaw2d_patch_callback_t pcb, void *user);

/** Iterate over all patches of all levels.
 * \param [in] domain	General domain structure.
 * \param [in] pcb	Function called for each patch in the domain.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw2d_domain_iterate_patches (fclaw2d_domain_t * domain,
                                     fclaw2d_patch_callback_t pcb,
                                     void *user);

/** Iterate over all families of sibling patches.
 * \param [in] domain	General domain structure.
 * \param [in] pcb	Function called for each family in the domain.
 *                      Its patch argument points to an array of four
 *                      valid patches that constitute a family of siblings.
 *                      Their patchnos are consecutive, blockno is the same.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw2d_domain_iterate_families (fclaw2d_domain_t * domain,
                                      fclaw2d_patch_callback_t pcb,
                                      void *user);

/** Determine physical boundary status as 1, or 0 for neighbor patches.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in,out] boundaries	Domain boundary boolean flags.
 *			The order is left, right, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int fclaw2d_patch_boundary_type (fclaw2d_domain_t * domain,
                                 int blockno, int patchno, int boundaries[4]);

typedef enum fclaw2d_face_neighbor
{
    FCLAW2D_PATCH_BOUNDARY,
    FCLAW2D_PATCH_HALFSIZE,
    FCLAW2D_PATCH_SAMESIZE,
    FCLAW2D_PATCH_DOUBLESIZE
}
fclaw2d_patch_relation_t;

/** Determine neighbor patch(es) and orientation across a given face.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] faceno   Number of the patch face: left, right, bottom, top.
 * \param [out] rproc   Processor number of neighbor patches.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch numbers for up to 2 neighbors.
 *                              The patch number is relative to its block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghost_patches[.
 * \param [out] rfaceno Neighbor face number and orientation.
 * \return              The relative patch size of the face neighbor.
 */
fclaw2d_patch_relation_t fclaw2d_patch_face_neighbors (fclaw2d_domain_t *
                                                       domain, int blockno,
                                                       int patchno,
                                                       int faceno,
                                                       int rproc[2],
                                                       int *rblockno,
                                                       int rpatchno[2],
                                                       int *rfaceno);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  faceno      The number of the originating face.
 * \param [in]  rfaceno     Encoded as rfaceno = r * 4 + nf, where nf = 0..3 is
 *                          the neigbbor's connecting face number and r = 0..1
 *                          is the relative orientation to the neighbor's face.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face,
 *                          the first referring to the tangential and the second
 *                          to the normal.  A permutation of (0, 1).
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reversal flag for tangential axis (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 *              [1,4,7]     0 (unused for compatibility with 3D).
 */
void fclaw2d_patch_face_transformation (int faceno, int rfaceno,
                                        int ftransform[]);

/** Determine neighbor patch(es) and orientation across a given corner.
 * The current version only supports one neighbor, i.e. no true multi-block.
 * A query across a corner in the middle of a longer face returns the boundary.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] cornerno	Number of the patch corner: 0=bl, 1=br, 2=tl, 3=tr.
 * \param [out] rproc   Processor number of neighbor patch.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch number relative to the block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghosts_patches[.
 * \param [out] neighbor_size   The relative patch size of the neighbor.
 * \return			True if at least one corner neighbor exists.
 */
int fclaw2d_patch_corner_neighbors (fclaw2d_domain_t * domain,
                                    int blockno, int patchno, int cornerno,
                                    int *rproc, int *rblockno, int *rpatchno,
                                    fclaw2d_patch_relation_t * neighbor_size);

/** Mark a patch for refinement.
 * It is safe to call this function from an iterator callback except
 * fclaw2d_domain_iterate_adapted.
 */
void fclaw2d_patch_mark_refine (fclaw2d_domain_t * domain,
                                int blockno, int patchno);

/** Mark a patch for coarsening.
 * Coarsening will only happen if the patch family is not further refined
 * and all sibling patches are marked as well.
 * It is safe to call this function from an iterator callback except
 * fclaw2d_domain_iterate_adapted.
 */
void fclaw2d_patch_mark_coarsen (fclaw2d_domain_t * domain,
                                 int blockno, int patchno);

/** Callback prototype used in fclaw2d_domain_iterate_adapted.
 * The newsize value informs on refine/coarsen/noop status.
 * If refined (new patch is HALFSIZE), the old patch is old_patch[0] and the 
 * new patches are given by new_patch[0] through new_patch[3]. The new_patchno
 * numbers are consecutive as well.
 * If noop (new patch is SAMESIZE), only old_patch[0] and new_patch[0] matter.
 * If coarsened (new patch is DOUBLESIZE), situation is the reverse of refine.
 */
typedef void (*fclaw2d_match_callback_t) (fclaw2d_domain_t * old_domain,
                                          fclaw2d_patch_t * old_patch,
                                          fclaw2d_domain_t * new_domain,
                                          fclaw2d_patch_t * new_patch,
                                          fclaw2d_patch_relation_t newsize,
                                          int blockno,
                                          int old_patchno, int new_patchno,
                                          void *user);

/** Iterate over the previous and the adapted domain simultaneously.
 * \param [in,out] old_domain   Domain before adaptation.
 * \param [in,out] new_domain   Domain after adaptation.
 * \param [in] mcb              Callback 
 */
void fclaw2d_domain_iterate_adapted (fclaw2d_domain_t * old_domain,
                                     fclaw2d_domain_t * new_domain,
                                     fclaw2d_match_callback_t mcb,
                                     void *user);

/** Exchange data for parallel ghost neighbors.
 * This function gathers data from parallel neighbor (ghost) patches.
 * We assume that the data size for all patches is the same.
 * \param [in] domain           Used to access forest and ghost metadata.
 *                              #(local patches) is domain->num_patches_all.
 *                              #(off-proc patches) is domain->num_ghost_patches.
 * \param [in] data_size        Number of bytes per patch to transfer.
 * \param [in] patch_data       One pointer per processor-local patch as input.
 *                              It is only ever dereferenced for patches that
 *                              touch the parallel boundary from the inside, i.e.,
 *                              if (flags & FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY).
 *                              The pointers to the others' data can be NULL.
 * \param [in,out] ghost_data   One pointer per off-processor patch as output.
 *                              The memory must be pre-allocated.
 */
void fclaw2d_domain_parallel_exchange (fclaw2d_domain_t * domain,
                                       size_t data_size,
                                       void **patch_data, void **ghost_data);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
