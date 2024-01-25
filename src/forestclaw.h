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

#ifndef FORESTCLAW_H
#define FORESTCLAW_H

#include <fclaw_base.h>
#include <sc_keyvalue.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/** 
 * @file
 * Main dimension independent ForestClaw structures and routines
 */

/* ---------------------------------------------------------------------- */
///                      @name Data Types
/* ---------------------------------------------------------------------- */
///@{

typedef struct fclaw_domain_persist fclaw_domain_persist_t;
/** Typedef for fclaw_domain */
typedef struct fclaw_domain fclaw_domain_t;
/** Typedef for fclaw_block */
typedef struct fclaw_block fclaw_block_t;
/** Typedef for fclaw_patch */
typedef struct fclaw_patch fclaw_patch_t;

typedef struct fclaw_domain_exchange fclaw_domain_exchange_t;
/* Data structure for storing indirect parallel neighbor information */
typedef struct fclaw_domain_indirect fclaw_domain_indirect_t;

/* forward declare dimensioned patch types */
struct fclaw2d_patch;
struct fclaw3d_patch;

/** 
 * @brief The metadata structure for a forest leaf, which is a forestclaw patch.
 * The patch may be either a process-local patch or a ghost patch.
 */
struct fclaw_patch
{
    int refine_dim;                    /**< dimension */
    struct fclaw2d_patch* d2;          /**< 2D patch */
    struct fclaw3d_patch* d3;          /**< 3D patch */

    /** @{ @brief left/right coordinate */
    double xlower, xupper;
    /** @} */
    /** @{ @brief front/back coordinate */
    double ylower, yupper;
    /** @} */
    /** @{ @brief bottom/top coordinate. For 2D refinement, these are always set to 0,1 respectively. */
    double zlower, zupper;
    /** @} */
    int level;                  /**< 0 is root, increases if refined */

    void *user;                 /**< User Pointer */
};

/* forward declare dimensioned block types */
struct fclaw2d_block;
struct fclaw3d_block;

/**
 * @brief Data Structure for a block
 */
typedef struct fclaw_block
{
    int refine_dim;             /**< dimension */
    struct fclaw2d_block *d2;   /**< 2D block */
    struct fclaw3d_block *d3;   /**< 3D block */

    int num_patches;            /**< local patches in this block */
    int num_patches_before;     /**< in all previous blocks */
    double* vertices;     /**< for each block corner, the xyz coordinates
                                     of the connectivity structure */
    fclaw_patch_t *patches;           /**< The patches for this block */
    void *user;                       /**< User pointer */
} fclaw_block_t;

/* forward declare dimensioned domain types */
struct fclaw2d_domain;
struct fclaw3d_domain;

/**
 * @brief The domain structure is a collection of blocks
 * 
 * The domain structure gives a processor local view of the grid hierarchy.
 * Unless explicitly noted otherwise, all variables are processor local,
 * i.e., they are generally different on each processor.
 * Variables that are synchronized and shared between processors
 * are denoted *global*.
 */
struct fclaw_domain
{
    int refine_dim;                    /**< dimension */
    struct fclaw2d_domain *d2;         /**< 2D domain */
    struct fclaw3d_domain *d3;         /**< 3D domain */

    sc_MPI_Comm mpicomm;        /**< MPI communicator */
    int mpisize;                /**< MPI size */
    int mpirank;                /**< MPI rank */

    int local_num_patches;      /**< sum of patches over all blocks on this proc */
    int local_max_patches;      /**< maximum over all procs of local_num_patches */
    /** @{ */
    /** Local to proc.  If this proc doesn't
        store any patches at all, we set
        local_maxlevel < 0 <= local_minlevel. */
    int local_minlevel;
    int local_maxlevel;
    /** @} */
    int64_t global_num_patches; /**< sum of local_num_patches over all procs */
    int64_t global_num_patches_before;  /**< Number of patches on lower procs */
    int global_minlevel;       /**< global min level */
    int global_maxlevel;       /**< global max level */

    int num_blocks;             /**< Total number of blocks. */
    fclaw_block_t *blocks;    /**< allocated storage */
    int num_exchange_patches;   /**< number my patches relevant to other procs.
                                   Identified by this expression to be true:
                                   (patch->flags &
                                   FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY) */
    fclaw_patch_t **exchange_patches; /**< explicitly store exchange patches */
    int num_ghost_patches;      /**< number of off-proc patches relevant to this proc */
    fclaw_patch_t *ghost_patches;     /**< array of off-proc patches */

    /* formerly in domain_data_t */

    int count_set_patch;
    int count_delete_patch;

    fclaw_domain_exchange_t* exchange;
    fclaw_domain_indirect_t* indirect;

    void *user; /**< user data pointer */
};

///@}
/* ---------------------------------------------------------------------- */
///                  @name Topological Properties
/* ---------------------------------------------------------------------- */
///@{

/** Return the space dimension. */
int fclaw_domain_dimension (const fclaw_domain_t * domain);

/** Return the number of siblings a node has in a domain.
 *  4 in 2d, 8 in 3d. */
int fclaw_domain_num_siblings (const fclaw_domain_t * domain);

/** Return the number of faces of a cube: 4 in 2D, 6 in 3D. */
int fclaw_domain_num_faces (const fclaw_domain_t * domain);

/** Return the number of edges of a cube: 0 in 2D, 12 in 3D. */
int fclaw_domain_num_edges (const fclaw_domain_t * domain);

/** Return the number of corners of a cube: 4 in 2D, 8 in 3D.
 * This is the same as the number of siblings in a refined tree. */
int fclaw_domain_num_corners (const fclaw_domain_t * domain);

/** Return the number of corners of a cube face: 2 in 2D, 4 in 3D.
 * This is the same as the number of refined (smaller) face neighbors. */
int fclaw_domain_num_face_corners (const fclaw_domain_t * domain);

/** Return the number of possible orientations of a cube face.
 * This is mostly used for internal encodings.
 */
int fclaw_domain_num_orientations (const fclaw_domain_t * domain);

/** Find the numbers of faces adjacent ot a cube edge: 0 in 2D, 3 in 3D */
void fclaw_domain_edge_faces(const fclaw_domain_t *domain, int iedge, int faces[2]);
/** Find the numbers of faces adjacent to a cube corner: 2 in 2D, 3 in 3D. */
void fclaw_domain_corner_faces (const fclaw_domain_t * domain,
                                int icorner, int faces[3]);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Patch Functions
/* ---------------------------------------------------------------------- */
///@{

/** Return the dimension of a corner.
 * This function is LEGAL to call for both local and ghost patches.
 * \param [in] patch    A patch with properly set member variables.
 * \param [in] cornerno A corner number in 0..3 for 2D or 0..7 for 3D.
 * \return              0 if the corner is always at a fourfold intersection,
 *                      1 if the corner would end up in the middle of a face
 *                      when there is a coarser neighbor.
 */
int fclaw_patch_corner_dimension (const fclaw_patch_t * patch,
                                    int cornerno);

/** Return the number of a patch with respect to its parent in the tree.
 * This function is LEGAL to call for both local and ghost patches.
 * \param [in] patch    A patch with properly set member variables.
 * \return              The child id is a number in 0..3 for 2D or 0..7 for 3D.
 */
int fclaw_patch_childid (const fclaw_patch_t * patch);

/** Check if a patch is the first in a family of four siblings.
 * For ghost patches, we always return false.
 * \param [in] patch    A patch with properly set member variables.
 * \return              True if patch is the first sibling.
 */
int fclaw_patch_is_first_sibling (const fclaw_patch_t * patch);

/** Check whether a patch is a parallel ghost patch.
 * \param [in] patch    A patch with properly set member variables.
 * \return              True if patch is off-processor, false if local.
 */
int fclaw_patch_is_ghost (const fclaw_patch_t * patch);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Patch Iterators
/* ---------------------------------------------------------------------- */
///@{

/** Callback prototype for the patch iterators.
 * We iterate over local patches only.
 * \param [in] domain	General domain structure.
 * \param [in] patch	The local patch currently processed by the iterator.
 * \param [in] blockno  Block number of processed patch.
 * \param [in] patchno  Patch number within block of processed patch.
 * \param [in,out] user	Data that was passed into the iterator functions.
 */
typedef void (*fclaw_patch_callback_t)
    (fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user);

/** Iterate over all local patches on a given level.
 * \param [in] domain	General domain structure.
 * \param [in] level	Level to iterate.  Ignore patches of other levels.
 * \param [in] pcb	Function called for each patch of matching level.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                                 fclaw_patch_callback_t pcb, void *user);

/** Iterate over all local patches of all levels.
 * \param [in] domain	General domain structure.
 * \param [in] pcb	Function called for each patch in the domain.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw_domain_iterate_patches (fclaw_domain_t * domain,
                                   fclaw_patch_callback_t pcb,
                                   void *user);

/** Iterate over all families of local sibling patches.
 * \param [in] domain	General domain structure.
 * \param [in] pcb	Function called for each family in the domain.
 *                      Its patch argument points to an array of four
 *                      valid patches that constitute a family of siblings.
 *                      Their patchnos are consecutive, blockno is the same.
 * \param [in,out] user	Data is passed to the pcb callback.
 */
void fclaw_domain_iterate_families (fclaw_domain_t * domain,
                                    fclaw_patch_callback_t pcb,
                                    void *user);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Patch Neighbors
/* ---------------------------------------------------------------------- */
///@{

/** Determine physical boundary status as 1, or 0 for neighbor patches.
 * This must ONLY be called for local patches.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in,out] boundaries	Domain boundary boolean flags 
            (length 4 for 2D, 6 for 3D).
 *			The order is left, right, front, back, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int fclaw_patch_boundary_type (fclaw_domain_t * domain,
                               int blockno, int patchno, int boundaries[]);

/** Determine whether the normal to a face neighbor align.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in] faceno   Number of the face of the patch.
 * \return		True if normals match, false for mismatch.
 */
int fclaw_patch_normal_match (fclaw_domain_t * domain,
                              int blockno, int patchno, int faceno);

/**
 * @brief The type of face neighbor
 */
typedef enum fclaw_face_neighbor
{
    /** Physical boundary */
    FCLAW_PATCH_BOUNDARY,
    /** Half-size (finer) neighbor */
    FCLAW_PATCH_HALFSIZE,
    /** Same-size neighbor */
    FCLAW_PATCH_SAMESIZE,
    /** Double-size (coarser) neighbor */
    FCLAW_PATCH_DOUBLESIZE
}
fclaw_patch_relation_t;

/** Determine neighbor patch(es) and orientation across a given face.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] faceno   Number of the patch face: left, right, bottom, top.
 * \param [out] rproc   Processor number of neighbor patches (length 2 in 2D, 4 in 3D).  
 *                      Exception:
 *                      If the neighbor is a bigger patch, rproc[1] contains
 *                      the number of the small patch as one of the smaller faces.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch numbers (length 2 in 2D, 4 in 3D)
 *                              The patch number is relative to its block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghost_patches[.
 * \param [out] rfaceno Neighbor face number and orientation.
 * \return              The relative patch size of the face neighbor.
 */
fclaw_patch_relation_t fclaw_patch_face_neighbors (fclaw_domain_t *
                                                   domain, int blockno,
                                                   int patchno,
                                                   int faceno,
                                                   int rproc[],
                                                   int *rblockno,
                                                   int rpatchno[],
                                                   int *rfaceno);

/** Change perspective across a face neighbor situation.
 * \param [in]     dim      Dimension [2, 3].
 * \param [in,out] faceno   On input, valid face number for a patch.
 *                          On output, valid face number seen from
 *                          faceno's neighbor patch.
 * \param [in,out] rfaceno  On input, encoded neighbor face number as returned
 *                          by fclaw_patch_face_neighbors.
 *                          On output, encoded neighbor face number seen from
 *                          faceno's neighbor patch.
 */
void fclaw_patch_face_swap (int dim, int *faceno, int *rfaceno);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  dim         Dimension [2, 3].
 * \param [in]  faceno      The number of the originating face.
 * \param [in]  rfaceno     In 2D: 
 *                          Encoded as rfaceno = r * 4 + nf, where nf = 0..3 is
 *                          the neigbbor's connecting face number and r = 0..1
 *                          is the relative orientation to the neighbor's face.
 *
 *                          In 3D:
 *                          Encoded as rfaceno = r * 6 + nf, where nf = 0..3 is
 *                          the neigbbor's connecting face number and r = 0..1
 *                          is the relative orientation to the neighbor's face.
 * \param [out] ftransform  This array holds 9 integers.
 *                          In 2D:
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
 *                          [8] & 4: Both patches are in the same block,
 *                                   the \a ftransform contents are ignored.
 *              [1,4,7]     0 (unused for compatibility with 3D).
 *
 *                          In 3D:
 *              [0..2]      The coordinate axis sequence of the origin face,
 *                          the first two referring to the tangentials and the
 *                          third to the normal.  A permutation of (0, 1, 2).
 *              [3..5]      The coordinate axis sequence of the target face.
 *              [6..8]      Edge reversal flags for tangential axes (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 *                          [8] & 4: Both patches are in the same block,
 *                                   the \a ftransform contents are ignored.
 */
void fclaw_patch_face_transformation (int dim, int faceno, int rfaceno,
                                      int ftransform[]);

/** Modify an existing face transformation depending on intra-block usage.
 * This function can be called any number of times on the same transform array.
 * \param [in] dim              Dimension [2, 3].
 * \param [in,out] ftransform   Array of values necessarily created by \ref
 *                              fclaw_patch_face_transformation.
 * \param [in] sameblock        Transformation supposed to work in same block?
 */
void fclaw_patch_face_transformation_block (int dim, int ftransform[],
                                            int sameblock);

/** Fill an array with the axis combination of a face neighbor transformation
 * that operates on two patches in the same block (the trivial case).
 * Use when there is no prior call to \ref fclaw_patch_face_transformation.
 * Don't postprocess the result any further -- it's only useful intra-block.
 * \param [in] dim              Dimension [2, 3].
 * \param [out] ftransform      Gets initialized to a same-block transform.
 */
void fclaw_patch_face_transformation_intra (int dim, int ftransform[]);

/** Return whether a face transformation is valid.
 * \param [in] dim              Dimension [2, 3].
 * \param [in] ftransform       Array of values as created by \ref
 *                              fclaw_patch_face_transformation,
 *                              possibly modified by \ref
 *                              fclaw_patch_face_transformation_block.
 * \return                      True if valid, false if not.
 */
int fclaw_patch_face_transformation_valid (int dim, const int ftransform[]);

/** Transform a 2D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the two patches are of the SAME size.
 * If the neighbor patch is in the same block we must set (ftransform[8] & 4).
 * Else we have an input patch in one block and on output patch across a face.
 * It is LEGAL to call this function for both local and ghost patches.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] ftransform   It must have room for NINE (9) integers and be
 *                          computed by \a fclaw_patch_face_transformation.
 *                          If \a ipatch and \a opatch are in the same block,
 *                          we require \a ftransform[8] |= 4.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 *                          The number of cells must match according to the face
 *                          transformation.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        Integer coordinate along x-axis in \a based .. \a mx.
 * \param [in,out] j        Integer coordinate along y-axis in \a based .. \a my.
 */
void fclaw_patch_2d_transform_face (fclaw_patch_t * ipatch,
                                    fclaw_patch_t * opatch,
                                    const int ftransform[],
                                    int mx, int my, int based, int *i, int *j);

/** Transform a 2D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size).
 * If the neighbor patch is in the same block we must set (ftransform[8] & 4).
 * Else we have an input patch in one block and on output patch across a face.
 * It is LEGAL to call this function for both local and ghost patches.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] ftransform   It must have room for NINE (9) integers and be
 *                          computed by \a fclaw_patch_face_transformation.
 *                          If \a ipatch and \a opatch are in the same block,
 *                          we require \a ftransform[8] |= 4.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 *                          The number of cells must match according to the face
 *                          transformation.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        FOUR (4) integer coordinates along x-axis in
 *                          \a based .. \a mx.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] j        FOUR (4) integer coordinates along y-axis in
 *                          \a based .. \a my.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 */
void fclaw_patch_2d_transform_face2 (fclaw_patch_t * ipatch,
                                     fclaw_patch_t * opatch,
                                     const int ftransform[],
                                     int mx, int my, int based, int i[],
                                     int j[]);

/** Transform a 3D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the two patches are of the SAME size.
 * If the neighbor patch is in the same block we must set (ftransform[8] & 4).
 * Else we have an input patch in one block and on output patch across a face.
 * It is LEGAL to call this function for both local and ghost patches.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] ftransform   It must have room for NINE (9) integers and be
 *                          computed by \a fclaw3d_patch_face_transformation.
 *                          If \a ipatch and \a opatch are in the same block,
 *                          we require \a ftransform[8] |= 4.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] mz           Number of cells along z direction of patch.
 *                          The number of cells must match according to the face
 *                          transformation.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        Integer coordinate along x-axis in \a based .. \a mx.
 * \param [in,out] j        Integer coordinate along y-axis in \a based .. \a my.
 * \param [in,out] k        Integer coordinate along z-axis in \a based .. \a mz.
 */
void fclaw_patch_3d_transform_face (fclaw_patch_t * ipatch,
                                    fclaw_patch_t * opatch,
                                    const int ftransform[],
                                    int mx, int my, int mz, int based,
                                    int *i, int *j, int *k);

/** Transform a 3D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size).
 * If the neighbor patch is in the same block we must set (ftransform[8] & 4).
 * Else we have an input patch in one block and on output patch across a face.
 * It is LEGAL to call this function for both local and ghost patches.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] ftransform   It must have room for NINE (9) integers and be
 *                          computed by \a fclaw3d_patch_face_transformation.
 *                          If \a ipatch and \a opatch are in the same block,
 *                          we require \a ftransform[8] |= 4.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] mz           Number of cells along z direction of patch.
 *                          The number of cells must match according to the face
 *                          transformation.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        EIGHT (8) integer coordinates along x-axis in
 *                          \a based .. \a mx.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] j        EIGHT (8) integer coordinates along y-axis in
 *                          \a based .. \a my.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] k        EIGHT (8) integer coordinates along z-axis in
 *                          \a based .. \a mz.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 */
void fclaw_patch_3d_transform_face2 (fclaw_patch_t * ipatch,
                                     fclaw_patch_t * opatch,
                                     const int ftransform[],
                                     int mx, int my, int mz, int based,
                                     int i[], int j[], int k[]);

/** Determine neighbor patch(es) and orientation across a given edge.
 * The current version only supports one neighbor, i.e., no true multi-block.
 * A query across an edge in the middle of a longer face returns the boundary.
 * We only return edge neighbors that are not already face neighbors.
 * Inter-tree edges are only returned if the number of meeting edges is
 * exactly four.  Five or more are currently not supported.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] edgeno	Number of the patch edge: 4 parallel to x axis,
                        then 4 parallel to y axis, then 4 parallel to z.
 * \param [out] rproc   Processor number of neighbor patch.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch number relative to the block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghosts_patches[.
 * \param [out] redge           Number of the edge from the other neighbor.
 * \param [out] neighbor_size   The relative patch size of the neighbor.
 * \return                      True if at least one edge neighbor exists
 *                              that is not already a face neighbor.
 */
int fclaw_patch_edge_neighbors (fclaw_domain_t * domain,
                                int blockno, int patchno, int edgeno,
                                int *rproc, int *rblockno, int *rpatchno,
                                int *redge,
                                fclaw_patch_relation_t * neighbor_size);

/** Change perspective across an edge neighbor situation.
 * \param [in,out] edgeno       On input, valid edge number for a patch.
 *                              On output, edge number seen from
 *                              the edge neighbor patch.
 * \param [in,out] redgeno      On input, valid edge number as returned
 *                              by fclaw_patch_edge_neighbors.
 *                              On output, edge number seen from
 *                              the edge neighbor patch.
 */
void fclaw_patch_edge_swap (int *edgeno, int *redgeno);

/** Transform a patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the two patches are of the SAME size and that the
 * patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * It is ILLEGAL to call this function for patches from face-neighboring blocks.
 * Use \ref fclaw_patch_3d_transform_face for such patches instead.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] iedge        Edge number of this patch to transform across.
 *                          This function assumes oedge == iedge ^ 3, so
 *                          oedge is the edge opposite of iedge.
 * \param [in] is_block_boundary      Set to true for a block edge.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] mz           Number of cells along z direction of patch.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        Integer coordinate along x-axis in \a based .. \a mx.
 * \param [in,out] j        Integer coordinate along y-axis in \a based .. \a my.
 * \param [in,out] k        Integer coordinate along z-axis in \a based .. \a mz.
 */
void fclaw_patch_3d_transform_edge (fclaw_patch_t * ipatch,
                                    fclaw_patch_t * opatch,
                                    int iedge, int is_block_boundary,
                                    int mx, int my, int mz,
                                    int based, int *i, int *j, int *k);

/** Transform a patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size) and that
 * the patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * It is ILLEGAL to call this function for patches from face-neighboring blocks.
 * Use \ref fclaw_patch_3d_transform_face for such patches instead.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] iedge        Edge number of this patch to transform across.
 *                          This function assumes oedge == iedge ^ 3, so
 *                          oedge is the edge opposite of iedge.
 * \param [in] is_block_boundary      Set to true for a block edge.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] mz           Number of cells along z direction of patch.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        EIGHT (8) integer coordinates along x-axis in
 *                          \a based .. \a mx.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] j        EIGHT (8) integer coordinates along y-axis in
 *                          \a based .. \a my.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] k        EIGHT (8) integer coordinates along y-axis in
 *                          \a based .. \a mz.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 */
void fclaw_patch_3d_transform_edge2 (fclaw_patch_t * ipatch,
                                     fclaw_patch_t * opatch,
                                     int iedge, int is_block_boundary,
                                     int mx, int my, int mz, int based,
                                     int i[], int j[], int k[]);

/** Determine neighbor patch(es) and orientation across a given corner.
 * The current version only supports one neighbor, i.e., no true multi-block.
 * A query across a corner in the middle of a longer face returns the boundary.
 * We only return corner neighbors that are not already face neighbors.
 * Inter-tree corners are only returned if the number of meeting corners is
 * exactly four.  Five or more are currently not supported.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] cornerno	Number of the patch corner:
 *                          In 2D: 0=bl, 1=br, 2=tl, 3=tr.
 *                          In 3D: 0=bfl, ..., 7=tbr.
 * \param [out] rproc   Processor number of neighbor patch.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch number relative to the block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghosts_patches[.
 * \param [out] rcorner         Number of the corner from the other neigbor.
 * \param [out] neighbor_size   The relative patch size of the neighbor.
 * \return                      True if at least one corner neighbor exists
 *                              that is not already a face neighbor.
 */
int fclaw_patch_corner_neighbors (fclaw_domain_t * domain,
                                  int blockno, int patchno, int cornerno,
                                  int *rproc, int *rblockno, int *rpatchno,
                                  int *rcorner,
                                  fclaw_patch_relation_t * neighbor_size);

/** Change perspective across a corner neighbor situation.
 * \param [in]     dim          Dimension [2, 3].
 * \param [in,out] cornerno     On input, valid corner number for a patch.
 *                              On output, corner number seen from
 *                              the corner neighbor patch.
 * \param [in,out] rcornerno    On input, valid corner number as returned
 *                              by fclaw_patch_face_neighbors.
 *                              On output, corner number seen from
 *                              the corner neighbor patch.
 */
void fclaw_patch_corner_swap (int dim, int *cornerno, int *rcornerno);

/** Transform a 2D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the two patches are of the SAME size and that the
 * patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * It is ILLEGAL to call this function for patches from face-neighboring blocks.
 * Use \ref fclaw2d_patch_transform_face for such patches instead.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] icorner      Corner number of this patch to transform across.
 *                          This function assumes ocorner == icorner ^ 3, so
 *                          ocorner is the opposite corner of icorner.
 * \param [in] is_block_boundary      Set to true for a block corner.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        Integer coordinate along x-axis in \a based .. \a mx.
 * \param [in,out] j        Integer coordinate along y-axis in \a based .. \a my.
 */
void fclaw_patch_2d_transform_corner (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my,
                                      int based, int *i, int *j);

/** Transform a 2D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size) and that
 * the patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * It is ILLEGAL to call this function for patches from face-neighboring blocks.
 * Use \ref fclaw2d_patch_transform_face2 for such patches instead.
 * \param [in] ipatch       The patch that the input coordinates are relative to.
 * \param [in] opatch       The patch that the output coordinates are relative to.
 * \param [in] icorner      Corner number of this patch to transform across.
 *                          This function assumes ocorner == icorner ^ 3, so
 *                          ocorner is the opposite corner of icorner.
 * \param [in] is_block_boundary      Set to true for a block corner.
 * \param [in] mx           Number of cells along x direction of patch.
 * \param [in] my           Number of cells along y direction of patch.
 * \param [in] based        Indices are 0-based for corners and 1-based for cells.
 * \param [in,out] i        FOUR (4) integer coordinates along x-axis in
 *                          \a based .. \a mx.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] j        FOUR (4) integer coordinates along y-axis in
 *                          \a based .. \a my.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 */
void fclaw_patch_2d_transform_corner2 (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int based,
                                       int i[], int j[]);

/** Transform a 3D patch coordinate into a neighbor patch's coordinate system.
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
                                      int based, int *i, int *j, int *k);

/** Transform a 3D patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size) and that
 * the patches lie in coordinate systems with the same orientation.
 * It is LEGAL to call this function for both local and ghost patches.
 * It is ILLEGAL to call this function for patches from face- or
 * edge-neighboring blocks. Use \ref fclaw3d_patch_transform_face2 or
 * \ref fclaw_patch_3d_transform_edge2 for such patches instead.
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
 * \param [in,out] i        EIGHT (8) integer coordinates along x-axis in
 *                          \a based .. \a mx.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] j        EIGHT (8) integer coordinates along y-axis in
 *                          \a based .. \a my.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 * \param [in,out] k        EIGHT (8) integer coordinates along y-axis in
 *                          \a based .. \a mz.  On input, only the first is used.
 *                          On output, they are relative to the fine patch and
 *                          stored in order of the children of the coarse patch.
 */
void fclaw_patch_3d_transform_corner2 (fclaw_patch_t * ipatch,
                                       fclaw_patch_t * opatch,
                                       int icorner, int is_block_boundary,
                                       int mx, int my, int mz, int based,
                                       int i[], int j[], int k[]);
///@}
/* ---------------------------------------------------------------------- */
///                         @name Adaptivity
/* ---------------------------------------------------------------------- */
///@{

/** Set parameters of refinement strategy in a domain.
 * This function only needs to be called once, and only for the first domain
 * created in the program.  The values of the parameters are automatically
 * transferred on adaptation and partitioning.
 * \param [in,out] domain       This domain's refinement strategy is set.
 * \param [in] smooth_refine    Activate or deactivete refinement smoothing.
 *                              A newly created domain has this set to false.
 * \param [in] smooth_level     If \b smooth_refine is true, denotes the
 *                              lowest level that activates the smoothing.
 *                              Use zero for smoothing across all levels.
 * \param [in] coarsen_delay    Non-negative number to set the delay for
 *                              coarsening after a patch has been last refined.
 *                              This number is a global threshold that is compared
 *                              against each patch's individual counter.
 */
void fclaw_domain_set_refinement (fclaw_domain_t * domain,
                                  int smooth_refine, int smooth_level,
                                  int coarsen_delay);

/** Mark a patch for refinement.
 * This must ONLY be called for local patches.
 * It is safe to call this function from an iterator callback except
 * \ref fclaw_domain_iterate_adapted.
 */
void fclaw_patch_mark_refine (fclaw_domain_t * domain,
                              int blockno, int patchno);

/** Mark a patch for coarsening.
 * This must ONLY be called for local patches.
 * Coarsening will only happen if the patch family is not further refined
 * and all sibling patches are marked as well.
 * It is safe to call this function from an iterator callback except
 * \ref fclaw_domain_iterate_adapted.
 */
void fclaw_patch_mark_coarsen (fclaw_domain_t * domain,
                               int blockno, int patchno);

/** Callback prototype used in fclaw_domain_iterate_adapted.
 * The newsize value informs on refine/coarsen/noop status.
 * If refined (new patch is HALFSIZE), the old patch is old_patch[0] and the
 * new patches are given by new_patch[0] through new_patch[3]. The new_patchno
 * numbers are consecutive as well.
 * If noop (new patch is SAMESIZE), only old_patch[0] and new_patch[0] matter.
 * If coarsened (new patch is DOUBLESIZE), situation is the reverse of refine.
 * We iterate over local patches only.
 */
typedef void (*fclaw_match_callback_t) (fclaw_domain_t * old_domain,
                                        fclaw_patch_t * old_patch,
                                        fclaw_domain_t * new_domain,
                                        fclaw_patch_t * new_patch,
                                        fclaw_patch_relation_t newsize,
                                        int blockno,
                                        int old_patchno, int new_patchno,
                                        void *user);

/** Iterate over the previous and the adapted domain simultaneously.
 * We iterate over local patches only.
 * \param [in,out] old_domain   Domain before adaptation.
 * \param [in,out] new_domain   Domain after adaptation.
 * \param [in] mcb              Callback.
 * \param [in,out] user         This pointer is passed to the callback.
 */
void fclaw_domain_iterate_adapted (fclaw_domain_t * old_domain,
                                   fclaw_domain_t * new_domain,
                                   fclaw_match_callback_t mcb,
                                   void *user);

///@}
/* ---------------------------------------------------------------------- */
///                         @name Parititon
/* ---------------------------------------------------------------------- */
///@{

/** Allocate data buffer for parallel transfer of all patches.
 * \param [in,out] domain       The memory lives inside this domain.
 * \param [in] data_size        Number of bytes per patch to transfer.
 * \param [in,out] patch_data   Address of an array of void pointers.
 *                              Data is allocated by this function.  After the
 *                              call, *patch_data holds one pointer per patch
 *                              that points to exactly data_size bytes of
 *                              memory that can be written to by forestclaw.
 *                              *patch_data must be NULL before the call.
 */
void fclaw_domain_allocate_before_partition (fclaw_domain_t * domain,
                                             size_t data_size,
                                             void ***patch_data);

/** Reallocate data buffer to reflect patch data after partition.
 * \param [in,out] domain       The memory lives inside this domain.
 * \param [in,out] patch_data   Address of an array of void pointers.
 *                              Data is reallocated by this function.  After the
 *                              call, *patch_data holds one pointer per patch
 *                              that points to exactly data_size bytes of
 *                              memory that can be read from by forestclaw.
 */
void fclaw_domain_retrieve_after_partition (fclaw_domain_t * domain,
                                            void ***patch_data);

/** Callback to iterate through the partitions.
 * We traverse every patch in the new partition.  If that patch was already
 * on the local processor before the partition, we identify its memory.
 * We iterate over local patches only.
 * \param [in,out] old_domain   Domain before partition.
 * \param [in,out] old_patch    If the patch stayed local, this is the pointer
 *                              in reference to the old domain and partition.
 *                              Otherwise, this patch pointer is set to NULL.
 * \param [in,out] new_domain   Domain after partition.
 * \param [in,out] new_patch    Patch in the new domain and partition.
 * \param [in] blockno          Number of the current block.
 * \param [in] old_patchno      Number of the patch that stayed local wrt. the
 *                              old domain and partition.  Minus one otherwise.
 * \param [in] new_patchno      Number of the iterated patch wrt. the new
 *                              domain and partition.
 * \param [in,out] user         Pointer passed to \ref
 *                              fclaw_domain_iterate_partitioned.
 */
typedef void (*fclaw_transfer_callback_t) (fclaw_domain_t * old_domain,
                                           fclaw_patch_t * old_patch,
                                           fclaw_domain_t * new_domain,
                                           fclaw_patch_t * new_patch,
                                           int blockno,
                                           int old_patchno, int new_patchno,
                                           void *user);

/** Iterate over the previous and partitioned domain simultaneously.
 * We iterate over local patches only.
 * \param [in,out] old_domain   Domain before partition.
 * \param [in,out] new_domain   Domain after partition.
 * \param [in] tcb              Callback.
 * \param [in,out] user         This pointer is passed to the callback.
 */
void fclaw_domain_iterate_partitioned (fclaw_domain_t * old_domain,
                                       fclaw_domain_t * new_domain,
                                       fclaw_transfer_callback_t tcb,
                                       void *user);

/** Free buffers that were used in transfering data during partition.
 * \param [in,out] domain       The memory lives inside this domain.
 * \param [in,out] patch_data   Address of an array of void pointers to free.
 *                              *patch_data will be NULL after the call.
 */
void fclaw_domain_free_after_partition (fclaw_domain_t * domain,
                                        void ***patch_data);

///@}
/* ---------------------------------------------------------------------- */
///                     @name Communication
/* ---------------------------------------------------------------------- */
///@{

/** Compute and return the maximum over all processors of a double value.
 * The minimum can be computed by using this function on the negative value.
 */
double fclaw_domain_global_maximum (fclaw_domain_t * domain, double d);

/** Compute and return the sum over all processors of a double value.
 */
double fclaw_domain_global_sum (fclaw_domain_t * domain, double d);

/** Synchronize all processes.  Avoid using if at all possible.
 */
void fclaw_domain_barrier (fclaw_domain_t * domain);

/** Serialize a section of code.
 * THIS IS NOT SCALABLE.
 * WILL BE HORRIBLY SLOW FOR LARGE NUMBERS OF PROCESSORS.
 * A processor returns from this function only after all lower-numbered
 * processors have called fclaw_domain_serialization_leave.
 * No collective communication routines must be called between the calls
 * to this function and fclaw_domain_serialization_leave.
 * \param [in] domain           The domain is not modified.
 */
void fclaw_domain_serialization_enter (fclaw_domain_t * domain);

/** Serialize a section of code.
 * THIS IS NOT SCALABLE.
 * WILL BE HORRIBLY SLOW FOR LARGE NUMBERS OF PROCESSORS.
 * A processor must call this function to allow all higher-numbered
 * processors to return from fclaw_domain_serialization_enter.
 * \param [in] domain           The domain is not modified.
 */
void fclaw_domain_serialization_leave (fclaw_domain_t * domain);

///@}
/* ---------------------------------------------------------------------- */
///                         @name Exchange
/* ---------------------------------------------------------------------- */
///@{

/**
 * @brief Simple wrapper struct for fclaw2d_exchange_info_t 
 *        and fclaw3d_exchange_info_t
 */
struct fclaw_domain_exchange
{
    int refine_dim; /**< refinement dimension */
    struct fclaw2d_domain_exchange *d2; /**< 2d exchange data */
    struct fclaw3d_domain_exchange *d3; /**< 3d exchange data */
};

/** Allocate buffer to hold the data from off-processor patches.
 * Free this by fclaw_domain_free_after_exchange before regridding.
 * \param [in] domain           The domain is not modified.
 * \param [in] data_size        Number of bytes per patch to exchange.
 * \return                      Allocated data structure.
 *                              The pointers in patch_data[i] need to be set
 *                              after this call by forestclaw.
 */
fclaw_domain_exchange_t *
fclaw_domain_allocate_before_exchange (fclaw_domain_t * domain,
                                       size_t data_size);

/** Exchange data for parallel ghost neighbors.
 * This function receives data from parallel neighbor (ghost) patches.
 * It can be called multiple times on the same allocated buffers.
 * We assume that the data size for all patches is the same.
 * \param [in] domain           Used to access forest and ghost metadata.
 *                              #(sent patches) is domain->num_exchange_patches.
 *                              #(received patches) is domain->num_ghost_patches.
 * \param [in] e                Allocated buffers whose e->patch_data[i] pointers
 *                              must have been set properly by forestclaw.
 * \param [in] exchange_minlevel The minimum quadrant level that is exchanged.
 * \param [in] exchange_maxlevel The maximum quadrant level that is exchanged.
 */
void fclaw_domain_ghost_exchange (fclaw_domain_t * domain,
                                  fclaw_domain_exchange_t * e,
                                  int exchange_minlevel,
                                  int exchange_maxlevel);

/** Start asynchronous exchange of parallel ghost neighbors.
 * The arguments are the same as for fclaw_domain_ghost_exchange.
 * It must be followed by a call to fclaw_domain_ghost_exchange_end.
 * Between begin and end, neither of \ref fclaw_domain_indirect_begin
 * and _end must be called.
 * \param [in,out] e            Its ghost_data member must survive and not
 *                              be written to until the completion of
 *                              fclaw_domain_ghost_exchange_end.
 *                              Its patch_data member may already be
 *                              overwritten after this function returns.
 */
void fclaw_domain_ghost_exchange_begin (fclaw_domain_t * domain,
                                        fclaw_domain_exchange_t * e,
                                        int exchange_minlevel,
                                        int exchange_maxlevel);

/** Complete asynchronous exchange of parallel ghost neighbors.
 * Must be called at some point after fclaw_domain_ghost_exchange_begin.
 * Between begin and end, neither of \ref fclaw_domain_indirect_begin
 * and _end must be called.
 * \param [in,out] e            Its ghost_data member must have survived.
 */
void fclaw_domain_ghost_exchange_end (fclaw_domain_t * domain,
                                      fclaw_domain_exchange_t * e);

/** Free buffers used in exchanging off-processor data during time stepping.
 * This should be done just before regridding.
 * \param [in] domain           The domain is not modified.
 * \param [in] e                Allocated buffers.
 */
void fclaw_domain_free_after_exchange (fclaw_domain_t * domain,
                                       fclaw_domain_exchange_t * e);

///@}
/* ---------------------------------------------------------------------- */
///                 @name Indirect Parallel Neighbors
/* ---------------------------------------------------------------------- */
///@{


/** Begin sending messages to determine neighbors of ghost patches.
 * This call must not be interleaved with any ghost_exchange calls.
 * \param [in] domain           This domain must remain valid until
 *                              \ref fclaw_domain_indirect_destroy.
 * \return                      A private data structure that will hold
 *                              the context for indirect ghost neighbors.
 */
fclaw_domain_indirect_t
 * fclaw_domain_indirect_begin (fclaw_domain_t * domain);

/** End sending messages to determine neighbors of ghost patches.
 * This call must not be interleaved with any ghost_exchange calls.
 * When this function returns, the necessary information is complete
 * and \ref fclaw_domain_indirect_neighbors may be called any number of times.
 * \param [in] domain           Must be the same domain used in the begin call.
 * \param [in,out] ind          Must be returned by an earlier call to
 *                              \ref fclaw_domain_indirect_begin
 *                              and will be completed with parallel information.
 */
void fclaw_domain_indirect_end (fclaw_domain_t * domain, 
                                fclaw_domain_indirect_t * ind);

/** Call this analogously to \ref fclaw_domain_face_neighbors.
 * Return an indirect face neighbor patch:  It is defined as a ghost patch
 * that is face neighbor to the calling ghost patch and belongs to a process
 * that is neither the owner of that ghost patch nor our own process.
 * \param [in] domain           Must be the same domain used in begin and end.
 * \param [in] ind              Must have been initialized by \ref
 *                              fclaw_domain_indirect_end.
 * \param [in] ghostno          Number of the ghost patch whose neighbors we seek.
 * \param [in] faceno           Number of the ghost patch's face to look across.
 * \param [out] rproc           Processor number of neighbor patches.  Exception 1:
 *                              If the neighbor is a bigger patch, rproc[1] contains
 *                              the number of the small patch as one the smaller faces,
 *                              but only if the neighbor fits the above criteria.
 *                              Exception 2: For non-indirect patches, set it to -1.
 * \param [out] rblockno        The number of the neighbor block.
 * \param [out] rpatchno        Only for indirect ghost patches, we store
 *                              the number relative to our ghost patch array.
 *                              For all other patches, this is -1.
 * \param [out] faceno          The face number and orientation of the neighbor(s).
 * \return                      Only for indirect ghost patches, the size of the
 *                              neighbor(s).  For all others, we set this to
 *                              \ref fclaw_PATCH_BOUNDARY.
 */
fclaw_patch_relation_t
fclaw_domain_indirect_neighbors (fclaw_domain_t * domain,
                                 fclaw_domain_indirect_t * ind,
                                 int ghostno, int faceno, int rproc[],
                                 int *rblockno, int rpatchno[],
                                 int *rfaceno);

/** Destroy all context data for indirect ghost neighbor patches.
 * \param [in] domain           Must be the same domain used in begin and end.
 * \param [in,out] ind          Memory will be freed.
 */
void fclaw_domain_indirect_destroy (fclaw_domain_t * domain,
                                    fclaw_domain_indirect_t * ind);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Meta Domains
/* ---------------------------------------------------------------------- */
///@{

/** Return true if \a domain is an artifical domain.
 *
 * This function can be used in \ref fclaw_interpolate_point_t callbacks to
 * distinguish domains that were created during a partition search (and only
 * contain some meta information) from real domains in a local search.
 */
int fclaw_domain_is_meta (fclaw_domain_t * domain);

///@}
#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FORESTCLAW_H */
