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
/**
 * @file
 * Dimension-independent wrapper of a forestclaw patch.
 */

#ifndef FORESTCLAW_H
#define FORESTCLAW_H

/*
 * Domain-independent header file should not include domain-specific headers.
 * The corresponding source file include the 2d and 3d domain-specific headers.
 */
#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* ---------------------------------------------------------------------- */
///                      @name Data Types
/* ---------------------------------------------------------------------- */
///@{

/**
 * @brief Enum for encoding patch information
 */
typedef enum
{
    /** Number relative to parent */
    FCLAW2D_PATCH_CHILDID = 0x7,
    /** Patch is the first sibling */
    FCLAW2D_PATCH_FIRST_SIBLING = 0x8,
    /** Has neighbor on different processor */
    FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY = 0x10,
    /** Patch is a ghost patch */
    FCLAW2D_PATCH_IS_GHOST = 0x20,
    /** Face 0 is on a block boundary */
    FCLAW2D_PATCH_ON_BLOCK_FACE_0 = 0x040,
    /** Face 1 is on a block boundary */
    FCLAW2D_PATCH_ON_BLOCK_FACE_1 = 0x080,
    /** Face 2 is on a block boundary */
    FCLAW2D_PATCH_ON_BLOCK_FACE_2 = 0x100,
    /** Face 3 is on a block boundary */
    FCLAW2D_PATCH_ON_BLOCK_FACE_3 = 0x200,
    /** Patch is on a block boundary */
    FCLAW2D_PATCH_ON_BLOCK_BOUNDARY = 0xFC0
}
fclaw2d_patch_flags_t;

/** For each of the four faces, the corresponding block boundary flag. */
extern const fclaw2d_patch_flags_t fclaw2d_patch_block_face_flags[4];

typedef enum
{
    /** Number relative to parent */
    FCLAW3D_PATCH_CHILDID = 0x7,
    /** Patch is the first sibling */
    FCLAW3D_PATCH_FIRST_SIBLING = 0x8,
    /** Has neighbor on different processor */
    FCLAW3D_PATCH_ON_PARALLEL_BOUNDARY = 0x10,
    /** Patch is a ghost patch */
    FCLAW3D_PATCH_IS_GHOST = 0x20,
    /** Face 0 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_0 = 0x040,
    /** Face 1 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_1 = 0x080,
    /** Face 2 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_2 = 0x100,
    /** Face 3 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_3 = 0x200,
    /** Face 4 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_4 = 0x400,
    /** Face 5 is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_FACE_5 = 0x800,
    /** Patch is on a block boundary */
    FCLAW3D_PATCH_ON_BLOCK_BOUNDARY = 0xFC0
}
fclaw3d_patch_flags_t;

/** For each of the 6 faces, the corresponding block boundary flag. */
extern const fclaw3d_patch_flags_t fclaw3d_patch_block_face_flags[6];


typedef struct fclaw_patch fclaw_patch_t;
typedef struct fclaw_patch_bounds_2d fclaw_patch_bounds_2d_t;
typedef struct fclaw_patch_bounds_3d fclaw_patch_bounds_3d_t;

struct fclaw_patch_bounds_2d
{
    double xlower, xupper;
    double ylower, yupper;
};

struct fclaw_patch_bounds_3d
{
    double xlower, xupper;
    double ylower, yupper;
    double zlower, zupper;
};

/** 
 * @brief The metadata structure for a forest leaf, which is a forestclaw patch.
 * The patch may be either a process-local patch or a ghost patch.
 */
struct fclaw_patch
{
    int dim;
    fclaw_patch_bounds_2d_t* d2;
    fclaw_patch_bounds_3d_t* d3;
    int level;                  /**< 0 is root, increases if refined */
    int target_level;           /**< level desired after adaptation */
    int flags;                  /**< flags that encode tree information */
    /** Union, If this is a local patch, it points to the next patch, otherwise it gives
     * the bock number of this patch */
    union
    {
        fclaw_patch_t *next;  /**< local: next patch same level same block */
        int blockno;            /**< off-proc: this patch's block number */
    }
    u;
    void *user;                 /**< User Pointer */
};

typedef struct fclaw_block_d2
{
    /** @{ @brief lower left coordinate */
    double xlower, xupper;
    /** @} */
    /** @{ @brief upper right coordinate */
    double ylower, yupper;
    /** @} */
    double vertices[4 * 3];     /**< for each block corner, the xyz coordinates
                                     of the p4est_connectivity structure */
    int is_boundary[4];         /**< physical boundary flag */
} fclaw_block_d2_t;

typedef struct fclaw_block_d3
{
    /** @{ @brief left/right coordinate */
    double xlower, xupper;
    /** @} */
    /** @{ @brief front/back coordinate */
    double ylower, yupper;
    /** @} */
    /** @{ @brief bottom/top coordinate */
    double zlower, zupper;
    /** @} */
    double vertices[8 * 3];     /**< for each block corner, the xyz coordinates
                                     of the p8est_connectivity structure */
    int is_boundary[6];         /**< physical boundary flag */
} fclaw_block_d3_t;

/**
 * @brief Data Structure for a block
 */
typedef struct fclaw_block
{
    int dim;
    fclaw_block_d2_t* d2;
    fclaw_block_d3_t* d3;
    int num_patches;            /**< local patches in this block */
    int num_patches_before;     /**< in all previous blocks */
    int num_exchange_patches;   /**< exchange patches in this block */
    /** @{ 
     * @brief min/max level
     * local over this block.  If this proc doesn't
     * store any patches in this block, we set
     * maxlevel < 0 <= minlevel. 
     */
    int minlevel;
    int maxlevel;
    /** @} */
    fclaw_patch_t *patches;           /**< The patches for this block */
    fclaw_patch_t **patchbylevel;     /**< Pointer to the first patch in each level **/
    fclaw_patch_t **exchange_patches; /**< Pointer for each exchange patch */
    void *user;                         /**< User pointer */
} fclaw_block_t;

/** This structure identify parameters that are copied from a domain
 * to a new domain derived by adaptation or partitioning. */
typedef struct fclaw_domain_persist
{
    int smooth_refine;          /**< Boolean tells us whether to communicate
                                     the desired refinement level to neighbors. */
    int smooth_level;           /**< The minimum level that refinement smoothing
                                     is enabled on.  Use 0 for al levels. */
}
fclaw_domain_persist_t;

typedef struct fclaw_domain_exchange fclaw_domain_exchange_t;
typedef struct fclaw2d_domain_indirect fclaw2d_domain_indirect_t;

typedef struct fclaw2d_domain_data
{
    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;

    fclaw_domain_exchange_t *domain_exchange;
    fclaw2d_domain_indirect_t *domain_indirect;

} fclaw2d_domain_data_t;

typedef struct fclaw3d_domain_indirect fclaw3d_domain_indirect_t;

typedef struct fclaw3d_domain_data
{
    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;

    fclaw_domain_exchange_t *domain_exchange;
    fclaw3d_domain_indirect_t *domain_indirect;

} fclaw3d_domain_data_t;


/**
 * @brief The domain structure is a collection of blocks
 * 
 * The domain structure gives a processor local view of the grid hierarchy.
 * Unless explicitly noted otherwise, all variables are processor local,
 * i.e., they are generally different on each processor.
 * Variables that are synchronized and shared between processors
 * are denoted *global*.
 */
typedef struct fclaw_domain
{
    int dim;
    fclaw2d_domain_data_t* d2;
    fclaw3d_domain_data_t* d3;

    sc_MPI_Comm mpicomm;        /**< MPI communicator */
    int mpisize;                /**< MPI size */
    int mpirank;                /**< MPI rank */
    int possible_maxlevel;      /**< theoretical maximum that can be reached */

    fclaw_domain_persist_t p;         /**< Parameters that carry over from
                                             one domain to a derived one. */

    int local_num_patches;      /**< sum of patches over all blocks on this proc */
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

    int just_adapted;           /**< true after non-trivial adaptation */
    int just_partitioned;       /**< true after non-trivial partition */
    int partition_unchanged_first;      /**< local index of first unchanged patch */
    int partition_unchanged_length;     /**< number of unchanged quadrants */
    int partition_unchanged_old_first;  /**< local index wrt. previous partition */

    int num_blocks;             /**< Total number of blocks. */
    fclaw_block_t *blocks;    /**< allocated storage */
    int num_exchange_patches;   /**< number my patches relevant to other procs.
                                   Identified by this expression to be true:
                                   (patch->flags &
                                   FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY) */
    fclaw_patch_t **exchange_patches; /**< explicitly store exchange patches */
    int num_ghost_patches;      /**< number of off-proc patches relevant to this proc */
    fclaw_patch_t *ghost_patches;     /**< array of off-proc patches */

    void **mirror_target_levels;  /**< Points to target level of each mirror. */
    int *ghost_target_levels;   /**< Contains target level for each ghost. */

    void *pp;                   /**< opaque backend data */
    int pp_owned;               /**< True if the pp member is owned by this domain */
    sc_keyvalue_t *attributes;  /**< Reserved to store domain attributes */

    void *user; /**< user data pointer */
} fclaw_domain_t;

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

///@}

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

/** Callback prototype used in fclaw3d_domain_iterate_adapted.
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
 *                              fclaw3d_domain_iterate_partitioned.
 */
typedef void (*fclaw_transfer_callback_t) (fclaw_domain_t * old_domain,
                                             fclaw_patch_t * old_patch,
                                             fclaw_domain_t * new_domain,
                                             fclaw_patch_t * new_patch,
                                             int blockno,
                                             int old_patchno, int new_patchno,
                                             void *user);
///* ---------------------------------------------------------------------- */
///                      @name Domain Attributes
/* ---------------------------------------------------------------------- */
///@{

/** Add a named attribute to the domain.
 * Attribute names starting with 'fclaw' are reserved.
 * \param [in] domain   This domain will get a new attribute.
 * \param [in] name     This name must not yet be used for another attribute.
 * \param [in] attribute        Arbitrary data stored under \a name.
 */
void fclaw_domain_attribute_add (fclaw_domain_t * domain,
                                 const char *name, void *attribute);

/** Access a named attribute of the domain.
 * \param [in] domain   The domain may or may not have the queried attribute.
 * \param [in] name     The attribute by this \a name is retrieved.
 * \param [in] default_attr     Returned if the attribute does not exist.
 * \return              The data that was previously stored under \a name,
 *                      or \a default_attr if the attribute does not exist.
 */
void *fclaw_domain_attribute_access (fclaw_domain_t * domain,
                                     const char *name, void *default_attr);

/** Remove a named attribute from the domain.
 * It is NOT necessary to call this function before domain destruction.
 * \param [in] domain   The domain must have the attribute \a name.
 * \param [in] name     An attribute of this name must exist.
 */
void fclaw_domain_attribute_remove (fclaw_domain_t * domain,
                                    const char *name);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Patch Functions
/* ---------------------------------------------------------------------- */
///@{

/** Return the dimension of a corner.
 * This function is LEGAL to call for both local and ghost patches.
 * \param [in] patch    A patch with properly set member variables.
 * \param [in] cornerno A corner number in 0..3.
 * \return              0 if the corner is always at a fourfold intersection,
 *                      1 if the corner would end up in the middle of a face
 *                      when there is a coarser neighbor.
 */
int fclaw_patch_corner_dimension (const fclaw_patch_t * patch,
                                  int cornerno);

/** Return the number of a patch with respect to its parent in the tree.
 * This function is LEGAL to call for both local and ghost patches.
 * \param [in] patch    A patch with properly set member variables.
 * \return              The child id is a number in 0..3.
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
///                      @name Patch Neighbors
/* ---------------------------------------------------------------------- */
///@{

/** Change perspective across a face neighbor situation.
 * \param [in]     dim      The space dimension.
 * \param [in,out] faceno   On input, valid face number for a patch.
 *                          On output, valid face number seen from
 *                          faceno's neighbor patch.
 * \param [in,out] rfaceno  On input, encoded neighbor face number as returned
 *                          by fclaw2d_patch_face_neighbors.
 *                          On output, encoded neighbor face number seen from
 *                          faceno's neighbor patch.
 */
void fclaw_patch_face_swap (int dim, int *faceno, int *rfaceno);

/** Determine physical boundary status as 1, or 0 for neighbor patches.
 * This must ONLY be called for local patches.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in,out] boundaries	Domain boundary boolean flags.
 *			The order is left, right, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int fclaw_patch_boundary_type (fclaw_domain_t * domain,
                               int blockno, int patchno, int boundaries[4]);

/** Determine whether the normal to a face neighbor align.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in] faceno   Number of the face of the patch.
 * \return		True if normals match, false for mismatch.
 */
int fclaw_patch_normal_match (fclaw_domain_t * domain,
                              int blockno, int patchno, int faceno);

/** Determine neighbor patch(es) and orientation across a given face.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] faceno   Number of the patch face: left, right, bottom, top.
 * \param [out] rproc   Processor number of neighbor patches.  Exception:
 *                      If the neighbor is a bigger patch, rproc[1] contains
 *                      the number of the small patch as one of two half faces.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch numbers for up to 2 neighbors.
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
                                                   int rproc[2],
                                                   int *rblockno,
                                                   int rpatchno[2],
                                                   int *rfaceno);

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
 * \param [in] cornerno	Number of the patch corner: 0=bl, 1=br, 2=tl, 3=tr.
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

///@}
/* ---------------------------------------------------------------------- */
///                      @name Patch Iterators
/* ---------------------------------------------------------------------- */
///@{

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
///                  @name Topological Properties
/* ---------------------------------------------------------------------- */
///@{


/** Return the space dimension. */
int fclaw_domain_dimension (const fclaw_domain_t * domain);

/** Return refine factor in tree: 2 in 2D, 4 in 3D. */
int fclaw_domain_refine_factor (const fclaw_domain_t * domain);

/** Return number of children in tree: 4 in 2D, 8 in 3D. */
int fclaw_domain_num_children (const fclaw_domain_t * domain);

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

/** Find the numbers of faces adjacent to a cube corner: 2 in 2D, 3 in 3D. */
void fclaw_domain_corner_faces (const fclaw_domain_t * domain,
                                int icorner, int * faces);

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
 * processors have called fclaw2d_domain_serialization_leave.
 * No collective communication routines must be called between the calls
 * to this function and fclaw2d_domain_serialization_leave.
 * \param [in] domain           The domain is not modified.
 */
void fclaw_domain_serialization_enter (fclaw_domain_t * domain);

/** Serialize a section of code.
 * THIS IS NOT SCALABLE.
 * WILL BE HORRIBLY SLOW FOR LARGE NUMBERS OF PROCESSORS.
 * A processor must call this function to allow all higher-numbered
 * processors to return from fclaw2d_domain_serialization_enter.
 * \param [in] domain           The domain is not modified.
 */
void fclaw_domain_serialization_leave (fclaw_domain_t * domain);

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
 * \ref fclaw2d_domain_iterate_adapted.
 */
void fclaw_patch_mark_refine (fclaw_domain_t * domain,
                              int blockno, int patchno);

/** Mark a patch for coarsening.
 * This must ONLY be called for local patches.
 * Coarsening will only happen if the patch family is not further refined
 * and all sibling patches are marked as well.
 * It is safe to call this function from an iterator callback except
 * \ref fclaw2d_domain_iterate_adapted.
 */
void fclaw_patch_mark_coarsen (fclaw_domain_t * domain,
                               int blockno, int patchno);

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
///                         @name Exchange
/* ---------------------------------------------------------------------- */
///@{

/** Data structure for storing allocated data for parallel exchange. */
typedef struct fclaw_domain_exchange
{
    size_t data_size; /**< The number of bytes per patch to exchange */

    /* These two members are for consistency checking */
    int num_exchange_patches; /**< Number of patches to send */
    int num_ghost_patches; /**< Number of patches to receive */
    /**
       One pointer per processor-local exchange patch in order, for a total
       count of domain->num_exchange_patches.  This applies precisely to local
       patches that touch the parallel boundary from the inside, i.e., if
       (flags & FCLAW2D_PATCH_ON_PARALLEL_BOUNDARY).
     */
    void **patch_data;
    /**
       Array of domain->num_ghost_patches many void pointers, each pointing to
       exactly data_size bytes of memory that can be read from by forestclaw
       after each fclaw2d_domain_parallel_exchange.
     */
    void **ghost_data;
    /**
       Memory where the ghost patch data actually lives.
       The above array ghost_data points into this memory.
       It will not be necessary to dereference this memory explicitly.
     */
    char *ghost_contiguous_memory;

    /** Temporary storage required for asynchronous ghost exchange.
     * It is allocated and freed by the begin/end calls below.
     */
    void *async_state;
    int inside_async;           /**< Between asynchronous begin and end? */
    int by_levels;              /**< Did we use levels on the inside? */
}
fclaw_domain_exchange_t;

/** Allocate buffer to hold the data from off-processor patches.
 * Free this by fclaw2d_domain_free_after_exchange before regridding.
 * \param [in] domain           The domain is not modified.
 * \param [in] data_size        Number of bytes per patch to exchange.
 * \return                      Allocated data structure.
 *                              The pointers in patch_data[i] need to be set
 *                              after this call by forestclaw.
 */
fclaw_domain_exchange_t
    * fclaw_domain_allocate_before_exchange (fclaw_domain_t * domain,
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
 * The arguments are the same as for fclaw2d_domain_ghost_exchange.
 * It must be followed by a call to fclaw2d_domain_ghost_exchange_end.
 * Between begin and end, neither of \ref fclaw2d_domain_indirect_begin
 * and _end must be called.
 * \param [in,out] e            Its ghost_data member must survive and not
 *                              be written to until the completion of
 *                              fclaw2d_domain_ghost_exchange_end.
 *                              Its patch_data member may already be
 *                              overwritten after this function returns.
 */
void fclaw_domain_ghost_exchange_begin (fclaw_domain_t * domain,
                                        fclaw_domain_exchange_t * e,
                                        int exchange_minlevel,
                                        int exchange_maxlevel);

/** Complete asynchronous exchange of parallel ghost neighbors.
 * Must be called at some point after fclaw2d_domain_ghost_exchange_begin.
 * Between begin and end, neither of \ref fclaw2d_domain_indirect_begin
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

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FORESTCLAW_H */
