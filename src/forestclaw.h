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

#ifndef FCLAW_PATCH_H
#define FCLAW_PATCH_H

/*
 * Domain-independent header file should not include domain-specific headers.
 * The corresponding source file include the 2d and 3d domain-specific headers.
 */
#include <fclaw_base.h>

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

#endif /* !FCLAW_PATCH_H */
