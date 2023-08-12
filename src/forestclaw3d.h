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

#ifndef FORESTCLAW3D_H
#define FORESTCLAW3D_H

#include <fclaw_base.h>
#include <forestclaw.h>
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
 * Main ForestClaw structures and routines
 */

/* ---------------------------------------------------------------------- */
///                      @name Patch Neighbors
/* ---------------------------------------------------------------------- */
///@{

/** Determine physical boundary status as 1, or 0 for neighbor patches.
 * This must ONLY be called for local patches.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in,out] boundaries	Domain boundary boolean flags.
 *			The order is left, right, front, back, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int fclaw3d_patch_boundary_type (fclaw_domain_t * domain,
                                 int blockno, int patchno, int boundaries[6]);

/** Determine whether the normal to a face neighbor align.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in] faceno   Number of the face of the patch.
 * \return		True if normals match, false for mismatch.
 */
int fclaw3d_patch_normal_match (fclaw_domain_t * domain,
                                int blockno, int patchno, int faceno);


/** Determine neighbor patch(es) and orientation across a given face.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] faceno   Number of the patch face: left, right, front, back, bottom, top.
 * \param [out] rproc   Processor number of neighbor patches.  Exception:
 *                      If the neighbor is a bigger patch, rproc[1] contains
 *                      the number of the small patch as one of four half faces.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch numbers for up to 4 neighbors.
 *                              The patch number is relative to its block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghost_patches[.
 * \param [out] rfaceno Neighbor face number and orientation.
 * \return              The relative patch size of the face neighbor.
 */
fclaw_patch_relation_t fclaw3d_patch_face_neighbors (fclaw_domain_t *
                                                       domain, int blockno,
                                                       int patchno,
                                                       int faceno,
                                                       int rproc[4],
                                                       int *rblockno,
                                                       int rpatchno[4],
                                                       int *rfaceno);

/** Change perspective across a face neighbor situation.
 * \param [in,out] faceno   On input, valid face number for a patch.
 *                          On output, valid face number seen from
 *                          faceno's neighbor patch.
 * \param [in,out] rfaceno  On input, encoded neighbor face number as returned
 *                          by fclaw3d_patch_face_neighbors.
 *                          On output, encoded neighbor face number seen from
 *                          faceno's neighbor patch.
 */
void fclaw3d_patch_face_swap (int *faceno, int *rfaceno);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  faceno      The number of the originating face.
 * \param [in]  rfaceno     Encoded as rfaceno = r * 6 + nf, where nf = 0..5 is
 *                          the neigbbor's connecting face number and r = 0..3
 *                          is the relative orientation to the neighbor's face.
 * \param [out] ftransform  This array holds 9 integers.
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
void fclaw3d_patch_face_transformation (int faceno, int rfaceno,
                                        int ftransform[]);

/** Modify an existing face transformation depending on intra-block usage.
 * This function can be called any number of times on the same transform array.
 * \param [in,out] ftransform   Array of values necessarily created by \ref
 *                              fclaw3d_patch_face_transformation.
 * \param [in] sameblock        Transformation supposed to work in same block?
 */
void fclaw3d_patch_face_transformation_block (int ftransform[],
                                              int sameblock);

/** Fill an array with the axis combination of a face neighbor transformation
 * that operates on two patches in the same block (the trivial case).
 * Use when there is no prior call to \ref fclaw3d_patch_face_transformation.
 * Don't postprocess the result any further -- it's only useful intra-block.
 * \param [out] ftransform      Gets initialized to a same-block transform.
 */
void fclaw3d_patch_face_transformation_intra (int ftransform[]);

/** Return whether a face transformation is valid.
 * \param [in] ftransform       Array of values as created by \ref
 *                              fclaw3d_patch_face_transformation,
 *                              possibly modified by \ref
 *                              fclaw3d_patch_face_transformation_block.
 * \return                      True if valid, false if not.
 */
int fclaw3d_patch_face_transformation_valid (const int ftransform[]);

/** Transform a patch coordinate into a neighbor patch's coordinate system.
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
void fclaw3d_patch_transform_face (fclaw_patch_t * ipatch,
                                   fclaw_patch_t * opatch,
                                   const int ftransform[],
                                   int mx, int my, int mz, int based,
                                   int *i, int *j, int *k);

/** Transform a patch coordinate into a neighbor patch's coordinate system.
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
void fclaw3d_patch_transform_face2 (fclaw_patch_t * ipatch,
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
int fclaw3d_patch_edge_neighbors (fclaw_domain_t * domain,
                                  int blockno, int patchno, int edgeno,
                                  int *rproc, int *rblockno, int *rpatchno,
                                  int *redge,
                                  fclaw_patch_relation_t * neighbor_size);

/** Change perspective across an edge neighbor situation.
 * \param [in,out] edgeno       On input, valid edge number for a patch.
 *                              On output, edge number seen from
 *                              the edge neighbor patch.
 * \param [in,out] redgeno      On input, valid edge number as returned
 *                              by fclaw3d_patch_edge_neighbors.
 *                              On output, edge number seen from
 *                              the edge neighbor patch.
 */
void fclaw3d_patch_edge_swap (int *edgeno, int *redgeno);

/** Determine neighbor patch(es) and orientation across a given corner.
 * The current version only supports one neighbor, i.e., no true multi-block.
 * A query across a corner in the middle of a longer face returns the boundary.
 * We only return corner neighbors that are not already face or edge neighbors.
 * Inter-tree corners are only returned if the number of meeting corners is
 * exactly eight.  Nine or more are currently not supported.
 * This must ONLY be called for local patches.
 * \param [in] domain   Valid domain structure.
 * \param [in] blockno  Number of the block within the domain.
 * \param [in] patchno  Number of the patch within the block.
 * \param [in] cornerno	Number of the patch corner: 0=bfl, ..., 7=tbr.
 * \param [out] rproc   Processor number of neighbor patch.
 * \param [out] rblockno        Neighbor block number.
 * \param [out] rpatchno        Neighbor patch number relative to the block.
 *                              If the neighbor is off-processor, this is not
 *                              a patch number but in [0, num_ghosts_patches[.
 * \param [out] rcorner         Number of the corner from the other neighbor.
 * \param [out] neighbor_size   The relative patch size of the neighbor.
 * \return                      True if at least one corner neighbor exists
 *                              that is not already a face neighbor.
 */
int fclaw3d_patch_corner_neighbors (fclaw_domain_t * domain,
                                    int blockno, int patchno, int cornerno,
                                    int *rproc, int *rblockno, int *rpatchno,
                                    int *rcorner,
                                    fclaw_patch_relation_t * neighbor_size);

/** Change perspective across a corner neighbor situation.
 * \param [in,out] cornerno     On input, valid corner number for a patch.
 *                              On output, corner number seen from
 *                              the corner neighbor patch.
 * \param [in,out] rcornerno    On input, valid corner number as returned
 *                              by fclaw3d_patch_corner_neighbors.
 *                              On output, corner number seen from
 *                              the corner neighbor patch.
 */
void fclaw3d_patch_corner_swap (int *cornerno, int *rcornerno);

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
void fclaw3d_patch_transform_corner (fclaw_patch_t * ipatch,
                                     fclaw_patch_t * opatch,
                                     int icorner, int is_block_boundary,
                                     int mx, int my, int mz,
                                     int based, int *i, int *j, int *k);

/** Transform a patch coordinate into a neighbor patch's coordinate system.
 * This function assumes that the neighbor patch is smaller (HALF size) and that
 * the patches lie in coordinate systems with the same orientation.
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
void fclaw3d_patch_transform_corner2 (fclaw_patch_t * ipatch,
                                      fclaw_patch_t * opatch,
                                      int icorner, int is_block_boundary,
                                      int mx, int my, int mz, int based,
                                      int i[], int j[], int k[]);

///@}
/* ---------------------------------------------------------------------- */
///                         @name Exchange
/* ---------------------------------------------------------------------- */
///@{

/** Data structure for storing allocated data for parallel exchange. */
typedef struct fclaw3d_domain_exchange
{
    size_t data_size; /**< The number of bytes per patch to exchange */

    /* These two members are for consistency checking */
    int num_exchange_patches; /**< Number of patches to send */
    int num_ghost_patches; /**< Number of patches to receive */
    /**
       One pointer per processor-local exchange patch in order, for a total
       count of domain->num_exchange_patches.  This applies precisely to local
       patches that touch the parallel boundary from the inside, i.e., if
       (flags & FCLAW3D_PATCH_ON_PARALLEL_BOUNDARY).
     */
    void **patch_data;
    /**
       Array of domain->num_ghost_patches many void pointers, each pointing to
       exactly data_size bytes of memory that can be read from by forestclaw
       after each fclaw3d_domain_parallel_exchange.
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
fclaw3d_domain_exchange_t;

/** Allocate buffer to hold the data from off-processor patches.
 * Free this by fclaw3d_domain_free_after_exchange before regridding.
 * \param [in] domain           The domain is not modified.
 * \param [in] data_size        Number of bytes per patch to exchange.
 * \return                      Allocated data structure.
 *                              The pointers in patch_data[i] need to be set
 *                              after this call by forestclaw.
 */
fclaw3d_domain_exchange_t
    * fclaw3d_domain_allocate_before_exchange (fclaw_domain_t * domain,
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
void fclaw3d_domain_ghost_exchange (fclaw_domain_t * domain,
                                    fclaw3d_domain_exchange_t * e,
                                    int exchange_minlevel,
                                    int exchange_maxlevel);

/** Start asynchronous exchange of parallel ghost neighbors.
 * The arguments are the same as for fclaw3d_domain_ghost_exchange.
 * It must be followed by a call to fclaw3d_domain_ghost_exchange_end.
 * Between begin and end, neither of \ref fclaw3d_domain_indirect_begin
 * and _end must be called.
 * \param [in,out] e            Its ghost_data member must survive and not
 *                              be written to until the completion of
 *                              fclaw3d_domain_ghost_exchange_end.
 *                              Its patch_data member may already be
 *                              overwritten after this function returns.
 */
void fclaw3d_domain_ghost_exchange_begin (fclaw_domain_t * domain,
                                          fclaw3d_domain_exchange_t * e,
                                          int exchange_minlevel,
                                          int exchange_maxlevel);

/** Complete asynchronous exchange of parallel ghost neighbors.
 * Must be called at some point after fclaw3d_domain_ghost_exchange_begin.
 * Between begin and end, neither of \ref fclaw3d_domain_indirect_begin
 * and _end must be called.
 * \param [in,out] e            Its ghost_data member must have survived.
 */
void fclaw3d_domain_ghost_exchange_end (fclaw_domain_t * domain,
                                        fclaw3d_domain_exchange_t * e);

/** Free buffers used in exchanging off-processor data during time stepping.
 * This should be done just before regridding.
 * \param [in] domain           The domain is not modified.
 * \param [in] e                Allocated buffers.
 */
void fclaw3d_domain_free_after_exchange (fclaw_domain_t * domain,
                                         fclaw3d_domain_exchange_t * e);

///@}
/* ---------------------------------------------------------------------- */
///                 @name Indirect Parallel Neighbors
/* ---------------------------------------------------------------------- */
///@{

/* Data structure for storing indirect parallel neighbor information */
typedef struct fclaw3d_domain_indirect fclaw3d_domain_indirect_t;

/** Begin sending messages to determine neighbors of ghost patches.
 * This call must not be interleaved with any ghost_exchange calls.
 * \param [in] domain           This domain must remain valid until
 *                              \ref fclaw3d_domain_indirect_destroy.
 * \return                      A private data structure that will hold
 *                              the context for indirect ghost neighbors.
 */
fclaw3d_domain_indirect_t
    * fclaw3d_domain_indirect_begin (fclaw_domain_t * domain);

/** End sending messages to determine neighbors of ghost patches.
 * This call must not be interleaved with any ghost_exchange calls.
 * When this function returns, the necessary information is complete
 * and \ref fclaw3d_domain_indirect_neighbors may be called any number of times.
 * \param [in] domain           Must be the same domain used in the begin call.
 * \param [in,out] ind          Must be returned by an earlier call to
 *                              \ref fclaw3d_domain_indirect_begin
 *                              and will be completed with parallel information.
 */
void fclaw3d_domain_indirect_end (fclaw_domain_t * domain,
                                  fclaw3d_domain_indirect_t * ind);

/** Call this analogously to \ref fclaw3d_domain_face_neighbors.
 * We only return an indirect ghost neighbor patch:  This is defined as a ghost
 * patch that is neighbor to the calling ghost patch and belongs to a processor
 * that is neither the owner of that ghost patch nor our own processor.
 * \param [in] domain           Must be the same domain used in begin and end.
 * \param [in] ind              Must have been initialized by \ref
 *                              fclaw3d_domain_indirect_end.
 * \param [in] ghostno          Number of the ghost patch whose neighbors we seek.
 * \param [in] faceno           Number of the ghost patch's face to look across.
 * \param [out] rproc           Processor number of neighbor patches.  Exception 1:
 *                              If the neighbor is a bigger patch, rproc[1] contains
 *                              the number of the small patch as one of two half faces,
 *                              but only if the neighbor fits the above criteria.
 *                              Exception 2: For non-indirect patches, set it to -1.
 * \param [out] rblockno        The number of the neighbor block.
 * \param [out] rpatchno        Only for indirect ghost patches, we store
 *                              the number relative to our ghost patch array.
 *                              For all other patches, this is -1.
 * \param [out] faceno          The face number and orientation of the neighbor(s).
 * \return                      Only for indirect ghost patches, the size of the
 *                              neighbor(s).  For all others, we set this to
 *                              \ref FCLAW3D_PATCH_BOUNDARY.
 */
fclaw_patch_relation_t
fclaw3d_domain_indirect_neighbors (fclaw_domain_t * domain,
                                   fclaw3d_domain_indirect_t * ind,
                                   int ghostno, int faceno, int rproc[2],
                                   int *rblockno, int rpatchno[2],
                                   int *rfaceno);

/** Destroy all context data for indirect ghost neighbor patches.
 * \param [in] domain           Must be the same domain used in begin and end.
 * \param [in,out] ind          Memory will be freed.
 */
void fclaw3d_domain_indirect_destroy (fclaw_domain_t * domain,
                                      fclaw3d_domain_indirect_t * ind);

///@}
/* ---------------------------------------------------------------------- */
///                     @name Communication
/* ---------------------------------------------------------------------- */
///@{

/** Compute and return the maximum over all processors of a double value.
 * The minimum can be computed by using this function on the negative value.
 */
double fclaw3d_domain_global_maximum (fclaw_domain_t * domain, double d);

/** Compute and return the sum over all processors of a double value.
 */
double fclaw3d_domain_global_sum (fclaw_domain_t * domain, double d);

/** Synchronize all processes.  Avoid using if at all possible.
 */
void fclaw3d_domain_barrier (fclaw_domain_t * domain);

/** Serialize a section of code.
 * THIS IS NOT SCALABLE.
 * WILL BE HORRIBLY SLOW FOR LARGE NUMBERS OF PROCESSORS.
 * A processor returns from this function only after all lower-numbered
 * processors have called fclaw3d_domain_serialization_leave.
 * No collective communication routines must be called between the calls
 * to this function and fclaw3d_domain_serialization_leave.
 * \param [in] domain           The domain is not modified.
 */
void fclaw3d_domain_serialization_enter (fclaw_domain_t * domain);

/** Serialize a section of code.
 * THIS IS NOT SCALABLE.
 * WILL BE HORRIBLY SLOW FOR LARGE NUMBERS OF PROCESSORS.
 * A processor must call this function to allow all higher-numbered
 * processors to return from fclaw3d_domain_serialization_enter.
 * \param [in] domain           The domain is not modified.
 */
void fclaw3d_domain_serialization_leave (fclaw_domain_t * domain);

///@}
/* ---------------------------------------------------------------------- */
///                      @name Meta Domains
/* ---------------------------------------------------------------------- */
///@{

/** Return true if \a domain is an artifical domain.
 *
 * This function can be used in \ref fclaw3d_interpolate_point_t callbacks to
 * distinguish domains that were created during a partition search (and only
 * contain some meta information) from real domains in a local search.
 */
int fclaw3d_domain_is_meta (fclaw_domain_t * domain);

/** Initialize a meta domain.
 *
 * Initializes \a domain in an artificial manner, where the entry mpirank is
 * used to store arbitrary context information. The remaining entries are
 * initialized to -1 or NULL.
 * The resulting domain can be passed to an \ref fclaw3d_interpolate_point_t
 * in case the domain to interpolate on is not available locally (also see
 * \ref fclaw3d_overlap_exchange for an example).
 */
void fclaw3d_domain_init_meta (fclaw_domain_t *domain, int mpirank);

///@}
#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FORESTCLAW3D_H */
