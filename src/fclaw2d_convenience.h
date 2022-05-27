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

#ifndef FCLAW2D_CONVENIENCE_H
#define FCLAW2D_CONVENIENCE_H

#include <forestclaw2d.h>
#include <fclaw2d_map.h>
#include <p4est_connectivity.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* TODO: The unit square is a special case of the brick.  Use that instead. */
fclaw2d_domain_t *fclaw2d_domain_new_unitsquare (sc_MPI_Comm mpicomm,
                                                 int initial_level);
/* TODO: The torus is a special case of the brick.  Use that instead. */
fclaw2d_domain_t *fclaw2d_domain_new_torus (sc_MPI_Comm mpicomm,
                                            int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_twosphere (sc_MPI_Comm mpicomm,
                                                int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_cubedsphere (sc_MPI_Comm mpicomm,
                                                  int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_disk (sc_MPI_Comm mpicomm,
                                           int initial_level);

/** Create a brick connectivity, that is, a rectangular grid of blocks.
 * The origin is in the lower-left corner of the brick.
 * \param [in] mpicomm          We expect sc_MPI_Init to be called earlier.
 * \param [in] blocks_in_x      Positive number of blocks in x direction.
 * \param [in] blocks_in_y      Positive number of blocks in y direction.
 * \param [in] periodic_in_x    True if the right side of the rightmost blocks
 *                              connect periodically to the left side of the
 *                              leftmost blocks.
 * \param [in] periodic_in_y    Periodicity along the vertical direction.
 * \param [in] initial_level    A non-negative integer <= P4EST_QMAXLEVEL.
 * \param [in] cont             We do NOT take ownership of the mapping.
 *                              It is legal to pass NULL for this parameter.
 * \return                      A fully initialized domain structure.
 */
fclaw2d_domain_t *fclaw2d_domain_new_brick (sc_MPI_Comm mpicomm,
                                            int blocks_in_x, int blocks_in_y,
                                            int periodic_in_x,
                                            int periodic_in_y,
                                            int initial_level,
                                            fclaw2d_map_context_t * cont);

/** Create a domain from a given forest connectivity and matching map.
 * \param [in] mpicomm          We expect sc_MPI_Init to be called earlier.
 * \param [in] initial_level    A non-negative integer <= P4EST_QMAXLEVEL.
 * \param [in] conn             We DO take ownership of the connectivity.
 * \param [in] cont             We do NOT take ownership of the mapping.
 * \return                      A fully initialized domain structure.
 */
fclaw2d_domain_t *fclaw2d_domain_new_conn_map (sc_MPI_Comm mpicomm,
                                               int initial_level,
                                               p4est_connectivity_t * conn,
                                               fclaw2d_map_context_t * cont);

void fclaw2d_domain_destroy (fclaw2d_domain_t * domain);

/** Create a new domain based on refine and coarsen marks set previously.
 * All refine and coarsen markers are cancelled when this function is done.
 * \param [in,out] domain       Current domain with set adaptation markers.
 *                              It stays alive because it is needed to
 *                              project numerical values to the adapted domain.
 *                              If adapted, no queries are allowed afterwards.
 * \return                      Adapted domain if refinement occurred, or NULL.
 *                              The return status is identical across all ranks.
 */
fclaw2d_domain_t *fclaw2d_domain_adapt (fclaw2d_domain_t * domain);

/** Create a repartitioned domain after fclaw2d_domain_adapt returned non-NULL.
 * All refine and coarsen markers are cancelled when this function is done.
 * \param [in,out] domain       Current domain that was adapted previously.
 *                              It stays alive because it is needed to
 *                              transfer numerical values to the new partition.
 *                              If partitioned, no queries allowed afterwards.
 * \param [in] weight_exponent  The patches are weighted with an integer factor
 *                              2 ** (level * exponent).  If the exponent is 0,
 *                              all patches have equal weight.  If it is 1,
 *                              smaller patches are considered more expensive
 *                              by a factor two per level increase.
 * \return                      Partitioned domain if different, or NULL.
 *                              The return status is identical across all ranks.
 */
fclaw2d_domain_t *fclaw2d_domain_partition (fclaw2d_domain_t * domain,
                                            int weight_exponent);

/** Query the window of patches that is not transferred on partition.
 * \param [in] domain           A domain after a non-trivial partition
 *                              and before calling \ref fclaw2d_domain_complete.
 * \param [out] unchanged_first         First still-local patch in the new partition.
 * \param [out] unchanged_length        Number of patches that not changed owners.
 * \param [out] unchanged_old_first     First stayed_local patch in the old partition.
 */
void fclaw2d_domain_partition_unchanged (fclaw2d_domain_t * domain,
                                         int *unchanged_first,
                                         int *unchanged_length,
                                         int *unchanged_old_first);

/** Clean up after fclaw2d_domain_partition returned non-NULL.
 * \param [in,out] domain       Current domain that was partitioned.
 */
void fclaw2d_domain_complete (fclaw2d_domain_t * domain);

/** Write VTK file(s) for a domain structure.
 *  Each patch is drawn as one rectangle.
 *  We ignore any geometric transformations.
 *  and use the vertex locations specified in the p4est's connectivity.
 * \param [in] domain           A valid domain structure.  Is not changed.
 * \param [in] basename         Filename prefix passed to p4est_vtk functions.
 */
void fclaw2d_domain_write_vtk (fclaw2d_domain_t * domain,
                               const char *basename);

/** Print patch number by level on all processors */
void fclaw2d_domain_list_levels (fclaw2d_domain_t * domain, int log_priority);

/** Print face neighbor status for each face */
void fclaw2d_domain_list_neighbors (fclaw2d_domain_t * domain,
                                    int log_priority);

/** Print information on adapted patches */
void fclaw2d_domain_list_adapted (fclaw2d_domain_t * old_domain,
                                  fclaw2d_domain_t * new_domain,
                                  int log_priority);

/** Search triples of (block number, x coordinate, y coordinate) in the mesh.
 * The x, y coordinates must be in [0, 1]^2.
 * The input data must be equal on every process: This is a collective call.
 * The results will also be equal on every process.
 *
 * A point is found correctly even if it is on a patch boundary.
 * We return the smallest patch number on the smallest processor touching it.
 * However, if a point is on a block boundary, it must be decided before
 * calling this function which tree shall be queried for it.
 *
 * \note Currently we do not find the smallest matching process, but instead
 *       instead a point on a parallel boundary may be found on multiple processes.
 *       This should be fixed in the near future.
 *
 * \param [in] domain           Must be valid domain structure.  Will not be changed.
 * \param [in] block_offsets    Array of (num_blocks + 1) int variables.
 *                              The points to search in block t in [0, num_blocks)
 *                              have indices [block_offsets[t], block_offsets[t + 1])
 *                              in the \b coordinates and results arrays.
 * \param [in] coordinates      An array of elem_size == 2 * sizeof (double) with
 *                              entries (x, y) in [0, 1]^2.  Of these entries,
 *                              there are \b block_offsets[num_blocks] many.
 * \param [in,out] results      On input, an array of type int and
 *                              \b block_offsets[num_blocks] many entries.
 *                              On output, each entry will be -1 if the point has
 *                              not been found, or the patch number within its block.
 */
void fclaw2d_domain_search_points (fclaw2d_domain_t * domain,
                                   sc_array_t * block_offsets,
                                   sc_array_t * coordinates,
                                   sc_array_t * results);

/** Callback function to compute the integral of a "ray" within a patch.
 *
 * This function can be passed to \ref fclaw2d_domain_integrate_rays to
 * eventually compute the integrals over the whole domain for an array of rays.
 *
 * \param [in] domain           The domain to integrate on.
 * \param [in] patch            The patch under consideration.
 *                              When on a leaf, this is a valid forestclaw patch.
 *                              Otherwise, this is a temporary artificial patch
 *                              containing all standard patch information except
 *                              for the pointer to the next patch and user-data.
 *                              Only the FCLAW2D_PATCH_CHILDID and the
 *                              FCLAW2D_PATCH_ON_BLOCK_FACE_* flags are set.
 *                              Artificial patches are generally ancestors of
 *                              valid forestclaw patches that are leaves.
 * \param [in] blockno          The block id of the patch under consideration.
 * \param [in] patchno          When on a leaf, this is a valid patch number,
 *                              as always relative to its block.  For a leaf,
 *                              this callback must set the integral value to
 *                              the local contribution of this patch and ray.
 *                              Otherwise, patchno is -1.  In this case, the
 *                              integral value is ignored.
 * \param [in] ray              Representation of a "ray"; user-defined.
 *                              Points to an array element of the rays passed
 *                              to \ref fclaw2d_domain_integrate_rays.
 * \param [in,out] integral     The integral value associated with the ray.
 *                              On input this is 0.
 *                              For leaves this callback must set it to the
 *                              exact integral contribution for this patch and
 *                              ray.
 * \param [in,out] user         Arbitrary data passed in earlier.
 * \return                      True if there is a possible/partial intersection of the
 *                              patch (which may be an ancestor) with the ray.
 *                              This may be a false positive; we'll be fine.
 *                              Return false if there is definitely no intersection.
 *                              Only for leaves, this function must compute
 *                              the exact integral contribution for this
 *                              patch by intersecting this ray and store it in
 *                              the \a integral output argument.
 *                              The integral value may well be 0. if the intersection
 *                              is, in fact, none (a false positive).
 */
typedef int (*fclaw2d_integrate_ray_t) (fclaw2d_domain_t * domain,
                                        fclaw2d_patch_t * patch,
                                        int blockno, int patchno,
                                        void *ray, double *integral,
                                        void *user);

/** Compute the integrals of an array of user-defined rays.
 * The integral for each ray and intersection quadrant is supplied by a callback.
 * We store the results in an array of integral values of type double.
 *
 * \param [in] domain           The domain to integrate on.
 * \param [in] intersect        Callback function that returns true if a ray
 *                              intersects a patch and -- when called for a leaf
 *                              -- shall output the integral of the ray segment.
 * \param [in] rays             Array containing the rays of user-defined type.
 *                              Each entry contains one item of arbitrary data.
 *                              We do not dereference, just pass pointers around.
 * \param [in,out] integrals    Array of double variables.  The number of entries
 *                              must equal the number of rays.  Input values ignored.
 *                              On output, we provide the final integral values.
 * \param [in,out] user         Arbitrary data to be passed to the callback.
 */
void fclaw2d_domain_integrate_rays (fclaw2d_domain_t * domain,
                                    fclaw2d_integrate_ray_t intersect,
                                    sc_array_t * rays,
                                    sc_array_t * integrals,
                                    void * user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CONVENIENCE_H */
