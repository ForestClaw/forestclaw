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
#include <p4est_extended.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#define FCLAW2D_FILE_USER_STRING_BYTES P4EST_FILE_USER_STRING_BYTES

/** Opaque context used for writing a fclaw2d data file. */
typedef struct fclaw2d_file_context fclaw2d_file_context_t;

/** Questions: Do we want to fix the glob for each writing process?
 * Do we want to store the user pointer of the p4est and/or the p4est_wrap
 * structure?
 */

/** Create and open a file that is associated with the given domain structure.
 *
 * This is a collective function call that overwrites the file if it exists
 * already. This function writes a header with metadata on the underlying
 * p4est of \b domain to the file.
 *
 * The opened file can be used to write to the file using the functions
 * \ref fclaw2d_file_write_glob_options, \ref fclaw2d_file_write_domain,
 * \ref fclaw2d_file_write_patch_data.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]   domain    The underlying p4est is used for the metadata of the
 *                         the created file. Later function calls, e.g. the
 *                         writing of the domain must have \b domain as passed
 *                         to this function.
 * \param [in] filename    Path to parallel file that is to be created.
 * \param [in] user_string A user string that is written to the file header.
 *                         Only \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                         bytes without NUL-termination are
 *                         written to the file. If the user gives less
 *                         bytes the user_string in the file header is padded
 *                         by spaces.
 * \param [out] errcode    An errcode that can be interpreted by
 *                         \ref fclaw2d_file_error_string.
 * \return                 Newly allocated context to continue writing and
 *                         eventually closing the file. NULL in case of error.
 */
fclaw2d_file_context_t *fclaw2d_file_open_create (fclaw2d_domain_t * domain,
                                                  const char *filename,
                                                  const char *user_string,
                                                  int *errcode);

/** Open a file for reading and read its user string on rank zero.
 * The user string is broadcasted to all ranks after reading.
 * The file must exist and be at least of the size of the file header.
 *
 * If the file has wrong metadata the function reports the error using
 * /ref P4EST_LERRORF, collectively close the file and deallocate
 * the file context. In this case the function returns NULL on all ranks.
 * The wrong file format or a wrong file header causes \ref P4EST_FILE_ERR_FORMAT
 * as errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  mpicomm       MPI communicator that is used to read the file and
 *                            must be used for potentially later read domain
 *                            and glob.
 * \param [in]  filename      The path to the file that is opened.
 * \param [out] user_string   At least \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \return                    Newly allocated context to continue reading
 *                            and eventually closing the file. NULL in
 *                            case of error.
 */
fclaw2d_file_context_t *fclaw2d_file_open_read (sc_MPI_Comm mpicomm,
                                                const char *filename,
                                                char *user_string,
                                                int *errcode);

/** Write a domain to an opened parallel file.
 *
 * This function writes a domain without the patch data to an opened parallel
 * file. One can only write the domain that was used by openeing the file using
 * \ref fclaw2d_file_open_create.
 *
 * The mesh data is written in parallel using the partition of the mesh, i.e.
 * the domain and the underlying p4est, repectivley.
 *
 * See \ref fclaw2d_file_write_patch_data to write the patch data of \b domain.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc          Context previously created by \ref
 *                              fclaw2d_file_open_create.  It keeps track
 *                              of the data sets written one after another.
 * \param [in]      domain      The domain that is written to the file. This
 *                              function does not write the patch data of
 *                              \b domain. See \ref
 *                              fclaw2d_file_write_patch_data to write the patch
 *                              data. \b domain must coincide with the domain
 *                              that was passed for opening the file.
 * \param [in]      user_string A user string that is written to the file.
 *                              Only \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                              bytes without NUL-termination are
 *                              written to the file. If the user gives less
 *                              bytes the user_string in the file header is padded
 *                              by spaces.
 * \param [in]      errcode     An errcode that can be interpreted by
 *                              \ref fclaw2d_file_error_string.
 * \return                      Return a pointer to input context or NULL in case
 *                              of errors that does not abort the program.
 *                              In case of error the file is tried to close
 *                              and \b fc is freed.
 */
fclaw2d_file_context_t *fclaw2d_file_write_domain (fclaw2d_file_context_t *
                                                   fc,
                                                   fclaw2d_domain_t * domain,
                                                   const char *user_string,
                                                   int *errcode);

/** Read a domain from an opened file using the MPI communicator of \b fc.
 *
 * The read domain does not have patch data, which can be read by
 * \ref fclaw2d_file_read_patch_data.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  fc            Context previously created by \ref
 *                            fclaw2d_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [in]  filename      The path to the file that is opened.
 * \param [out] user_string   At least \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in] replace_fn     Qudrant replace callback that is passed to
 *                            \ref p4est_wrap_new_p4est. This callback can be
 *                            NULL.
 * \param [in] attributes     Attributes of the read domain. Currently the
 *                            attributes of the domain are not written to the
 *                            file by \ref fclaw2d_file_write_domain and
 *                            therefore must be passed to this function. Can be
 *                            NULL.
 * \param [in] wrap_user_pointer  A pointer to anonymous user data that is
 *                            passed to \ref p4est_wrap_new_p4est to create
 *                            the p4est_wrap structure that is used as element
 *                            of the read domain.
 * \param [out] domain        Newly allocated domain that is read from the file.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
fclaw2d_file_context_t *fclaw2d_file_read_domain (fclaw2d_file_context_t * fc,
                                                  const char *filename,
                                                  char *user_string,
                                                  p4est_replace_t replace_fn,
                                                  sc_keyvalue_t * attributes,
                                                  void *wrap_user_pointer,
                                                  fclaw2d_domain_t ** domain,
                                                  int *errcode);

/** Write patch data to an opened parallel file.
 *
 * This function writes the patch data of a given domain to the opened file.
 * \b domain must coincide with the domain that was passed in for opening the
 * file by \ref fclaw2d_file_open_create.
 *
 * The patch data is written in parallel according to the partition of the
 * domain and the underlying p4est, respectively.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc          Context previously created by \ref
 *                              fclaw2d_file_open_create.  It keeps track
 *                              of the data sets written one after another.
 * \param [in]      user_string A user string that is written to the file.
 *                              Only \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                              bytes without NUL-termination are
 *                              written to the file. If the user gives less
 *                              bytes the user_string in the file header is padded
 *                              by spaces.
 * \param [in]      domain      The domain that points to the patch data that
 *                              is written to the opened parallel file.
 * \param [in] patch_data_size  The patch data size in number of bytes.
 * \param [out]     errcode     An errcode that can be interpreted by
 *                              \ref fclaw2d_file_error_string.
 * \return                      Return a pointer to input context or NULL in case
 *                              of errors that does not abort the program.
 *                              In case of error the file is tried to close
 *                              and \b fc is freed.
 */
fclaw2d_file_write_t *fclaw2d_file_write_patch_data (fclaw2d_file_context_t *
                                                     fc, char *user_string,
                                                     fclaw2d_domain_t *
                                                     domain,
                                                     size_t patch_data_size,
                                                     int *errcode);

/** Read a patch data from an opened file using the MPI communicator of \b fc.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  fc            Context previously created by \ref
 *                            fclaw2d_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [out] user_string   At least \ref FCLAW2D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in, out] domain    The domain for that the patch data is read.
 *                            The patch data will be assigned to the domain
 *                            after a successful call of this function.
 * \param [in] patch_data_size The patch data size in number of bytes.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 * \param [out]     errcode   An errcode that can be interpreted by
 *                            \ref fclaw2d_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and \b fc is freed.
 * \note                      The patch data size can not be read from
 *                            the file since it is not written by \ref
 *                            fclaw2d_file_write_patch_data. In practice
 *                            one can use \ref
 *                            fclaw2d_file_write_global_opt and then
 *                            retrieve from this information the patch
 *                            data size and call this function. See also
 *                            the implementation of \ref
 *                            fclaw2d_file_read_global.
 */
fclaw2d_file_write_t *fclaw2d_file_read_patch_data (fclaw2d_file_context_t *
                                                    fc, char *user_string,
                                                    fclaw2d_domain_t * domain,
                                                    size_t patch_data_size,
                                                    int *errcode);

/** Close a file opened for parallel write/read and free the context.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in,out] fc       Context previously created by \ref
 *                          fclaw2d_file_open_create or \ref
 *                          fclaw2d_file_open_read.  Is freed.
 * \param [out] errcode     An errcode that can be interpreted by \ref
 *                          fclaw2d_file_error_string.
 * \return                  0 for a successful call and -1 in case of
 *                          an error. See also errcode argument.
 */
int fclaw2d_file_close (fclaw2d_file_context_t * fc, int *errcode);

/** Turn fclaw2d_file errcode into a string.
 *
 * \param [in] errcode      An errcode that is output by a
 *                          fclaw_file function.
 * \param [in,out] string   At least sc_MPI_MAX_ERROR_STRING bytes.
 * \param [out] resultlen   Length of string on return.
 * \return                  \ref 1 on success or
 *                          something else on invalid arguments.
 */
int
fclaw2d_file_error_string (int errcode, char *string, int *resultlen);

fclaw2d_domain_t *fclaw2d_domain_new_unitsquare (sc_MPI_Comm mpicomm,
                                                 int initial_level);

fclaw2d_domain_t *fclaw2d_domain_new_torus (sc_MPI_Comm mpicomm,
                                            int initial_level);

fclaw2d_domain_t *fclaw2d_domain_new_twosphere (sc_MPI_Comm mpicomm,
#if 0
     int fclaw2d_domain_write (const char *filename,
                               fclaw2d_domain_t * domain);

     fclaw2d_domain_t *fclaw2d_domain_read (sc_MPI_Comm mpicomm,
                                            const char *filename,
                                            p4est_replace_t replace_fn,
                                            sc_keyvalue_t * attributes,
                                            void *wrap_user_pointer);
#endif

                                                int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_cubedsphere (sc_MPI_Comm mpicomm,
                                                  int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_disk (sc_MPI_Comm mpicomm,
                                           int periodic_in_x,
                                           int periodic_in_y,
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
 * \return                      A fully initialized domain structure.
 */
fclaw2d_domain_t *fclaw2d_domain_new_brick (sc_MPI_Comm mpicomm,
                                            int blocks_in_x, int blocks_in_y,
                                            int periodic_in_x,
                                            int periodic_in_y,
                                            int initial_level);

/** Create a domain from a given forest connectivity.
 * \param [in] mpicomm          We expect sc_MPI_Init to be called earlier.
 * \param [in] initial_level    A non-negative integer <= P4EST_QMAXLEVEL.
 * \param [in] conn             We DO take ownership of the connectivity.
 * \return                      A fully initialized domain structure.
 */
fclaw2d_domain_t *fclaw2d_domain_new_conn (sc_MPI_Comm mpicomm,
                                           int initial_level,
                                           p4est_connectivity_t * conn);

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

/** Search triples of (block number, x, y coordinates) in the mesh.
 * The x, y coordinates must be in [0, 1]^2.
 * The input data must be equal on every process: This is a collective call.
 *
 * A point is found at most once even if it is on a patch boundary.
 * We return the smallest patch number on the smallest processor touching it.
 * However, if a point is on a block boundary, it must be decided before
 * calling this function which tree shall be queried for it.
 *
 * \param [in] domain           Must be valid domain structure.  Will not be changed.
 * \param [in] block_offsets    Monotonous array of (num_blocks + 1) int variables.
 *                              The points to search in block t in [0, num_blocks)
 *                              have indices [block_offsets[t], block_offsets[t + 1])
 *                              in the \b coordinates and results arrays.
 * \param [in] coordinates      An array of elem_size == 2 * sizeof (double) with
 *                              entries (x, y) in [0, 1]^2.  Of these entries,
 *                              there are \b block_offsets[num_blocks] many.
 *                              We do not enforce the x and y ranges
 *                              and simply do not find any point outside its block.
 * \param [in,out] results      On input, an array of type int and an element
 *                              count of \b block_offsets[num_blocks].
 *                              The data in \b results is ignored on input.
 *                              On output, an entry will be -1 if the point has
 *                              not been found on this process, or the patch
 *                              number within its block otherwise.
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

/** Callback function to compute the interpolation data for a point and a patch.
 *
 * This function can be passed to \ref fclaw2d_overlap_exchange to eventually
 * compute the interpolation data over the whole producer domain for an
 * array of points.
 * It will be called both in a partition search and a local search of the
 * producer domain. Use \ref fclaw2d_domain_is_meta, to determine which is the
 * case.
 *
 * \param [in] domain           The domain we interpolate on.
 *                              On the producer side, this is a valid forestclaw
 *                              domain.
 *                              On the consumer side, this is a temporary
 *                              artifical domain. Only the mpi-information
 *                              (mpicomm, mpisize and mpirank) as well as the
 *                              backend data (pp, pp_owned and attributes) are
 *                              set. The backend data is not owned and shall
 *                              not be changed by the callback. The mpirank is
 *                              set to a valid rank only when we are at a leaf
 *                              (a patch that belongs to exactly one process)
 *                              of the partition search, else it will be -1.
 * \param [in] patch            The patch under consideration.
 *                              When on a leaf on the producer side, this is a
 *                              valid patch from the producer domain.
 *                              Otherwise, this is a temporary artificial patch
 *                              containing all standard patch information except
 *                              for the pointer to the next patch and user-data.
 *                              Only the FCLAW2D_PATCH_CHILDID and the
 *                              FCLAW2D_PATCH_ON_BLOCK_FACE_* flags are set.
 *                              Artificial patches are generally ancestors of
 *                              valid forestclaw patches that are leaves.
 * \param [in] blockno          The block id of the patch under consideration.
 * \param [in] patchno          If patchno is -1, we are on an artifical patch.
 *                              Otherwise, this is a valid patchno from the
 *                              producer domain.
 * \param [in,out] point        Representation of a point; user-defined.
 *                              Points to an array element of the query points
 *                              passed to \ref fclaw2d_overlap_exchange.
 *                              If patchno is non-negative, the points
 *                              interpolation data should be updated by the
 *                              local patch's contribution.
 * \param [in, out] user        Arbitrary data passed in earlier.
 * \return                      True, if there is a possible contribution of the
 *                              patch or one of its ancestors to the point
 *                              interpolation data.
 *                              Return false if there is definitely no
 *                              contribution.
 *                              If we are on a leaf on the producer side
 *                              (patchno is non-negative) or the consumer side
 *                              (domain_is_meta and mpirank is non-negative)
 *                              this callback should do an exact test for
 *                              contribution.
 *                              Else, the return value may be a false positive,
 *                              we'll be fine.
 */
typedef int (*fclaw2d_interpolate_point_t) (fclaw2d_domain_t * domain,
                                            fclaw2d_patch_t * patch,
                                            int blockno, int patchno,
                                            void *point, void *user);

/** Exchange interpolation data of query points between two domains.
 *
 * We compute the user-defined interpolation data of an array of user-defined
 * query points, which originate from the so-called consumer side.
 * The interpolation data will be computed on the domain of the so-called
 * producer-side based on a \ref fclaw2d_interpolate_point_t callback function.
 * Afterwards, the results will be collected and combined on the consumer side.
 *
 * \param [in] domain           The producer domain to interpolate on.
 * \param [in,out] query_points Array containing points of user-defined type.
 *                              Each entry contains one item of arbitrary data.
 *                              We do not dereference, just pass pointers around.
 *                              The points will be sent via MPI, so they may not
 *                              contain pointers to further data.
 *                              The array is defined processor-local and may
 *                              contain different points on different processes.
 *                              The query points are supposed to be computed
 *                              (and transformed to the producer space by an
 *                              inverse mapping) locally on the consumer side.
 *                              On output, the points will contain collected
 *                              interpolation data according to \b interpolate.
 * \param [in] interpolate      Callback function that returns true if a point
 *                              intersects a patch and -- when called for a leaf
 *                              on the producer side -- shall write the
 *                              interpolation data for the current
 *                              point-patch-combination into the user-defined
 *                              point structure.
 * \param [in,out] user         Arbitrary data to be passed to the callback.
 */
void fclaw2d_overlap_exchange (fclaw2d_domain_t * domain,
                               sc_array_t * query_points,
                               fclaw2d_interpolate_point_t interpolate,
                               void *user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CONVENIENCE_H */
