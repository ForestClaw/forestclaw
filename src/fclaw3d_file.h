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

/** \file fclaw3d_file.h
 * Routines for parallel IO for \ref fclaw3d_domain_t and generic data.
 *
 * Principle: All files are associated with a single domain, i.e. one can
 * write only one domain to one file.
 */

#ifndef FCLAW3D_FILE_H
#define FCLAW3D_FILE_H

#include <forestclaw3d.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#define FCLAW3D_FILE_USER_STRING_BYTES 48 /**< number of user string bytes */
#define FCLAW3D_FILE_MAX_BLOCK_SIZE ((1000 * 1000 * 1000 * 1000 * 10) - 1) /**< maximal data size of a block */
#define FCLAW3D_FILE_MAX_FIELD_ENTRY_SIZE ((1000 * 1000 * 1000 * 1000 * 10) - 1) /**< maximal data size per field entry*/

/** Error values for fclaw2d_file functions.
 * In the future we may add further error codes. Therefore, the error codes
 * should only be used by the enum but not by the explicit numeric values since
 * these explicit numeric values may change.
 */
typedef enum fclaw3d_file_error
{
    FCLAW3D_FILE_ERR_SUCCESS = 0, /**< file function completed with success */
    FCLAW3D_FILE_ERR_FILE = sc_MPI_ERR_LASTCODE, /**< invalid file handle */
    FCLAW3D_FILE_ERR_NOT_SAME, /**< collective arg not identical */
    FCLAW3D_FILE_ERR_AMODE, /**< access mode error */
    FCLAW3D_FILE_ERR_NO_SUCH_FILE, /**< file does not exist */
    FCLAW3D_FILE_ERR_FILE_EXIST, /**< file exists already */
    FCLAW3D_FILE_ERR_BAD_FILE, /**< invalid file name */
    FCLAW3D_FILE_ERR_ACCESS, /**< permission denied */
    FCLAW3D_FILE_ERR_NO_SPACE, /**< not enough space */
    FCLAW3D_FILE_ERR_QUOTA, /**< quota exceeded */
    FCLAW3D_FILE_ERR_READ_ONLY, /**< read only file (system) */
    FCLAW3D_FILE_ERR_IN_USE, /**< file currently open by other process */
    FCLAW3D_FILE_ERR_IO, /**< other I/O error */
    FCLAW3D_FILE_ERR_FORMAT, /**< read file has a wrong format */
    FCLAW3D_FILE_ERR_SECTION_TYPE, /**< a valid non-matching section type */
    FCLAW3D_FILE_ERR_CONN, /**< invalid serialized connectivity data */
    FCLAW3D_FILE_ERR_P4EST, /**< invalid p4est data */
    FCLAW3D_FILE_ERR_IN_DATA, /**< input data of file function is invalid */
    FCLAW3D_FILE_ERR_COUNT, /**< read or write count error that was not
                                 classified as a format error */
    FCLAW3D_FILE_ERR_DIM, /**< file has wrong dimension */
    FCLAW3D_FILE_ERR_UNKNOWN, /**< unknown error */
    FCLAW3D_FILE_ERR_NOT_IMPLEMENTED, /**< functionality is not implemented */
    FCLAW3D_FILE_ERR_LASTCODE /**< to define own error codes for
                                  a higher level application
                                  that is using fclaw3d_file
                                  functions */
}
fclaw3d_file_error_t;

/** Opaque context used for writing a fclaw3d data file. */
typedef struct fclaw3d_file_context fclaw3d_file_context_t;

/** Create and open a file that is associated with the given domain structure.
 *
 * This is a collective function call that overwrites the file if it exists
 * already. This function writes a header with metadata of the underlying
 * p8est of the passed \b domain to the file. In addition, the passed \b domain
 * is written to the file.
 *
 * The opened file can be used to write to the file using the functions
 * \ref fclaw3d_file_write_block, \ref fclaw3d_file_write_array and
 * potentially other functions operating on \ref fclaw3d_file_context_t that
 * may be implemented on the top of the current fclaw3d_file routines.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in] filename    Path to parallel file base name that is to be created.
 *                         This means that the user shall not pass the file
 *                         extension since the file extension is specified by
 *                         fclaw3d_file to '.f3d'.
 * \param [in] user_string A NUL-terminated user string that is written to the
 *                         file header having FCLAW3D_FILE_USER_STRING_BYTES
 *                         bytes including the NUL-termination. Shorter
 *                         user strings are padded by spaces in the file header.
 *                         Too long user strings result in an error with the
 *                         error code \ref FCLAW3D_FILE_ERR_IN_DATA.
 * \param [in]   domain    The underlying p8est is used for the metadata of the
 *                         the created file and the \b domain is written to the
 *                         file.
 * \param [out] errcode    An errcode that can be interpreted by
 *                         \ref fclaw3d_file_error_string.
 * \return                 Newly allocated context to continue writing and
 *                         eventually closing the file. NULL in case of error.
 */
fclaw3d_file_context_t *fclaw3d_file_open_write (const char *filename,
                                                 const char *user_string,
                                                 fclaw3d_domain_t * domain,
                                                 int *errcode);

#if 0
/** Write the partition of the domain associated with \b fc to disk.
 *
 * This function derives the partition filename from the base filename 'fn' of
 * the given file context \b fc such that the function creates a file named
 * 'fn_partition.f3d' that will contain the partition of the domain that is
 * associated with the given file context \b fc.
 *
 * This means in particular that the given \b fc stays unchanged and that the
 * internally creates a new file context that will be freed and deallocated at
 * the end of the function in any case. This is important to notice since the
 * \b errcode output parameter of this function refers to this internally created
 * file context and not to the \b fc.
 *
 * \note In contrast to the not-partition-related writing functions this
 * function always returns the passed \b fc, which stays unchanged, even if an
 * error occurred since this error then relates to the internal file context
 * for writing the separate partition file.
 * \b errcode always relates to the internal file context and not to the passed
 * \b fc.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]     filename     Path to partition file base name that is to be
 *                              created. This means that the user shall not pass
 *                              the file extension since the partition file
 *                              extension is specified by fclaw3d_file as '.fp3d'.
 * \param [in]     user_string  A NUL-terminated user string that is written to
 *                              the partition file header having \ref
 *                              FCLAW3D_FILE_USER_STRING_BYTES bytes including
 *                              the NUL-termination. Shorter user strings are
 *                              padded by spaces in the file header.
 *                              Too long user strings result in an error with the
 *                              error code \ref FCLAW3D_FILE_ERR_IN_DATA.
 * \param [out]    errcode      An errcode that can be interpreted by
 *                              \ref fclaw3d_file_error_string. \b errcode
 *                              does not refer to the passed but it refers to
 *                              the internal file context used for for the
 *                              partition file.
 * \return                      0 for a successful write of the partition file.
 *                              -1 in case of an error. See also \b errcode
 *                              to examine the error of the partition file
 *                              I/O process.
 */
int fclaw3d_file_write_partition (const char *filename,
                                  const char *user_string,
                                  fclaw3d_domain_t * domain, int *errcode);
#endif

/** Write a serial data block to an opened parallel file.
 *
 * This is a collective function.
 * This function writes a serial data block to the opened file.
 *
 * The number of block bytes must be less or equal \ref
 * FCLAW3D_FILE_MAX_BLOCK_SIZE.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc          Context previously created by \ref
 *                              fclaw3d_file_open_write.  It keeps track
 *                              of the data sets written one after another.
 * \param [in]      user_string A NUL-terminated user string that is written to
 *                              the section header having
 *                              FCLAW3D_FILE_USER_STRING_BYTES bytes including
 *                              the NUL-termination. Shorter user strings are
 *                              padded by spaces in the section header. Too long
 *                              user strings result in an error with the error
 *                              code \ref FCLAW3D_FILE_ERR_IN_DATA.
 * \param [in]      block_size  The size of the block in bytes. May be equal to
 *                              0. In this case the section header and the padding
 *                              is still written. This function returns the passed
 *                              \b fc parameter and sets errcode to \ref
 *                              FCLAW3D_FILE_ERR_SUCCESS if it is called for
 *                              \b block_size == 0.
 * \param [in]      block_data  A sc_array with one element and element size
 *                              equal to \b block_size.
 *                              The array points to the block data. The user is
 *                              responsible for the validality of the block
 *                              data. block_data can be NULL if \b block_size == 0.
 * \param [out]     errcode     An errcode that can be interpreted by
 *                              \ref fclaw3d_file_error_string.
 * \return                      Return a pointer to input context or NULL in case
 *                              of errors that does not abort the program.
 *                              In case of error the file is tried to close
 *                              and \b fc is freed.
 */
fclaw3d_file_context_t *fclaw3d_file_write_block (fclaw3d_file_context_t *
                                                  fc, const char *user_string,
                                                  size_t block_size,
                                                  sc_array_t * block_data,
                                                  int *errcode);

/** Write per-patch data to a parallel output file.
 *
 * This is a collective function.
 * This function writes the per-patch data with respect the domain that was passed
 * to the \ref fclaw3d_file_open_write for opening the file.
 *
 * The per-patch data is written in parallel according to the partition of the
 * domain and the underlying p8est, respectively.
 * The data size per patch must be fixed.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc          Context previously created by \ref
 *                              fclaw3d_file_open_write.  It keeps track
 *                              of the data sets written one after another.
 * \param [in]      user_string A NUL-terminated user string that is written to
 *                              the section header having
 *                              FCLAW3D_FILE_USER_STRING_BYTES bytes including
 *                              the NUL-termination. Shorter user strings are
 *                              padded by spaces in the section header. Too long
 *                              user strings result in an error with the error
 *                              code \ref FCLAW3D_FILE_ERR_IN_DATA.
 * \param [in]      patch_size  The number of bytes per patch. This number
 *                              must coincide with \b patch_data->elem_size.
 * \param [in]      patch_data  An array of the length number of local patches
 *                              with the element size equal to sizeof (sc_array_t).
 *                              This means each element of \b patch_data must be
 *                              a sc_array. These sc_arrays must have an element
 *                              size of \b patch_size and store the actual
 *                              patch data. This memory layout enables the user
 *                              to write by indirect addressing, i.e. to write
 *                              non-contigous patch data. The patch data is
 *                              expected to be stored according to the Morton
 *                              order of the patches. For \b patch_size == 0
 *                              the function writes an empty array. The section
 *                              header and the padding is still written.
 *                              In this case the function passes successfully and
 *                              \b errcode is set to \ref FCLAW3D_FILE_ERR_SUCCESS.
 * \param [out]     errcode     An errcode that can be interpreted by
 *                              \ref fclaw3d_file_error_string.
 * \return                      Return a pointer to input context or NULL in case
 *                              of errors that does not abort the program.
 *                              In case of error the file is tried to close
 *                              and \b fc is freed.
 */
fclaw3d_file_context_t *fclaw3d_file_write_array (fclaw3d_file_context_t *
                                                  fc, const char *user_string,
                                                  size_t patch_size,
                                                  sc_array_t * patch_data,
                                                  int *errcode);

/** Open a file for reading and read the stored domain.
 * The user string is broadcasted to all ranks after reading.
 * The file must exist and be at least of the size of the file header.
 *
 * This is a collective function.
 * If the file has wrong metadata the function reports the error using
 * \ref fclaw_errorf, collectively closes the file and deallocates
 * the file context. In this case the function returns NULL on all ranks.
 * The wrong file format or a wrong file header causes \ref FCLAW3D_FILE_ERR_FORMAT
 * as errcode.
 *
 * After calling this function the user can continue reading the opened file
 * by calling \ref fclaw3d_file_read_block, \ref fclaw3d_file_read_array and
 * potentially other functions operating on \ref fclaw3d_file_context_t that
 * may be implemented on the top of the current fclaw3d_file routines.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  filename      The path to the base name file that is opened.
 *                            This means that the user shall not pass the file
 *                            extension since the extension ('f3d') is added by
 *                            fclaw3d_file.
 * \param [out] user_string   At least \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in]  mpicomm       MPI communicator that is used to read the file and
 *                            is used for potential other reading operations of
 *                            MPI communicator dependent objects.
 * \param [in]  par_filename  The path to the base name file that is used to
 *                            load the parititon. If \b par_filename is NULL,
 *                            a uniform parititon with respect to the patch
 *                            count is computed and used. If the partition is
 *                            read, it is used for the parallel I/O operations
 *                            and stored in the returned \b domain. If the MPI
 *                            size and the partition size do not coincide a
 *                            partition for the current MPI size is computed.
 *                            The function call results in an error if there
 *                            is no partition to read.
 *                            Currently, this parameter does not have any effect.
 * \param [out] domain        Newly allocated domain that is read from the file.
 * \param [out] errcode       An errcode that can be interpreted by
 *                            \ref fclaw3d_file_error_string.
 * \return                    Newly allocated context to continue reading
 *                            and eventually closing the file. NULL in
 *                            case of error.
 */
fclaw3d_file_context_t *fclaw3d_file_open_read (const char *filename,
                                                char *user_string,
                                                sc_MPI_Comm mpicomm,
                                                const char *par_filename,
                                                fclaw3d_domain_t ** domain,
                                                int *errcode);

/** Read a serial data block from an opened file.
 *
 * This is a collective function.
 *
 * This function requires an opened file context.
 *
 * The passed \b block_size is compared to the block size stored in the file.
 * If the values do not equal each other, the function reports details, closes
 * and deallocates the file context. The return value in this case is NULL.
 * If the data block section header information is not matching the passed
 * parameters the function sets \ref FCLAW3D_FILE_ERR_FORMAT for errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  fc            Context previously created by \ref
 *                            fclaw3d_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [out] user_string   At least \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in]  block_size    The size of the header that is read.
 * \param [in, out] block_data \b block_size allocated bytes in an sc_array
 *                            with one element and \b block_size as element
 *                            size. This data will be filled with the block
 *                            data from the file. If this is NULL it means that
 *                            the current header block is skipped and the
 *                            internal file pointer of the file context is
 *                            set to the next data section. If the current data
 *                            section is not a data block, the file is closed
 *                            and the file context is deallocated. Furthermore,
 *                            in this case the function returns NULL and sets
 *                            errcode to \ref FCLAW3D_FILE_ERR_FORMAT. In case
 *                            of skipping the data block section \b block_size
 *                            needs also to coincide with the data block size
 *                            given in the file.
 * \param [out]     errcode   An errcode that can be interpreted by
 *                            \ref fclaw3d_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and \b fc is freed.
 */
fclaw3d_file_context_t *fclaw3d_file_read_block (fclaw3d_file_context_t *
                                                 fc, char *user_string,
                                                 size_t block_size,
                                                 sc_array_t * block_data,
                                                 int *errcode);

/** Read per-patch data from a parallel output file.
 *
 * This is a collective function.
 * The function closes and deallocates the file context and returns NULL
 * if the bytes the user wants to read exceed the given file and/or
 * the element size of the array given by \b patch_data->elem_size does not
 * coincide with the element size according to the array metadata given in
 * the file.
 *
 * The data is read in parallel using the partition of the domain (and the
 * underlying p8est) that is associated with the passed file context in the
 * sense that the file context was created using the respective domain.
 * The data size per patch must be fixed.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  fc            Context previously created by \ref
 *                            fclaw3d_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [out] user_string   At least \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in] patch_size     The number of bytes per patch. This number
 *                            must coincide with \b patch_data->elem_size.
 * \param [in,out] patch_data An array of the length number of local patches
 *                            with the element size equal to sizeof (sc_array_t).
 *                            This means each element of \b patch_data must be
 *                            a sc_array. These sc_arrays must have an element
 *                            size of \b patch_size. This memory layout enables
 *                            the user to read by indirect addressing, i.e. to
 *                            read into non-contigous patch data. The patch data
 *                            is read according to the Morton order of the patches.
 *                            \b patch_size must coincide with the section data
 *                            size in the file. \b patch_data == NULL means that
 *                            the data is skipped and the internal file pointer
 *                            is incremented. In the case of skipping
 *                            \b patch_size is still checked using the
 *                            corresponding value read from the file. The data
 *                            is read using the partition of patches given by the
 *                            domain that is associated to \b fc.
 * \param [out]     errcode   An errcode that can be interpreted by
 *                            \ref fclaw3d_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and \b fc is freed.
 */
fclaw3d_file_context_t *fclaw3d_file_read_array (fclaw3d_file_context_t *
                                                 fc, char *user_string,
                                                 size_t patch_size,
                                                 sc_array_t * patch_data,
                                                 int *errcode);

/** Close a file opened for parallel write/read and free the context.
 *
 * This is a collective function.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in,out] fc       Context previously created by \ref
 *                          fclaw3d_file_open_write or \ref
 *                          fclaw3d_file_open_read.  Is freed.
 * \param [out] errcode     An errcode that can be interpreted by \ref
 *                          fclaw3d_file_error_string.
 * \return                  0 for a successful call and -1 in case of
 *                          an error. See also errcode argument.
 */
int fclaw3d_file_close (fclaw3d_file_context_t * fc, int *errcode);

/** Turn fclaw3d_file errcode into a string.
 *
 * \param [in] errcode      An errcode that is output by a
 *                          fclaw3d_file function.
 * \param [in,out] string   At least sc_MPI_MAX_ERROR_STRING bytes.
 * \param [out] resultlen   Length of string on return.
 * \return                  0 on success or
 *                          something else on invalid arguments.
 */
int fclaw3d_file_error_string (int errcode, char *string, int *resultlen);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW3D_FILE_H */
