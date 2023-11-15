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

#ifndef P4_TO_P8
#include <fclaw2d_file.h>
#include <fclaw2d_convenience.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_wrap.h>
#else
#include <fclaw3d_file.h>
#include <fclaw3d_convenience.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_wrap.h>
#endif

#define FCLAW2D_FILE_NAME_BYTES BUFSIZ /**< maximal number of filename bytes */
#ifndef P4_TO_P8
#define FCLAW2D_FILE_EXT "f2d" /**< file extension of fclaw2d data files */
#else
#define FCLAW2D_FILE_EXT    FCLAW3D_FILE_EXT
#define FCLAW3D_FILE_EXT "f3d" /**< file extension of fclaw3d data files */
#endif

/* For legacy and compatibility reasons we provide here the first version
 * of the file format that was originally implemented in p4est in p4est_io.c.
 * The origin of the these function is also the reason why we still use
 * "p4data0" as magic string for the data format.
 */

/* Begin of legacy file functions */

/** p4est_v1 data file format
 * All p4est data files have 64 bytes file header section at the beginning of the file.
 * The file header section is written to the file as string without NUL-termination
 * (called string*) and is therefore readable in a text editor.
 *
 * File Header (96 bytes):
 * 7 bytes magic number (p{4,8}data0) and 1 byte new line char.
 * 23 bytes p4est file version string* and 1 byte new line char.
 * 47 bytes user string*  and 1 byte new line char.
 * 16 bytes number of global quadrants.
 *
 * The file header section is padded by 16 bytes consisting of 1 byte
 * new line char succeeded by 14 bytes of spaces and 1 trailing byte
 * new line char.
 *
 * The actual data is stored in arrays corresponding to a mesh of a p4est
 * or in block sections that have a fixed user-defined size. The block
 * sections are written and read on rank 0.
 * One data field stores a fixed number of bytes of user-
 * defined data per quadrant of a certain p4est. Therefore, one user-defined
 * data field is of the size p4est->global_num_quadrants * data_size, where
 * data_size is set by the user. The file format is partition independent.
 * The data fields are padded such that the number of bytes for
 * an array is divisible by 16. The padding also enforced for data blocks
 * that have a size that is divisble by 16.
 * The p4est data file consists of a variable number (including 0) of
 * these two types of data sections.
 * Every data section includes 64 bytes of section header written at the beginning
 * by p4est. These 64 bytes are again written to the file as string* and can
 * be read using a text editor.
 *
 * Data section Header (64 bytes):
 * One byte data section type specific character (B for a block section and F for
 * a data field), 1 byte space and 13 bytes size in number of bytes for a
 * block section and data size per element in byte for a field section
 * and one trailing byte new line char.
 * 47 bytes user-defined string* and 1 byte new line char.
 *
 * The structure of p4est and p8est data files differs only by the magic number.
 *
 * The p4est metadata of a p4est data file can be accessed by \ref
 * fclaw2d_file_info_v1().
 */

/* non-public macros */
#define FCLAW2D_FILE_MAGIC_NUMBER_V1 "p4data0" /**< magic string for fclaw2d data files */
#define FCLAW2D_FILE_METADATA_BYTES_V1 96 /**< number of file metadata bytes */
#define FCLAW2D_FILE_MAGIC_BYTES_V1 8 /**< number of bytes of the magic number */
#define FCLAW2D_FILE_VERSION_STR_BYTES_V1 24 /**< number of bytes of the version string*/
#define FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 14 /**< number of array metadata bytes */
/* subtract 2 for '\n' at the beginning and end of the array metadata */
#define FCLAW2D_FILE_ARRAY_METADATA_CHARS_V1 (FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 - 2) /**< number of array metadata chars */
#define FCLAW2D_FILE_BYTE_DIV_V1 16 /**< All data blocks are padded to be divisible by this. */
#define FCLAW2D_FILE_MAX_NUM_PAD_BYTES_V1 (FCLAW2D_FILE_BYTE_DIV_V1 + 1) /**< We enforce to pad in any
                                                               case and the padding string
                                                               needs to contain two
                                                               newline characters and
                                                               therefore this is the
                                                               maximal number of pad
                                                               bytes. */
#define FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 (2 + FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + FCLAW2D_FILE_USER_STRING_BYTES_V1)
                                     /**< number of bytes of one field header */
#define FCLAW2D_FILE_USER_STRING_BYTES_V1 48 /**< number of user string bytes */
#define FCLAW2D_FILE_MAX_BLOCK_SIZE_V1 9999999999999 /**< maximal data size of a block */
#define FCLAW2D_FILE_MAX_FIELD_ENTRY_SIZE_V1 9999999999999 /**< maximal data size per field entry*/
#define FCLAW2D_FILE_MAX_GLOBAL_QUAD_V1 9999999999999999 /**< maximal number of global quadrants */

#define FCLAW2D_FILE_HEADER_STRING_V1 "p4est data file v1"

#ifndef P4_TO_P8
#define FCLAW2D_FILE_STRING_V1 "fclaw2d_file_v1" /**< fclaw2d_file_v1 string */
#else
#define FCLAW2D_FILE_STRING_V1     FCLAW3D_FILE_STRING_V1
#define FCLAW3D_FILE_STRING_V1 "fclaw3d_file_v1" /**< fclaw3d_file_v1 string */
#endif

/** Error values for fclaw2d_file_v1 functions. These values are dimension
 * inependent and therefore we do not translate them since they are used only
 * internal.
 */
typedef enum fclaw2d_file_error_v1
{
    FCLAW2D_FILE_ERR_SUCCESS_V1 = sc_MPI_ERR_LASTCODE, /**< file function completed with success */
    FCLAW2D_FILE_ERR_FILE_V1, /**< invalid file handle */
    FCLAW2D_FILE_ERR_NOT_SAME_V1, /**< collective arg not identical */
    FCLAW2D_FILE_ERR_AMODE_V1, /**< access mode error */
    FCLAW2D_FILE_ERR_NO_SUCH_FILE_V1, /**< file does not exist */
    FCLAW2D_FILE_ERR_FILE_EXIST_V1, /**< file exists already */
    FCLAW2D_FILE_ERR_BAD_FILE_V1, /**< invalid file name */
    FCLAW2D_FILE_ERR_ACCESS_V1, /**< permission denied */
    FCLAW2D_FILE_ERR_NO_SPACE_V1, /**< not enough space */
    FCLAW2D_FILE_ERR_QUOTA_V1, /**< quota exceeded */
    FCLAW2D_FILE_ERR_READ_ONLY_V1, /**< read only file (system) */
    FCLAW2D_FILE_ERR_IN_USE_V1, /**< file currently open by other process */
    FCLAW2D_FILE_ERR_IO_V1, /**< other I/O error */
    FCLAW2D_FILE_ERR_FORMAT_V1, /**< read file has a wrong format */
    FCLAW2D_FILE_ERR_SECTION_TYPE_V1, /**< a valid non-matching section type */
    FCLAW2D_FILE_ERR_CONN_V1, /**< invalid serialized connectivity data */
    FCLAW2D_FILE_ERR_P4EST_V1, /**< invalid p4est data */
    FCLAW2D_FILE_ERR_IN_DATA_V1, /**< input data of file function is invalid */
    FCLAW2D_FILE_ERR_COUNT_V1, /**< read or write count error that was not
                                 classified as a format error */
    FCLAW2D_FILE_ERR_DIM_V1, /**< file has wrong dimension */
    FCLAW2D_FILE_ERR_UNKNOWN_V1, /**< unknown error */
    FCLAW2D_FILE_ERR_LASTCODE_V1 /**< to define own error codes for
                                  a higher level application
                                  that is using fclaw2d_file_v1
                                  functions */
}
fclaw2d_file_error_v1_t;

/* Avoid redefinition in translation file */
#ifdef P4_TO_P8
#define fclaw2d_file_context_p4est_v1            fclaw3d_file_context_p8est_v1
#endif

/** Translate dimension-dependent static functions. These IO functions are
 * dimension dependent in the sense that their API depends directly or
 * indirectly on p4est/p8est
 */
#ifdef P4_TO_P8
#define fclaw2d_file_open_create_v1   fclaw3d_file_open_create_v1
#define fclaw2d_file_open_read_ext_v1 fclaw3d_file_open_read_ext_v1
#define fclaw2d_file_open_read_v1     fclaw3d_file_open_read_v1
#define fclaw2d_file_write_block_v1   fclaw3d_file_write_block_v1
#define fclaw2d_file_read_block_metadata_v1 fclaw3d_file_read_block_metadata_v1
#define fclaw2d_file_read_block_v1    fclaw3d_file_read_block_v1
#define fclaw2d_file_write_field_v1   fclaw3d_file_write_field_v1
#define fclaw2d_file_read_field_ext_v1 fclaw3d_file_read_field_ext_v1
#define fclaw2d_file_read_field_v1    fclaw2d_file_read_field_v1
#define fclaw2d_file_info_v1          fclaw3d_file_info_v1
#define fclaw2d_file_write_p4est_v1   fclaw3d_file_write_p4est_v1
#define fclaw2d_file_read_p4est_v1    fclaw3d_file_read_p4est_v1
#define fclaw2d_file_write_connectivity_v1 fclaw3d_file_write_connectivity_v1
#define fclaw2d_file_read_connectivity_v1 fclaw3d_file_read_connectivity_v1
#define fclaw2d_file_close_v1         fclaw3d_file_close_v1
#endif

/** The file context for for fclaw2d data files associated to a p4est. */
typedef struct fclaw2d_file_context_p4est_v1
{
    sc_MPI_Comm mpicomm;                  /**< corresponding MPI communicator */
    p4est_locidx_t local_num_quadrants;    /**< number of local quadrants */
    p4est_gloidx_t global_num_quadrants;    /**< number of global quadrants */
    p4est_gloidx_t *global_first_quadrant;   /**< represents the partition */
    int gfq_owned;                        /**< Boolean to indicate if global_first_quadrant
                                               is owned. */
    size_t num_calls;                   /**< redundant but for convenience;
                                            counts the number of calls of
                                            write and read, respectively */
    sc_MPI_File file;                   /**< file object */
    sc_MPI_Offset accessed_bytes;       /**< count only array data bytes and
                                           array metadata bytes */
}
fclaw2d_file_context_p4est_v1_t;

#define FCLAW2D_FILE_COMPRESSED_QUAD_SIZE_V1 ((P4EST_DIM + 1) *\
                                        sizeof (p4est_qcoord_t))
                                        /**< size of a compressed quadrant */

/* error checking macros for fclaw2d_file functions */

#define FCLAW2D_FILE_IS_SUCCESS_V1(errcode) ((errcode == sc_MPI_SUCCESS)\
                                         || (errcode == FCLAW2D_FILE_ERR_SUCCESS_V1))

/** Declare error code function that is required for the error handling macros. */
static int fclaw2d_file_error_code_v1 (int errcode, int *fclaw2d_errcode);
/** Declare the error cleanup function that is required for the error handling
 * macros.
 */
static int fclaw2d_file_error_cleanup_v1 (sc_MPI_File * file);
/** Declare error string function that is required for the error handling macros. */
static int fclaw2d_file_error_string_v1 (int errclass, char *string,
                                         int *resultlen);

/** Examine the fclaw2d file return value and print an error if there is one.
 * The message passed is appended to fclaw2d file, file and line information.
 */
#define FCLAW2D_FILE_CHECK_VERBOSE_V1(errcode,user_msg) do {          \
  char fclaw2d_msg[sc_MPI_MAX_ERROR_STRING];                       \
  int fclaw2d_msglen;                                              \
  if (!FCLAW2D_FILE_IS_SUCCESS_V1 (errcode)) {                        \
    fclaw2d_file_error_code_v1 (errcode, &errcode);                   \
    fclaw2d_file_error_string_v1 (errcode, fclaw2d_msg, &fclaw2d_msglen); \
    SC_GLOBAL_LERRORF ("%s at %s:%d: %s\n",                      \
                      (user_msg), __FILE__, __LINE__, fclaw2d_msg);\
  }} while (0)

/** This macro performs a clean up in the case of a MPI I/O open error.
 * We make use of the fact that sc_mpi_open is always called collectively.
 */
#define FCLAW2D_FILE_CHECK_OPEN_V1(errcode, fc, user_msg, cperrcode) do {\
                                            FCLAW2D_FILE_CHECK_VERBOSE_V1 (errcode, user_msg);\
                                            *cperrcode = errcode;                       \
                                            if (!FCLAW2D_FILE_IS_SUCCESS_V1 (errcode)) {     \
                                            fclaw2d_file_error_cleanup_v1 (&fc->file);       \
                                            FCLAW_FREE (fc);                            \
                                            fclaw2d_file_error_code_v1 (errcode, cperrcode); \
                                            return NULL;}} while (0)

/** The same as \ref FCLAW2D_FILE_CHECK_OPEN_V1 but returns -1 instead of NULL */
#define FCLAW2D_FILE_CHECK_INT_V1(errcode, user_msg, cperrcode) do {\
                                            FCLAW2D_FILE_CHECK_VERBOSE_V1 (errcode, user_msg);   \
                                            *cperrcode = errcode;                       \
                                            if (!FCLAW2D_FILE_IS_SUCCESS_V1 (errcode)) {     \
                                            fclaw2d_file_error_code_v1 (errcode, cperrcode); \
                                            return -1;}} while (0)

/** This macro prints the MPI error for sc_mpi_{read,write}_all and return NULL.
 * This means that this macro is appropriate to call it after a collective
 * read or write.
 */
#define FCLAW2D_FILE_CHECK_NULL_V1(errcode, fc, user_msg, cperrcode) do {\
                                            FCLAW2D_FILE_CHECK_VERBOSE_V1 (errcode, user_msg);\
                                            *cperrcode = errcode;                       \
                                            if (!FCLAW2D_FILE_IS_SUCCESS_V1 (errcode)) {     \
                                            fclaw2d_file_error_cleanup_v1 (&fc->file);       \
                                            FCLAW_FREE (fc);                            \
                                            fclaw2d_file_error_code_v1 (errcode, cperrcode);\
                                            return NULL;}} while (0)

/** This macro prints the MPI error for sc_mpi_{read,write}.
 * This means that this macro is appropriate to call it after a non-collective
 * read or write. For a correct error handling it is required to skip the rest
 * of the non-collective code and then broadcast the error flag.
 * Can be used only multiple times in a function but will always jump to the
 * same label. This leads to correct error managing.
 */
#define FCLAW2D_FILE_CHECK_MPI_V1(errcode, user_msg) do {FCLAW2D_FILE_CHECK_VERBOSE_V1 (errcode, user_msg);\
                                                        if (!FCLAW2D_FILE_IS_SUCCESS_V1 (mpiret)) {\
                                                        goto fclaw2d_read_write_error;}} while (0)

/** Use this macro after \ref FCLAW2D_FILE_CHECK_MPI_V1 *directly* after the end of
 * non-collective statements.
 * Can be only used once in a function.
 */
/* Remark: Since we use a declaration after the label we need an empty statement. */
#define FCLAW2D_HANDLE_MPI_ERROR_V1(mpiret,fc,comm,cperrcode) do {fclaw2d_read_write_error: ;             \
                                                    int fclaw2d_mpiret_handle_error =                \
                                                    sc_MPI_Bcast (&mpiret, 1, sc_MPI_INT, 0, comm);\
                                                    SC_CHECK_MPI (fclaw2d_mpiret_handle_error);      \
                                                    *cperrcode = mpiret;                           \
                                                    if (!FCLAW2D_FILE_IS_SUCCESS_V1 (mpiret)) {         \
                                                    fclaw2d_file_error_cleanup_v1 (&fc->file);          \
                                                    FCLAW_FREE (fc);                               \
                                                    fclaw2d_file_error_code_v1 (mpiret, cperrcode);     \
                                                    return NULL;}} while (0)

/** A macro to check for file write related count errors.
 * These errors are handled as fatal errors. The macro is only applicable for
 * collective calls.
 */
#define FCLAW2D_FILE_CHECK_COUNT_V1(icount,ocount,fc,cperrcode) do { int fclaw2d_count_error_global, fclaw2d_mpiret,\
                                                 fclaw2d_rank;                                               \
                                                 int fclaw2d_file_check_count = ((int) icount != ocount);    \
                                                 fclaw2d_mpiret = sc_MPI_Allreduce (&fclaw2d_file_check_count, \
                                                 &fclaw2d_count_error_global, 1, sc_MPI_INT, sc_MPI_LOR,     \
                                                 fc->mpicomm);                                             \
                                                 SC_CHECK_MPI (fclaw2d_mpiret);                              \
                                                 fclaw2d_mpiret = sc_MPI_Comm_rank (fc->mpicomm, &fclaw2d_rank);\
                                                 SC_CHECK_MPI (fclaw2d_mpiret);                              \
                                                 *cperrcode = (fclaw2d_file_check_count) ?                   \
                                                 FCLAW2D_FILE_ERR_COUNT_V1 : sc_MPI_SUCCESS;                    \
                                                 if (fclaw2d_count_error_global)                             \
                                                 { if (fclaw2d_rank == 0) {                                  \
                                                  SC_LERRORF ("Count error at %s:%d.\n",__FILE__,          \
                                                 __LINE__);}                                               \
                                                 fclaw2d_file_error_cleanup_v1 (&fc->file);                     \
                                                 FCLAW_FREE (fc);                                          \
                                                 return NULL;}} while (0)

/** A macro to check for file write related count errors. This macro is
 * only applicable for serial calls. The errors are handled as fatal errors.
 * We assume that the macro is called on rank 0.
 */
#define FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1(icount, ocount) do {if (((int) icount) != ocount) {                        \
                                                        SC_LERRORF ("Count error on rank 0 at %s:%d.\n",__FILE__,\
                                                        __LINE__);                                               \
                                                        goto fclaw2d_write_count_error;}} while (0)

/** A macro to handle a file write error that occurred on rank 0 but need to be
 * handled collectivly. We need count_error as input since we need a variable to
 * broadcast the count error status. count_error is true if there is a count error
 * and false otherwise.
 */
/* Remark: Since we use a declaration after the label we need an empty statement. */
#define FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1(count_error,fc,cperrcode) do {fclaw2d_write_count_error: ;\
                                                    int fclaw2d_mpiret_handle = sc_MPI_Bcast (&count_error, 1, sc_MPI_INT, 0,\
                                                    fc->mpicomm);\
                                                    SC_CHECK_MPI (fclaw2d_mpiret_handle);\
                                                    *cperrcode = (count_error) ? FCLAW2D_FILE_ERR_COUNT_V1 : sc_MPI_SUCCESS;\
                                                    if (count_error) {\
                                                    fclaw2d_file_error_cleanup_v1 (&fc->file);\
                                                    FCLAW_FREE (fc);\
                                                    return NULL;}} while (0)

/** Declare closing function since it is needed before its definition. */
static int fclaw2d_file_close_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                                  int *errcode);

/** This function calculates a padding string consisting of spaces.
 * We require an already allocated array pad or NULL.
 * The number of bytes in pad must be at least divisor + 1!
 * For NULL the function calculates only the number of padding bytes.
 */
static void
fcaw2d_file_get_padding_string_v1 (size_t num_bytes, size_t divisor,
                                   char *pad, size_t *num_pad_bytes)
{
    FCLAW_ASSERT (divisor != 0 && num_pad_bytes != NULL);

    *num_pad_bytes = (divisor - (num_bytes % divisor)) % divisor;
    if (*num_pad_bytes == 0 || *num_pad_bytes == 1)
    {
        /* In these cases there is no space to add new line characters
         * but this is necessary to ensure a consistent layout in a text editor
         */
        *num_pad_bytes += divisor;
    }

    FCLAW_ASSERT (*num_pad_bytes > 1);
    FCLAW_ASSERT (*num_pad_bytes <= divisor);
    if (pad != NULL)
    {
        snprintf (pad, *num_pad_bytes + 1, "\n%-*s\n",
                  (int) *num_pad_bytes - 2, "");
    }
}

static int
fclaw2d_file_check_file_metadata_v1 (sc_MPI_Comm mpicomm,
                                     const char *filename,
                                     char
                                     user_string
                                     [FCLAW2D_FILE_USER_STRING_BYTES_V1],
                                     char *metadata,
                                     p4est_gloidx_t * global_num_quadrants)
{
    long read_global_num_quads;
    int mpiret, rank;
    int error_flag;

    FCLAW_ASSERT (metadata != NULL);

    error_flag = 0;

    mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    /* check magic number */
    if (metadata[FCLAW2D_FILE_MAGIC_BYTES_V1 - 1] != '\n')
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", "Error reading. Wrong file header format.\n");
        }
        return FCLAW2D_FILE_ERR_FORMAT_V1;
    }

    metadata[FCLAW2D_FILE_MAGIC_BYTES_V1 - 1] = '\0';
    if (strcmp (metadata, FCLAW2D_FILE_MAGIC_NUMBER_V1))
    {
        /* This check breaks if we change the magic string but this code is
         * only dedicated to implement version 1 so there should be no
         * version bump in \ref FCLAW2D_FILE_MAGIC_NUMBER_V1. */
        if (!strcmp (metadata, "p8data0") || !strcmp (metadata, "p4data0")) {
            /* check if it is just a dimension mismatch */
            return FCLAW2D_FILE_ERR_DIM_V1;
        }
        /* TODO: check for wrong endianness */
        fclaw_errorf
            ("Error reading <%s>. Wrong magic number (in file = %s, magic number = %s).\n",
             filename, metadata, FCLAW2D_FILE_MAGIC_NUMBER_V1);
        error_flag = 1;
    }

    /* check format of version string line */
    if (metadata
        [FCLAW2D_FILE_MAGIC_BYTES_V1 + FCLAW2D_FILE_VERSION_STR_BYTES_V1 -
         1] != '\n')
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", "Error reading. Wrong file header format.\n");
        }
        return FCLAW2D_FILE_ERR_FORMAT_V1;
    }

    metadata[FCLAW2D_FILE_MAGIC_BYTES_V1 + FCLAW2D_FILE_VERSION_STR_BYTES_V1 -
             1] = '\0';
    if (strlen (&metadata[FCLAW2D_FILE_MAGIC_BYTES_V1]) !=
        FCLAW2D_FILE_VERSION_STR_BYTES_V1 - 1)
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          " : Error reading. Wrong file header format.\n");
        }
        return FCLAW2D_FILE_ERR_FORMAT_V1;
    }

    /* check the format of the user string */
    if (metadata
        [FCLAW2D_FILE_MAGIC_BYTES_V1 + FCLAW2D_FILE_VERSION_STR_BYTES_V1 +
         FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1] != '\n')
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          " : Error reading. Wrong file header format.\n");
        }
        return FCLAW2D_FILE_ERR_FORMAT_V1;
    }
    /* the content of the user string is not checked */
    sc_strcopy (user_string, FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1,
                &metadata[FCLAW2D_FILE_MAGIC_BYTES_V1 +
                          FCLAW2D_FILE_VERSION_STR_BYTES_V1]);
    user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1] = '\0';

    /* check number of global quadrants */
    /* there is no \n at the end of this line */
    metadata[FCLAW2D_FILE_METADATA_BYTES_V1] = '\0';

    if (strlen
        (&metadata
         [FCLAW2D_FILE_MAGIC_BYTES_V1 + FCLAW2D_FILE_VERSION_STR_BYTES_V1 +
          FCLAW2D_FILE_USER_STRING_BYTES_V1]) != 16)
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          " : Error reading. Wrong file header format.\n");
        }
        return FCLAW2D_FILE_ERR_FORMAT_V1;
    }

    read_global_num_quads =
        sc_atol (&metadata
                 [FCLAW2D_FILE_MAGIC_BYTES_V1 +
                  FCLAW2D_FILE_VERSION_STR_BYTES_V1 +
                  FCLAW2D_FILE_USER_STRING_BYTES_V1]);
    *global_num_quadrants = (p4est_gloidx_t) read_global_num_quads;
    if (read_global_num_quads < 0)
    {
        fclaw_errorf (FCLAW2D_FILE_STRING_V1
                      " : Error reading <%s>. Negative global number of quadrants.\n",
                      filename);
        error_flag = 1;
    }

    return (error_flag) ? FCLAW2D_FILE_ERR_FORMAT_V1 : sc_MPI_SUCCESS;
}

/** Close an MPI file or its libsc-internal replacement in case of an error.
 * \param [in,out]  file    A sc_MPI_file
 * \return                  Always -1 since this function is only called
 *                          if an error already occurred.
 */
static int
fclaw2d_file_error_cleanup_v1 (sc_MPI_File * file)
{
    /* no error checking since we are called under an error condition */
    FCLAW_ASSERT (file != NULL);
#ifdef P4EST_ENABLE_MPIIO
    if (*file != sc_MPI_FILE_NULL)
    {
#else
    if ((*file)->file != sc_MPI_FILE_NULL)
    {
#endif
        /* We do not use here the libsc closing function since we do not perform
         * error checking in this function that is only called if we had already
         * an error.
         */
#ifdef P4EST_ENABLE_MPIIO
        MPI_File_close (file);
#else
        {
#ifdef P4EST_ENABLE_MPI
            int rank, mpiret;
#endif

#ifdef P4EST_ENABLE_MPI
            mpiret = sc_MPI_Comm_rank ((*file)->mpicomm, &rank);
            SC_CHECK_MPI (mpiret);

            if (rank == 0)
            {
#endif
                fclose ((*file)->file);
                (*file)->file = NULL;
#ifdef P4EST_ENABLE_MPI
            }
#endif
        }
#endif
    }
#ifndef P4EST_ENABLE_MPIIO
    SC_FREE (*file);
#endif
    return -1;
}

static int
fclaw2d_file_check_user_string (const char *user_string)
{
    int terminal_nul;
    char copied_user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1];

    /* copy user string */
    sc_strcopy (copied_user_string, FCLAW2D_FILE_USER_STRING_BYTES_V1,
                user_string);

    if (strlen (copied_user_string) < FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1)
    {
        /* user_string is nul-terminated */
        return 0;
    }

    FCLAW_ASSERT (strlen (copied_user_string) ==
                  FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1);

    /* check for nul at last byte position of user string */
    terminal_nul = user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1] == '\0';

    if (!terminal_nul)
    {
        /* user_string is not nul-terminated */
        return -1;
    }
    else
    {
        /* user_string is nul-terminated */
        return 0;
    }
}

/** Begin writing file header and saving data blocks into a parallel file.
 *
 * This function creates a new file or overwrites an existing one.
 * It is collective and creates the file on a parallel file system.
 * It takes an (optional) pointer to write a header of given size.
 * This function leaves the file open if MPI I/O is available.
 * It is necessary to call \ref
 * fclaw2d_file_close_v1 (possibly after writing one or more data sets).
 * The file is opened in a write-only mode.
 *
 * We add some basic metadata to the file.
 * The file written contains the file header and data sets
 * as specified by the open/write functions called.
 * The file header consists of the metadata specified by ForestClaw.
 *
 * The number of global quadrants must be less or equal
 * \ref FCLAW2D_FILE_MAX_GLOBAL_QUAD_V1.
 *
 * It is the application's responsibility to write sufficient header
 * information (cf. \ref fclaw2d_file_write_block_v1) to determine the number and
 * size of the data sets if such information is not recorded and maintained
 * externally.
 * However, ForestClaw makes some metadata accessible via
 * \ref fclaw2d_file_info_v1.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in] p4est          Valid forest.
 * \param [in] filename       Path to parallel file that is to be created.
 * \param [in] user_string    A user string that is written to the file header.
 *                            Only \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 -1
 *                            bytes without NUL-termination are
 *                            written to the file. If the user gives less
 *                            bytes the user_string in the file header is padded
 *                            by spaces.
 * \param [out] errcode       An errcode that can be interpreted by
 *                            \ref fclaw2d_file_error_string_v1.
 * \return                    Newly allocated context to continue writing
 *                            and eventually closing the file. NULL in
 *                            case of error.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_open_create_v1 (p4est_t * p4est, const char *filename,
                             const char *user_string, int *errcode)
{
    int mpiret, count, count_error, mpisize;
    /* We enforce the padding of the file header. */
    char metadata[FCLAW2D_FILE_METADATA_BYTES_V1 +
                  FCLAW2D_FILE_BYTE_DIV_V1 + 1];
    fclaw2d_file_context_p4est_v1_t *file_context;

    FCLAW_ASSERT (p4est_is_valid (p4est));
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (errcode != NULL);

    if (fclaw2d_file_check_user_string (user_string))
    {
        /* invalid user string */
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        /* We do not use fclaw2d file error macro since there is no
         * file context to clean up.
         */
        FCLAW2D_FILE_CHECK_VERBOSE_V1 (*errcode,
                                       FCLAW2D_FILE_STRING_V1
                                       " open_create: Invalid user string");
        return NULL;
    }

    if (!(p4est->global_num_quadrants <= FCLAW2D_FILE_MAX_GLOBAL_QUAD_V1))
    {
        /* number of global quadrant can not be written to the file header */
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        /* We do not use fclaw2d file error macro since there is no
         * file context to clean up.
         */
        FCLAW2D_FILE_CHECK_VERBOSE_V1 (*errcode,
                                       FCLAW2D_FILE_STRING_V1
                                       " open_create: Invalid number of global quadrants");
        return NULL;
    }

    file_context = FCLAW_ALLOC (fclaw2d_file_context_p4est_v1_t, 1);

    /* Open the file and create a new file if necessary */
    mpiret =
        sc_io_open (p4est->mpicomm, filename,
                    SC_IO_WRITE_CREATE, sc_MPI_INFO_NULL,
                    &file_context->file);
    FCLAW2D_FILE_CHECK_OPEN_V1 (mpiret, file_context, "File open create",
                                errcode);

    if (p4est->mpirank == 0)
    {
        /* write padded fclaw2d-defined header */
        snprintf (metadata,
                  FCLAW2D_FILE_METADATA_BYTES_V1 + FCLAW2D_FILE_BYTE_DIV_V1 +
                  1, "%.7s\n%-23s\n%-47s\n%.16lld\n%-14s\n",
                  FCLAW2D_FILE_MAGIC_NUMBER_V1, FCLAW2D_FILE_HEADER_STRING_V1,
                  user_string, (long long) p4est->global_num_quadrants, "");
        mpiret =
            sc_io_write_at (file_context->file, 0, metadata,
                            FCLAW2D_FILE_METADATA_BYTES_V1 +
                            FCLAW2D_FILE_BYTE_DIV_V1, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Writing the file header");
        count_error =
            (FCLAW2D_FILE_METADATA_BYTES_V1 + FCLAW2D_FILE_BYTE_DIV_V1 !=
             count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (FCLAW2D_FILE_METADATA_BYTES_V1 +
                                            FCLAW2D_FILE_BYTE_DIV_V1, count);
    }

    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, file_context, p4est->mpicomm,
                                 errcode);

    /* initialize the file context */
    file_context->mpicomm = p4est->mpicomm;
    file_context->local_num_quadrants = p4est->local_num_quadrants;
    file_context->global_num_quadrants = p4est->global_num_quadrants;
    mpiret = sc_MPI_Comm_size (p4est->mpicomm, &mpisize);
    file_context->global_first_quadrant =
        FCLAW_ALLOC (p4est_gloidx_t, mpisize + 1);
    memcpy (file_context->global_first_quadrant, p4est->global_first_quadrant,
            (mpisize + 1) * sizeof (p4est_gloidx_t));
    file_context->gfq_owned = 1;

    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, file_context, errcode);

    file_context->accessed_bytes = 0;
    file_context->num_calls = 0;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return file_context;
}

/** Open a file for reading without knowing the p4est that is associated
 * with the mesh-related data in the file (cf. \ref fclaw2d_file_open_read_v1).
 * For more general comments on open_read see the documentation of
 * \ref fclaw2d_file_open_read_v1.
 * The parameters that are not documented are the same as in \ref
 * fclaw2d_file_open_read_v1.
 *
 * \param [in]  mpicomm   The MPI communicator that is used to read the file.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_open_read_ext_v1 (sc_MPI_Comm mpicomm, const char *filename,
                               char *user_string,
                               p4est_gloidx_t * global_num_quadrants,
                               int *errcode)
{
    int mpiret, rank;
    int count, count_error;
    char metadata[FCLAW2D_FILE_METADATA_BYTES_V1 + 1];
    fclaw2d_file_context_p4est_v1_t *file_context =
        FCLAW_ALLOC (fclaw2d_file_context_p4est_v1_t, 1);

    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (global_num_quadrants != NULL);
    FCLAW_ASSERT (errcode != NULL);

    /* Open the file in the reading mode */
    mpiret =
        sc_io_open (mpicomm, filename, SC_IO_READ,
                    sc_MPI_INFO_NULL, &file_context->file);
    FCLAW2D_FILE_CHECK_OPEN_V1 (mpiret, file_context, "File open read",
                                errcode);

    file_context->mpicomm = mpicomm;
    file_context->local_num_quadrants = 0;      /* not set for read calls */
    file_context->global_first_quadrant = NULL;
    file_context->gfq_owned = 0;
    file_context->accessed_bytes = 0;
    file_context->num_calls = 0;

    /* get the MPI rank */
    mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    /* read metadata and deallocate in case of error */
    if (rank == 0)
    {
        /* read metadata on rank 0 */
        mpiret =
            sc_io_read_at (file_context->file, 0, metadata,
                           FCLAW2D_FILE_METADATA_BYTES_V1, sc_MPI_BYTE,
                           &count);

        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Reading metadata");
        count_error = (FCLAW2D_FILE_METADATA_BYTES_V1 != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (FCLAW2D_FILE_METADATA_BYTES_V1,
                                            count);

        metadata[FCLAW2D_FILE_METADATA_BYTES_V1] = '\0';
        /* parse metadata; we do not use file_info because we do not want a Bcast */
        mpiret =
            fclaw2d_file_check_file_metadata_v1 (mpicomm, filename,
                                                 user_string, metadata,
                                                 global_num_quadrants);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Check file header");
    }

    /* error checking */
    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, file_context, mpicomm, errcode);
    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, file_context, errcode);

    /* broadcast the user string of the file */
    mpiret =
        sc_MPI_Bcast (user_string, FCLAW2D_FILE_USER_STRING_BYTES_V1,
                      sc_MPI_BYTE, 0, mpicomm);
    SC_CHECK_MPI (mpiret);

    /* broadcast the number of global quadrants */
    mpiret =
        sc_MPI_Bcast (global_num_quadrants, sizeof (p4est_gloidx_t),
                      sc_MPI_BYTE, 0, mpicomm);
    SC_CHECK_MPI (mpiret);

    file_context->global_num_quadrants = *global_num_quadrants;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return file_context;
}

/** currently unused */
#if 0
/** Open a file for reading and read its user string on rank zero.
 * The user string is broadcasted to all ranks after reading.
 * The file must exist and be at least of the size of the file header.
 *
 * If the file has wrong metadata the function reports the error using
 * \ref fclaw_errorf, collectively close the file and deallocate
 * the file context. In this case the function returns NULL on all ranks.
 * The wrong file format or a wrong file header causes \ref FCLAW2D_FILE_ERR_FORMAT_V1
 * as errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in] p4est            The forest must be of the same refinement
 *                              pattern as the one used for writing the file.
 *                              Its global number of quadrants must match.
 *                              It is possible, however, to use a different
 *                              partition or number of ranks from writing it.
 * \param [in] filename         The path to the file that is opened.
 * \param [in,out] user_string  At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1
 *                              bytes. The user string is written
 *                              to the passed array including padding spaces
 *                              and a trailing NUL-termination.
 * \param [out] errcode         An errcode that can be interpreted by \ref
 *                              fclaw2d_file_error_string_v1.
 * \return                      Newly allocated context to continue reading
 *                              and eventually closing the file. NULL in
 *                              case of error.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_open_read_v1 (p4est_t * p4est, const char *filename,
                           char *user_string, int *errcode)
{
    p4est_gloidx_t global_num_quadrants;
    fclaw2d_file_context_p4est_v1_t *fc;

    FCLAW_ASSERT (p4est_is_valid (p4est));
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (errcode != NULL);

    fc = fclaw2d_file_open_read_ext_v1 (p4est->mpicomm, filename, user_string,
                                        &global_num_quadrants, errcode);

    /* check global number of quadrants */
    if (fc != NULL && p4est->global_num_quadrants != global_num_quadrants)
    {
        if (p4est->mpirank == 0)
        {
            fclaw_errorf (FCLAW2D_FILE_STRING_V1
                          " open_read: global number of "
                          "quadrants mismatch (in file = %lld,"
                          " by parameter = %lld)\n",
                          (long long) global_num_quadrants,
                          (long long) p4est->global_num_quadrants);
        }
        fclaw2d_file_close_v1 (fc, errcode);
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " open_read: close file", errcode);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return NULL;
    }

    if (fc != NULL)
    {
        /* use the partition of the given p4est */
        fc->global_first_quadrant = p4est->global_first_quadrant;
        fc->gfq_owned = 0;
    }

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return fc;
}
#endif

/** Write a block section to an opened file.
 * This function requires an opened file context.
 * The block data and its metadata are written on rank 0.
 * The number of block bytes must be less or equal
 * \ref FCLAW2D_FILE_MAX_BLOCK_SIZE_V1.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [out] fc            Context previously created by \ref
 *                            fclaw2d_file_open_create_v1.
 * \param [in]  block_size    The size of block in bytes.
 *                            May be equal to 0. In this case the
 *                            section header and the padding
 *                            is still written.
 *                            This function returns the passed fc
 *                            parameter and sets errcode to
 *                            \ref FCLAW2D_FILE_ERR_SUCCESS_V1 if it is called
 *                            for block_size == 0.
 * \param [in]  block_data    A sc_array with one element and element size
 *                            equal to \a block_size.
 *                            The array points to the block data. The user is
 *                            responsible for the validality of the block
 *                            data. block_data can be NULL if
 *                            block_size == 0.
 * \param [in]  user_string   Maximal \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                            These chars are written to the block
 *                            header and padded to
 *                            \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1 chars
 *                            by adding spaces. The '\0' is not written
 *                            to the file.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return the input context to continue writing
 *                            and eventually closing the file. The return
 *                            value is NULL in case of error, then
 *                            it also holds errcode != 0 and the file is
 *                            tried to close and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_write_block_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                             size_t block_size, sc_array_t * block_data,
                             const char *user_string, int *errcode)
{
    size_t num_pad_bytes;
    char header_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + 1],
        pad[FCLAW2D_FILE_MAX_NUM_PAD_BYTES_V1];
    int mpiret, count, count_error, rank;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (fc->global_first_quadrant != NULL);
    FCLAW_ASSERT (block_data != NULL);
    FCLAW_ASSERT (block_size == 0 || block_data->array != NULL);
    FCLAW_ASSERT (block_size == block_data->elem_size);
    FCLAW_ASSERT (block_data->elem_count == 1);
    FCLAW_ASSERT (errcode != NULL);

    if (fclaw2d_file_check_user_string (user_string))
    {
        /* invalid user string */
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " write_block: Invalid user string",
                                    errcode);
    }

    if (!(block_size <= FCLAW2D_FILE_MAX_BLOCK_SIZE_V1))
    {
        /* invalid header size */
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " write_block: Invalid block size",
                                    errcode);
    }

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

#ifdef P4EST_ENABLE_MPIIO
    /* set the file size */
    mpiret = MPI_File_set_size (fc->file,
                                FCLAW2D_FILE_METADATA_BYTES_V1 +
                                FCLAW2D_FILE_BYTE_DIV_V1 + block_size +
                                FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                fc->accessed_bytes);
    FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Set file size", errcode);
#else
    /* We do not perform this optimization without MPI I/O */
#endif

    num_pad_bytes = 0;
    if (rank == 0)
    {
        /* header-dependent metadata */
        snprintf (header_metadata,
                  FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                  1, "B %.13llu\n%-47s\n", (unsigned long long) block_size,
                  user_string);

        /* write header-dependent metadata */
        mpiret =
            sc_io_write_at (fc->file,
                            fc->accessed_bytes +
                            FCLAW2D_FILE_METADATA_BYTES_V1 +
                            FCLAW2D_FILE_BYTE_DIV_V1, header_metadata,
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1, sc_MPI_BYTE,
                            &count);

        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Writing header metadata");
        count_error = (FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1
            (FCLAW2D_FILE_FIELD_HEADER_BYTES_V1, count);
    }
    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, fc, fc->mpicomm, errcode);
    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, fc, errcode);

    /*write header data */
    if (rank == 0)
    {
        mpiret =
            sc_io_write_at (fc->file,
                            fc->accessed_bytes +
                            FCLAW2D_FILE_METADATA_BYTES_V1 +
                            FCLAW2D_FILE_BYTE_DIV_V1 +
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1,
                            block_data->array, block_size, sc_MPI_BYTE,
                            &count);

        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Writing block data");
        count_error = ((int) block_size != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (block_size, count);

        /* write padding bytes */
        fcaw2d_file_get_padding_string_v1 (block_size,
                                           FCLAW2D_FILE_BYTE_DIV_V1, pad,
                                           &num_pad_bytes);
        mpiret =
            sc_io_write_at (fc->file,
                            fc->accessed_bytes +
                            FCLAW2D_FILE_METADATA_BYTES_V1 +
                            FCLAW2D_FILE_BYTE_DIV_V1 +
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + block_size,
                            pad, num_pad_bytes, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret,
                                   "Writing padding bytes for header data");
        count_error = ((int) num_pad_bytes != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (num_pad_bytes, count);
    }
    else
    {
        fcaw2d_file_get_padding_string_v1 (block_size,
                                           FCLAW2D_FILE_BYTE_DIV_V1, NULL,
                                           &num_pad_bytes);
    }

    /* This is *not* the processor local value */
    fc->accessed_bytes +=
        block_size + FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + num_pad_bytes;
    ++fc->num_calls;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return fc;
}

/** Collectivly read and check block metadata.
 * If user_string == NULL data_size is not compared to
 * read_data_size.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_block_metadata_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                                     size_t *read_data_size, size_t data_size,
                                     char block_type,
                                     char *user_string, int *errcode)
{
    int mpiret, count, count_error, rank;
    int bytes_to_read;
    int err_flag, invalid_block;
    char block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1];
    size_t data_block_size, num_pad_bytes;

    FCLAW_ASSERT (read_data_size != NULL);
    FCLAW_ASSERT (errcode != NULL);

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    bytes_to_read =
        (user_string !=
         NULL) ? FCLAW2D_FILE_FIELD_HEADER_BYTES_V1
        : (FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + 2);
    if (rank == 0)
    {
        mpiret = sc_io_read_at (fc->file,
                                fc->accessed_bytes +
                                FCLAW2D_FILE_METADATA_BYTES_V1 +
                                FCLAW2D_FILE_BYTE_DIV_V1, block_metadata,
                                bytes_to_read, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret,
                                   "Reading data section-wise metadata");
        count_error = (bytes_to_read != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (bytes_to_read, count);
    }
    /* In the case of error the return value is still NULL */
    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, fc, fc->mpicomm, errcode);
    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, fc, errcode);

    /* broadcast block metadata to calculate correct internals on each rank */
    mpiret =
        sc_MPI_Bcast (block_metadata, bytes_to_read, sc_MPI_BYTE, 0,
                      fc->mpicomm);
    SC_CHECK_MPI (mpiret);

    /* check for given block specifying character */
    invalid_block = 0;
    if (block_metadata[0] != block_type)
    {
        invalid_block = block_metadata[0] != 'F' && block_metadata[0] != 'B';
        if (rank == 0)
        {
            if (invalid_block)
            {
                fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                              ": Error reading. Invalid data section type.\n");
            }
            else
            {
                fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                              ": Error reading. Wrong data section type.\n");
            }
        }
        fclaw2d_file_error_cleanup_v1 (&fc->file);
        FCLAW_FREE (fc);
        if (invalid_block)
        {
            *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        }
        else
        {
            *errcode = FCLAW2D_FILE_ERR_SECTION_TYPE_V1;
        }
        return NULL;
    }

    /* check '\n' to check the format */
    if (block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + 1] != '\n')
    {
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          ": Error reading. Wrong section header format.\n");
        }
        fclaw2d_file_error_cleanup_v1 (&fc->file);
        FCLAW_FREE (fc);
        *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        return NULL;
    }

    /* process the block metadata */
    block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + 1] = '\0';
    /* we cut off the block type specifier */
    *read_data_size = sc_atol (&block_metadata[2]);

    if (user_string != NULL && *read_data_size != data_size)
    {
        if (rank == 0)
        {
            fclaw_errorf (FCLAW2D_FILE_STRING_V1
                          ": Error reading. Wrong section data size (in file = %ld, by parameter = %ld).\n",
                          *read_data_size, data_size);
        }
        fclaw2d_file_error_cleanup_v1 (&fc->file);
        FCLAW_FREE (fc);
        *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        return NULL;
    }

    if (user_string != NULL)
    {
        /* check '\n' to check the format */
        if (block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 - 1] != '\n')
        {
            if (rank == 0)
            {
                fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                              ": Error reading. Wrong section header format.\n");
            }
            fclaw2d_file_error_cleanup_v1 (&fc->file);
            FCLAW_FREE (fc);
            *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
            return NULL;
        }

        /* null-terminate the user string of the current block */
        block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 - 1] = '\0';

        /* copy the user string, '\0' was already set above */
        sc_strcopy (user_string, FCLAW2D_FILE_USER_STRING_BYTES_V1,
                    &block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 +
                                    2]);
        FCLAW_ASSERT (user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1] ==
                      '\0');
    }

    /* check the padding structure */
    err_flag = 0;
    if (rank == 0)
    {
        /* calculate number of padding bytes */
        if (block_metadata[0] == 'F')
        {
            data_block_size = *read_data_size * fc->global_num_quadrants;
        }
        else if (block_metadata[0] == 'B')
        {
            data_block_size = *read_data_size;
        }
        else
        {
            /* We assume that this function is called for a valid block type. */
            SC_ABORT_NOT_REACHED ();
        }
        fcaw2d_file_get_padding_string_v1 (data_block_size,
                                           FCLAW2D_FILE_BYTE_DIV_V1, NULL,
                                           &num_pad_bytes);
        /* read padding bytes */
        mpiret = sc_io_read_at (fc->file,
                                fc->accessed_bytes +
                                FCLAW2D_FILE_METADATA_BYTES_V1 +
                                FCLAW2D_FILE_BYTE_DIV_V1 +
                                FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                data_block_size, block_metadata,
                                num_pad_bytes, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Reading padding bytes");
        count_error = ((int) num_pad_bytes != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (num_pad_bytes, count);
        /* check '\n' in padding bytes */
        if (block_metadata[0] != '\n'
            || block_metadata[num_pad_bytes - 1] != '\n')
        {
            err_flag = 1;
        }
    }
    /* broadcast error status */
    mpiret = sc_MPI_Bcast (&err_flag, 1, sc_MPI_INT, 0, fc->mpicomm);
    SC_CHECK_MPI (mpiret);

    if (err_flag)
    {
        /* wrong padding format */
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          ": Error reading. Wrong padding format.\n");
        }
        fclaw2d_file_error_cleanup_v1 (&fc->file);
        FCLAW_FREE (fc);
        *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        return NULL;
    }

    return fc;
}

/** Read a header block from an opened file.
 * This function requires an opened file context.
 * The header data is read on rank 0.
 *
 * If the user does not have the header_size to call this function, the user
 * can user \ref fclaw2d_file_info_v1 to obtain the required information.
 *
 * The passed header_size is compared to the header_size stored in the file.
 * If the values do not equal each other, the function reports details via
 * \ref fclaw_errorf and closes and deallocate the file context. The return
 * value in this case is NULL.
 * If the block header information is not matching the passed parameters
 * the function sets \ref FCLAW2D_FILE_ERR_FORMAT_V1 for errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [out] fc              Context previously created by \ref
 *                              fclaw2d_file_open_create_v1.
 * \param [in]  header_size     The size of the header that is read.
 * \param [in, out] header_data \a header_size allocated bytes in an sc_array
 *                              with one element and \a header_size as element
 *                              size. This data will be filled with the header
 *                              data from file. If this is NULL it means that
 *                              the current header block is skipped and the
 *                              internal file pointer of the file context is
 *                              set to the next data block. If current data
 *                              block is not a header block, the file is closed
 *                              and the file context is deallocated. Furthermore,
 *                              in this case the function returns NULL and sets
 *                              errcode to \ref FCLAW2D_FILE_ERR_FORMAT_V1. In case
 *                              of skipping the header section \a header_size
 *                              needs also to coincide with the header size
 *                              given in the file.
 * \param [in,out] user_string  At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                              Filled by the padded user string and
 *                              a trailing NUL-termination char.
 * \param [out] errcode         An errcode that can be interpreted by \ref
 *                              fclaw2d_file_error_string_v1.
 * \return                      Return the input context to continue reading
 *                              and eventually closing the file. The return value
 *                              is NULL if the function was called for
 *                              header_size == 0. The return
 *                              value is also NULL in case of error but then
 *                              it also holds errcode != 0 and the file is
 *                              tried to close and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_block_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                            size_t block_size, sc_array_t * block_data,
                            char *user_string, int *errcode)
{
    int mpiret, count, count_error, rank;
    size_t num_pad_bytes, read_data_size;
#ifdef P4EST_ENABLE_MPIIO
    sc_MPI_Offset size;
#endif

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (block_data == NULL || block_size == block_data->elem_size);
    FCLAW_ASSERT (block_data == NULL || block_data->elem_count == 1);
    FCLAW_ASSERT (errcode != NULL);
    FCLAW_ASSERT (user_string != NULL);

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    num_pad_bytes = 0;
    /* calculate the padding bytes for this header data */
    fcaw2d_file_get_padding_string_v1 (block_size, FCLAW2D_FILE_BYTE_DIV_V1,
                                       NULL, &num_pad_bytes);
    if (block_data == NULL)
    {
        /* Nothing to read but we shift our own file pointer */
        if (fclaw2d_file_read_block_metadata_v1
            (fc, &read_data_size, block_size, 'B', user_string,
             errcode) == NULL)
        {
            fclaw2d_file_error_code_v1 (*errcode, errcode);
            return NULL;
        }

        fc->accessed_bytes +=
            block_size + FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + num_pad_bytes;
        ++fc->num_calls;
        *errcode = sc_MPI_SUCCESS;
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return fc;
    }

    FCLAW_ASSERT (block_size == block_data->elem_size);

#ifdef P4EST_ENABLE_MPIIO
    /* check file size; no sync required because the file size does not
     * change in the reading mode.
     */
    mpiret = MPI_File_get_size (fc->file, &size);
    FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Get file size for read",
                                errcode);
    if (((size_t) size) - FCLAW2D_FILE_METADATA_BYTES_V1 -
        FCLAW2D_FILE_BYTE_DIV_V1 - FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 <
        block_size)
    {
        /* report wrong file size, collectively close the file and deallocate fc */
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          ": Error reading. File has less bytes than the user wants to read.\n");
        }
        mpiret = fclaw2d_file_close_v1 (fc, &mpiret);
        FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " read_data: close file", errcode);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return NULL;
    }
#else
    /* There is no C-standard functionality to get the file size */
#endif

    /* check the header metadata */
    if (fclaw2d_file_read_block_metadata_v1
        (fc, &read_data_size, block_size, 'B', user_string, errcode) == NULL)
    {
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return NULL;
    }

    if (rank == 0)
    {
        mpiret = sc_io_read_at (fc->file, fc->accessed_bytes +
                                FCLAW2D_FILE_METADATA_BYTES_V1 +
                                FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                FCLAW2D_FILE_BYTE_DIV_V1, block_data->array,
                                (int) block_size, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Reading header data");
        count_error = ((int) block_size != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (block_size, count);
    }
    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, fc, fc->mpicomm, errcode);
    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, fc, errcode);

    /* broadcast the header data */
    mpiret =
        sc_MPI_Bcast (block_data->array, block_size, sc_MPI_BYTE, 0,
                      fc->mpicomm);
    SC_CHECK_MPI (mpiret);

    fc->accessed_bytes +=
        block_size + FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + num_pad_bytes;
    ++fc->num_calls;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return fc;
}

/** Write one (more) per-quadrant data set to a parallel output file.
 *
 * This function requires an opened file context.
 * The data set is appended to the header/previously written data sets.
 * This function writes a block of the size number of quadrants * data_size.
 *
 * The number of bytes per field entry must be less or equal
 * \ref FCLAW2D_FILE_MAX_FIELD_ENTRY_SIZE_V1.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [out] fc            Context previously created by \ref
 *                            fclaw2d_file_open_create_v1.
 * \param [in] quadrant_size  The number of bytes per quadrant. This number
 *                            must coincide with \a quadrant_data->elem_size.
 * \param [in] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            written per quadrant. The quadrant data is expected
 *                            to be stored according to the Morton order of
 *                            the quadrants. For quadrant_data->elem_size == 0
 *                            the function writes an empty field. The section
 *                            header and the padding is still written.
 *                            In this case errcode is set
 *                            to \ref FCLAW2D_FILE_ERR_SUCCESS_V1.
 * \param [in] user_string    An array of maximal \ref
 *                            FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes that
 *                            is written without the NUL-termination
 *                            after the array-dependent metadata and before
 *                            the actual data. If the char array is shorter the
 *                            written char array will be padded to the
 *                            right by spaces. The user_string is
 *                            written on rank 0 and therefore also only
 *                            required on rank 0. Can be NULL for other
 *                            ranks.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return the input context to continue writing
 *                            and eventually closing the file. The return value
 *                            is NULL if the function was called for
 *                            quadrant_data->elem_size == 0. The return
 *                            value is also NULL in case of error but then
 *                            it also holds errcode != 0 and the file is
 *                            tried to close and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_write_field_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                             size_t quadrant_size, sc_array_t * quadrant_data,
                             const char *user_string, int *errcode)
{
    size_t bytes_to_write, num_pad_bytes, array_size;
    char array_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + 1],
        pad[FCLAW2D_FILE_MAX_NUM_PAD_BYTES_V1];
    sc_MPI_Offset write_offset;
    int mpiret, count, count_error, rank;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (quadrant_data != NULL
                  && (quadrant_data->elem_count == 0
                      || quadrant_data->elem_count ==
                      (size_t) fc->local_num_quadrants));
    FCLAW_ASSERT (quadrant_size == quadrant_data->elem_size);
    FCLAW_ASSERT (errcode != NULL);

    if (fclaw2d_file_check_user_string (user_string))
    {
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " write_field: Invalid user string",
                                    errcode);
    }

    if (!(quadrant_data->elem_size <= FCLAW2D_FILE_MAX_FIELD_ENTRY_SIZE_V1))
    {
        *errcode = FCLAW2D_FILE_ERR_IN_DATA_V1;
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " write_field: Invalid byte number per field entry",
                                    errcode);
    }

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    /* Check how many bytes we write to the disk */
    bytes_to_write = quadrant_data->elem_count * quadrant_data->elem_size;

    /* rank-dependent byte offset */
    write_offset = FCLAW2D_FILE_METADATA_BYTES_V1 + FCLAW2D_FILE_BYTE_DIV_V1 +
        fc->global_first_quadrant[rank] * quadrant_data->elem_size;

#ifdef P4EST_ENABLE_MPIIO
    /* set the file size */
    mpiret = MPI_File_set_size (fc->file,
                                FCLAW2D_FILE_METADATA_BYTES_V1 +
                                FCLAW2D_FILE_BYTE_DIV_V1 +
                                fc->global_num_quadrants *
                                quadrant_data->elem_size +
                                FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                fc->accessed_bytes);
    FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Set file size", errcode);
#else
    /* We do not perform this optimization without MPI I/O */
#endif

    num_pad_bytes = 0;
    if (rank == 0)
    {
        /* array-dependent metadata */
        snprintf (array_metadata,
                  FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                  1, "F %.13llu\n%-47s\n",
                  (unsigned long long) quadrant_data->elem_size, user_string);

        /* write array-dependent metadata */
        mpiret =
            sc_io_write_at (fc->file, fc->accessed_bytes + write_offset,
                            array_metadata,
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1, sc_MPI_BYTE,
                            &count);

        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret, "Writing array metadata");
        count_error = (FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1
            (FCLAW2D_FILE_FIELD_HEADER_BYTES_V1, count);
    }
    FCLAW2D_HANDLE_MPI_ERROR_V1 (mpiret, fc, fc->mpicomm, errcode);
    FCLAW2D_HANDLE_MPI_COUNT_ERROR_V1 (count_error, fc, errcode);

    /* write array data */
    mpiret =
        sc_io_write_at_all (fc->file,
                            fc->accessed_bytes + write_offset +
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1,
                            quadrant_data->array, bytes_to_write, sc_MPI_BYTE,
                            &count);
    FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Writing quadrant-wise", errcode);
    FCLAW2D_FILE_CHECK_COUNT_V1 (bytes_to_write, count, fc, errcode);

  /** We place the padding bytes write here because for the sequential
   * IO operations the order of fwrite calls plays a role.
   */
    /* write padding bytes */
    if (rank == 0)
    {
        /* Calculate and write padding bytes for array data */
        array_size = fc->global_num_quadrants * quadrant_data->elem_size;
        fcaw2d_file_get_padding_string_v1 (array_size,
                                           FCLAW2D_FILE_BYTE_DIV_V1, pad,
                                           &num_pad_bytes);

        mpiret =
            sc_io_write_at (fc->file,
                            fc->accessed_bytes +
                            FCLAW2D_FILE_METADATA_BYTES_V1 +
                            FCLAW2D_FILE_BYTE_DIV_V1 + array_size +
                            FCLAW2D_FILE_FIELD_HEADER_BYTES_V1, pad,
                            num_pad_bytes, sc_MPI_BYTE, &count);
        FCLAW2D_FILE_CHECK_MPI_V1 (mpiret,
                                   "Writing padding bytes for a data array");
        /* We do not need to call FCLAW2D_FILE_HANDLE_MPI_ERROR in the next
         * collective line of code since FCLAW2D_FILE_HANDLE_MPI_ERROR was already
         * called in this scope.
         */
        count_error = ((int) num_pad_bytes != count);
        FCLAW2D_FILE_CHECK_COUNT_SERIAL_V1 (num_pad_bytes, count);
    }
    else
    {
        array_size = fc->global_num_quadrants * quadrant_data->elem_size;
        fcaw2d_file_get_padding_string_v1 (array_size,
                                           FCLAW2D_FILE_BYTE_DIV_V1, NULL,
                                           &num_pad_bytes);
    }

    /* This is *not* the processor local value */
    fc->accessed_bytes +=
        quadrant_data->elem_size * fc->global_num_quadrants +
        FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + num_pad_bytes;
    ++fc->num_calls;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return fc;
}

/** Read a data field and specify the partition for reading in parallel.
 * See also the documentation of \ref fclaw2d_file_read_field_v1.
 *
 * \param [in]  gfq   An array of the size mpisize + 1 that contains the global
 *                    first quadrants per rank and
 *                    gfq[mpisize] == global_num_quadrants. This defines
 *                    partition that is used to read the data field in parallel.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_field_ext_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                                p4est_gloidx_t * gfq, size_t quadrant_size,
                                sc_array_t * quadrant_data, char *user_string,
                                int *errcode)
{
    int count;
    size_t bytes_to_read, num_pad_bytes, array_size, read_data_size;
#ifdef P4EST_ENABLE_MPIIO
    sc_MPI_Offset size;
#endif
    int mpiret, rank, mpisize;

    FCLAW_ASSERT (fc != NULL);

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (fc->mpicomm, &mpisize);
    SC_CHECK_MPI (mpiret);

    FCLAW_ASSERT (gfq != NULL);
    FCLAW_ASSERT (errcode != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (quadrant_data == NULL
                  || quadrant_size == quadrant_data->elem_size);

    /* check gfq in the debug mode */
    FCLAW_ASSERT (gfq[0] == 0);
    FCLAW_ASSERT (gfq[mpisize] == fc->global_num_quadrants);

    if (quadrant_data != NULL)
    {
        sc_array_resize (quadrant_data, (size_t) (gfq[rank + 1] - gfq[rank]));
    }

    /* check how many bytes we read from the disk */
    bytes_to_read = ((size_t) (gfq[rank + 1] - gfq[rank])) * quadrant_size;

#ifdef P4EST_ENABLE_MPIIO
    /* check file size; no sync required because the file size does not
     * change in the reading mode.
     */
    mpiret = MPI_File_get_size (fc->file, &size);
    FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Get file size for read",
                                errcode);
    if (((size_t) size) - FCLAW2D_FILE_METADATA_BYTES_V1 -
        FCLAW2D_FILE_BYTE_DIV_V1 - FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 <
        bytes_to_read)
    {
        /* report wrong file size, collectively close the file and deallocate fc */
        if (rank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          ": Error reading. File has less bytes than the user wants to read.\n");
        }
        mpiret = fclaw2d_file_close_v1 (fc, &mpiret);
        FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " read_data: close file", errcode);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return NULL;
    }
#else
    /* There is no C-standard functionality to get the file size */
#endif

    /* check the array metadata */
    if (fclaw2d_file_read_block_metadata_v1
        (fc, &read_data_size, quadrant_size, 'F', user_string,
         errcode) == NULL)
    {
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return NULL;
    }

    /* calculate the padding bytes for this data array */
    array_size = fc->global_num_quadrants * quadrant_size;
    fcaw2d_file_get_padding_string_v1 (array_size, FCLAW2D_FILE_BYTE_DIV_V1,
                                       NULL, &num_pad_bytes);

    if (quadrant_data != NULL)
    {
        mpiret = sc_io_read_at_all (fc->file,
                                    fc->accessed_bytes +
                                    FCLAW2D_FILE_METADATA_BYTES_V1 +
                                    FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                    FCLAW2D_FILE_BYTE_DIV_V1 + gfq[rank]
                                    * quadrant_data->elem_size,
                                    quadrant_data->array, bytes_to_read,
                                    sc_MPI_BYTE, &count);

        FCLAW2D_FILE_CHECK_NULL_V1 (mpiret, fc, "Reading quadrant-wise",
                                    errcode);
        FCLAW2D_FILE_CHECK_COUNT_V1 (bytes_to_read, count, fc, errcode);
    }

    fc->accessed_bytes +=
        quadrant_size * fc->global_num_quadrants +
        FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + num_pad_bytes;
    ++fc->num_calls;

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return fc;
}

/** Read one (more) per-quadrant data set from a parallel input file.
 * This function requires an opened file context.
 * This function requires the appropriate number of readable bytes.
 * In practice, the data size to read should match the size written.
 * This function reports an error if the number of bytes to read is
 * bigger than the dataset that corresponds to the processor.
 * The data size to read is encoded by the element size of quadrant_data
 * It is legal to close a file before all data sets have been read.
 *
 * The function closes and deallocates the file context and returns NULL
 * if the bytes the user wants to read exceed the given file and/or
 * the element size of the array given by quadrant_data->elem_size does not
 * coincide with the element size according to the array metadata given in
 * the file.
 *
 * If the block header information is not matching the passed parameters
 * the function sets \ref FCLAW2D_FILE_ERR_FORMAT_V1 for errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            fclaw2d_file_open_read_v1 (_ext).  It keeps track
 *                            of the data sets read one after another.
 * \param [in] quadrant_size  The number of bytes per quadrant. This number
 *                            must coincide with \a quadrant_data->elem_size.
 * \param [in,out] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            read per quadrant. The quadrant data is read
 *                            according to the Morton order of the quadrants.
 *                            \a quadrant_data->elem_size must coincide with
 *                            the section data size in the file.
 *                            quadrant_data == NULL means that the data is
 *                            skipped and the internal file pointer is incremented.
 *                            In the case of skipping \a quadrant_size is still
 *                            checked using the corresponding value read from
 *                            the file. If fc was opened by
 *                            \ref fclaw2d_file_open_read_ext_v1 and
 *                            \a fc->global_first_quadrant was not set by the
 *                            user, the function uses a uniform partition to read
 *                            the data field in parallel.
 *                            \a quadrant_data is resized by \ref sc_array_resize.
 * \param [in,out]  user_string At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                            The user string is read on rank 0 and internally
 *                            broadcasted to all ranks.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program or if
 *                            the function was called with quadrant_data == NULL.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_field_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                            size_t quadrant_size, sc_array_t * quadrant_data,
                            char *user_string, int *errcode)
{
    int mpiret, mpisize, rank;
    p4est_gloidx_t *gfq = NULL;
    fclaw2d_file_context_p4est_v1_t *retfc;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (errcode != NULL);
    FCLAW_ASSERT (quadrant_data == NULL
                  || quadrant_size == quadrant_data->elem_size);

    /* If this function is used on a file context obtained by
     * \ref fclaw2d_file_open_read_v1 the global_first_quadrant
     * array is set to the corresponding one of the p4est
     * given in the call \ref fclaw2d_file_open_read_v1.
     * Otherwise the file context was obtained by calling
     * \ref fclaw2d_file_open_read_ext_v1. This means the
     * global_first_quadrant array is not set since
     * there is no given p4est. In this case we
     * compute a uniform partition.
     */
    if (fc->global_first_quadrant == NULL)
    {
        /* there is no partition set in the file context */
        mpiret = sc_MPI_Comm_size (fc->mpicomm, &mpisize);
        SC_CHECK_MPI (mpiret);

        gfq = FCLAW_ALLOC (p4est_gloidx_t, mpisize + 1);

        /* calculate gfq for a uniform partition */
        p4est_comm_global_first_quadrant (fc->global_num_quadrants, mpisize,
                                          gfq);
    }
    else
    {
        gfq = fc->global_first_quadrant;
    }

    mpiret = sc_MPI_Comm_rank (fc->mpicomm, &rank);
    SC_CHECK_MPI (mpiret);

    if (quadrant_data != NULL)
    {
        /* allocate the memory for the quadrant data */
        sc_array_resize (quadrant_data, (size_t) (gfq[rank + 1] - gfq[rank]));
    }

    retfc =
        fclaw2d_file_read_field_ext_v1 (fc, gfq, quadrant_size, quadrant_data,
                                        user_string, errcode);
    if (fc->global_first_quadrant == NULL)
    {
        FCLAW_FREE (gfq);
    }

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return retfc;
}

/** A data type that encodes the metadata of one data block in a fclaw2d data file.
 */
typedef struct fclaw2d_file_section_metadata_v1
{
    char block_type;              /**< 'H' (header) or 'F' (data file) */
    size_t data_size;             /**< data size in bytes per array element ('F')
                                       or of the header section ('H') */
    char user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1];              /**< user string of the data section */
}
fclaw2d_file_section_metadata_v1_t;

/** currently unused */
#if 0
/** Read metadata information of a file written by a matching forest.
 * Matching refers to the global count of quadrants; partition is irrelevant.
 *
 * This function parses the given file on rank 0 and broadcasts the information
 * on the number of data fields contained to all other ranks.  Collective call.
 *
 * This function catches all I/O and file format errors and returns a valid MPI
 * error class related to file handling.  Errors are collectively synchronized.
 *
 * If the number of bytes that the user intend to read is larger than the number
 * bytes left in the file, the function prints out an information about this
 * situation using \ref fclaw_errorf. In this case the function reads the bytes
 * that are possible to read but returns NULL to indicate an error.
 * If the file or block header information is not matching the passed parameters
 * the function sets \ref FCLAW2D_FILE_ERR_FORMAT_V1 for errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  p4est               A p4est that is only required for the
 *                                  MPI communicator, and to verify the
 *                                  global quadrant count found in the file.
 * \param [in]  filename            Path to parallel file.
 * \param [in,out] user_string      At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1
 *                                  bytes. This array will be filled with the
 *                                  user string of the file after a successful
 *                                  call of this function.
 * \param [in,out] data_sections    After a successful function call this
 *                                  variable holds an array with a length
 *                                  corresponding to the number of data section
 *                                  in the file that are successfully found
 *                                  and seeked. The values in the array are the
 *                                  number of bytes of stored data per quadrant.
 *                                  Require elem_size->elem_size
 *                                  == sizeof (fclaw2d_file_section_metadata_v1_t)
 *                                  on input and preserve it on output.
 *                                  See fclaw2d_file_section_metadata_v1_t to obtain
 *                                  detailed information about the data blocks
 *                                  of the file.
 * \param [out] errcode             An errcode that can be interpreted by \ref
 *                                  fclaw2d_file_error_string_v1.
 * \return                          0 for a successful call and -1 in case of
 *                                  an error. See also errcode argument.
 */
static int
fclaw2d_file_info_v1 (p4est_t * p4est, const char *filename,
                      char *user_string, sc_array_t * data_sections,
                      int *errcode)
{
    int mpiret, eclass;
    int retval;
    int count, count_error;
    long long_header;
    size_t current_size, num_pad_bytes;
    char metadata[FCLAW2D_FILE_METADATA_BYTES_V1 + 1];
    char block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + 1];
    p4est_gloidx_t global_num_quadrants;
    fclaw2d_file_section_metadata_v1_t *current_member;
    sc_MPI_Offset current_position;
    sc_MPI_File file;

    FCLAW_ASSERT (p4est != NULL);
    FCLAW_ASSERT (p4est_is_valid (p4est));
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (data_sections != NULL);
    FCLAW_ASSERT (data_sections->elem_size ==
                  sizeof (fclaw2d_file_section_metadata_v1_t));
    FCLAW_ASSERT (errcode != NULL);

    /* set default output values */
    sc_array_reset (data_sections);

    /* open the file in reading mode */
    *errcode = eclass = sc_MPI_SUCCESS; /* MPI defines MPI_SUCCESS to equal 0. */
    file = sc_MPI_FILE_NULL;

    if ((retval =
         sc_io_open (p4est->mpicomm, filename, SC_IO_READ,
                     sc_MPI_INFO_NULL, &file)) != sc_MPI_SUCCESS)
    {
    }

    if (!FCLAW2D_FILE_IS_SUCCESS_V1 (eclass))
    {
        *errcode = eclass;
        SC_FREE (file);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return -1;
    }

    /* read file metadata on root rank */
    FCLAW_ASSERT (FCLAW2D_FILE_IS_SUCCESS_V1 (eclass));
    if (p4est->mpirank == 0)
    {
        if ((retval = sc_io_read_at (file, 0, metadata,
                                     FCLAW2D_FILE_METADATA_BYTES_V1,
                                     sc_MPI_BYTE, &count)) != sc_MPI_SUCCESS)
        {
            *errcode = eclass;
            /* There is no count error for a non-successful read. */
            count_error = 0;
        }
        else
        {
            count_error = (FCLAW2D_FILE_METADATA_BYTES_V1 != count);
        }
    }
    mpiret = sc_MPI_Bcast (&eclass, 1, sc_MPI_INT, 0, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
    if (!FCLAW2D_FILE_IS_SUCCESS_V1 (eclass))
    {
        fclaw2d_file_error_cleanup_v1 (&file);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return -1;
    }
    mpiret = sc_MPI_Bcast (&count_error, 1, sc_MPI_INT, 0, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
    if (count_error)
    {
        if (p4est->mpirank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          " file_info: read count error for file metadata reading");
        }
        *errcode = FCLAW2D_FILE_ERR_COUNT_V1;
        fclaw2d_file_error_cleanup_v1 (&file);
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return -1;
    }

    /* broadcast file metadata to all ranks and null-terminate it */
    mpiret =
        sc_MPI_Bcast (metadata, FCLAW2D_FILE_METADATA_BYTES_V1, sc_MPI_BYTE,
                      0, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
    metadata[FCLAW2D_FILE_METADATA_BYTES_V1] = '\0';

    if ((mpiret =
         fclaw2d_file_check_file_metadata_v1 (p4est->mpicomm, filename,
                                              user_string, metadata,
                                              &global_num_quadrants)) !=
        sc_MPI_SUCCESS)
    {
        *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return fclaw2d_file_error_cleanup_v1 (&file);
    }

    /* check global number of quadrants */
    if (p4est->global_num_quadrants != global_num_quadrants)
    {
        if (p4est->mpirank == 0)
        {
            fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                          " file_info: global number of quadrant mismatch");
        }
        *errcode = FCLAW2D_FILE_ERR_FORMAT_V1;
        fclaw2d_file_error_code_v1 (*errcode, errcode);
        return fclaw2d_file_error_cleanup_v1 (&file);
    }

    current_position =
        (sc_MPI_Offset) (FCLAW2D_FILE_METADATA_BYTES_V1 +
                         FCLAW2D_FILE_BYTE_DIV_V1);

    /* read all data headers that we find and skip the data itself */
    if (p4est->mpirank == 0)
    {
        for (;;)
        {
            /* read block metadata for current record */
            mpiret = sc_io_read_at (file, current_position, block_metadata,
                                    FCLAW2D_FILE_FIELD_HEADER_BYTES_V1,
                                    sc_MPI_BYTE, &count);
            *errcode = eclass;
            if (!FCLAW2D_FILE_IS_SUCCESS_V1 (eclass))
            {
                fclaw2d_file_error_code_v1 (*errcode, errcode);
                return fclaw2d_file_error_cleanup_v1 (&file);
            }
            if (FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 != count)
            {
                /* we did not read the correct number of bytes */
                break;
            }

            /* parse and store the element size, the block type and the user string */
            current_member =
                (fclaw2d_file_section_metadata_v1_t *)
                sc_array_push (data_sections);
            if (block_metadata[0] == 'B' || block_metadata[0] == 'F')
            {
                /* we want to read the block type */
                current_member->block_type = block_metadata[0];
            }
            else
            {
                /* the last entry is incomplete and is therefore removed */
                sc_array_rewind (data_sections,
                                 data_sections->elem_count - 1);
                /* current_member is freed if the whole array is freed */
                break;
            }

            /* check format */
            if (block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + 1] !=
                '\n')
            {
                /* the last entry is incomplete and is therefore removed */
                sc_array_rewind (data_sections,
                                 data_sections->elem_count - 1);
                break;
            }

            /* process block metadata */
            block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 + 1] = '\0';
            /* we cut off the block type specifier to read the data size */
            current_member->data_size = (size_t) sc_atol (&block_metadata[2]);

            /* read the user string */
            /* check '\n' to check the format */
            if (block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 - 1] !=
                '\n')
            {
                /* the last entry is incomplete and is therefore removed */
                sc_array_rewind (data_sections,
                                 data_sections->elem_count - 1);
                break;
            }

            /* null-terminate the user string of the current block */
            block_metadata[FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 - 1] = '\0';

            /* copy the user string, '\0' was already set above */
            sc_strcopy (current_member->user_string,
                        FCLAW2D_FILE_USER_STRING_BYTES_V1,
                        &block_metadata[FCLAW2D_FILE_ARRAY_METADATA_BYTES_V1 +
                                        2]);
            FCLAW_ASSERT (current_member->user_string
                          [FCLAW2D_FILE_USER_STRING_BYTES_V1 - 1] == '\0');

            /* get padding bytes of the current block */
            if (current_member->block_type == 'F')
            {
                current_size =
                    (size_t) (p4est->global_num_quadrants *
                              current_member->data_size);
            }
            else if (current_member->block_type == 'B')
            {
                current_size = current_member->data_size;
            }
            else
            {
                /* \ref fclaw2d_file_read_block_metadata_v1 checks for valid block type */
                SC_ABORT_NOT_REACHED ();
            }
            fcaw2d_file_get_padding_string_v1 (current_size,
                                               FCLAW2D_FILE_BYTE_DIV_V1, NULL,
                                               &num_pad_bytes);
            /* read padding bytes */
            mpiret = sc_io_read_at (file,
                                    current_position +
                                    FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 +
                                    current_size, block_metadata,
                                    num_pad_bytes, sc_MPI_BYTE, &count);
            *errcode = eclass;
            if (!FCLAW2D_FILE_IS_SUCCESS_V1 (eclass))
            {
                return fclaw2d_file_error_cleanup_v1 (&file);
            }
            /* check '\n' in padding bytes */
            if (block_metadata[0] != '\n'
                || block_metadata[num_pad_bytes - 1] != '\n')
            {
                /* the last entry is incomplete and is therefore removed */
                fclaw_errorf ("%s", FCLAW2D_FILE_STRING_V1
                              " file_info: stop parsing file and discard last element "
                              "due to wrong padding format.\n");
                sc_array_rewind (data_sections,
                                 data_sections->elem_count - 1);
                break;
            }
            current_position +=
                FCLAW2D_FILE_FIELD_HEADER_BYTES_V1 + current_size +
                num_pad_bytes;
        }
    }

    /* replicate block metadata in parallel */
    long_header = (long) data_sections->elem_count;     /* 0 on non-root */
    mpiret = sc_MPI_Bcast (&long_header, 1, sc_MPI_LONG, 0, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
    if (p4est->mpirank != 0)
    {
        sc_array_resize (data_sections, (size_t) long_header);
    }
    mpiret = sc_MPI_Bcast (data_sections->array,
                           data_sections->elem_count *
                           data_sections->elem_size, sc_MPI_BYTE, 0,
                           p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    FCLAW_ASSERT (FCLAW2D_FILE_IS_SUCCESS_V1 (eclass));
    /* close the file with error checking */
    fclaw2d_file_error_cleanup_v1 (&file);

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return 0;
}
#endif

/** Converts a error code (MPI or libsc error) into a fclaw2d_file error code.
 * This function turns MPI error codes into MPI error classes if
 * MPI IO is enabled.
 * If errcode is already a fclaw2d errorcode, it just copied to
 * fclaw2d_errcode.
 * If MPI IO is not enabled, the function processes the errors outside
 * of MPI but passes version 1.1 errors to MPI_Error_class.
 * Furthermore, fclaw2d_file functions can create \ref FCLAW2D_FILE_ERR_COUNT_V1
 * as errcode what is also processed by this function.
 * \param [in]  errcode     An errcode from a fclaw2d_file function.
 * \param [out] fclaw2d_errcode Non-NULL pointer. Filled with matching
 *                          error code on success.
 * \return                  FCLAW2D_FILE_ERR_SUCCESS_V1 on successful conversion.
 *                          Other MPI error code otherwise.
 */
static int
fclaw2d_file_error_code_v1 (int errcode, int *fclaw2d_errcode)
{
    FCLAW_ASSERT (fclaw2d_errcode != NULL);
    /* assertion on range of error code input */
    FCLAW_ASSERT ((sc_MPI_SUCCESS <= errcode && errcode < sc_MPI_ERR_LASTCODE)
                  || (FCLAW2D_FILE_ERR_SUCCESS_V1 <= errcode
                      && errcode < FCLAW2D_FILE_ERR_LASTCODE_V1));

    /* copy fclaw2d_file error codes that are not equal to
     * libsc error codes */
    if (FCLAW2D_FILE_ERR_SUCCESS_V1 <= errcode
        && errcode < FCLAW2D_FILE_ERR_LASTCODE_V1)
    {
        /* errcode is a fclaw2d_file errorcode */
        *fclaw2d_errcode = errcode;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    }

    /* map sc_io error codes to fclaw2d error codes */
    switch (errcode)
    {
        /* translate sc_io error codes that are the same as in fclaw2d */
    case sc_MPI_SUCCESS:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_SUCCESS_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_FILE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_FILE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_NOT_SAME:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_NOT_SAME_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_AMODE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_AMODE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_NO_SUCH_FILE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_NO_SUCH_FILE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_FILE_EXISTS:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_FILE_EXIST_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_BAD_FILE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_BAD_FILE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_ACCESS:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_ACCESS_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_NO_SPACE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_NO_SPACE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_QUOTA:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_QUOTA_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_READ_ONLY:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_READ_ONLY_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_FILE_IN_USE:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_IN_USE_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    case sc_MPI_ERR_UNKNOWN:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_UNKNOWN_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;

        /* map sc_io error codes that are summarized in fclaw2d */
    case sc_MPI_ERR_UNSUPPORTED_DATAREP:
    case sc_MPI_ERR_UNSUPPORTED_OPERATION:
    case sc_MPI_ERR_DUP_DATAREP:
    case sc_MPI_ERR_CONVERSION:
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_IO_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;

    default:
        /* errcode may be MPI version 1.1 error code */
        *fclaw2d_errcode = FCLAW2D_FILE_ERR_UNKNOWN_V1;
        return FCLAW2D_FILE_ERR_SUCCESS_V1;
    }
}

/** Turn fclaw2d_file errcode into a string.
 *
 * \param [in] errclass     An errcode that is output by a
 *                          fclaw2d_file function.
 * \param [in,out] string   At least sc_MPI_MAX_ERROR_STRING bytes.
 * \param [out] resultlen   Length of string on return.
 * \return                  \ref FCLAW2D_FILE_ERR_SUCCESS_V1 on success or
 *                          something else on invalid arguments.
 */
static int
fclaw2d_file_error_string_v1 (int errclass, char *string, int *resultlen)
{
    int retval;
    const char *tstr = NULL;

    FCLAW_ASSERT (resultlen != NULL);
    FCLAW_ASSERT (FCLAW2D_FILE_ERR_SUCCESS_V1 <= errclass
                  && errclass < FCLAW2D_FILE_ERR_LASTCODE_V1);

    if (string == NULL || resultlen == NULL)
    {
        return sc_MPI_ERR_ARG;
    }

    /* handle fclaw2d_file_v1-defined error codes */
    switch (errclass)
    {
    case FCLAW2D_FILE_ERR_SUCCESS_V1:
        tstr = "No file error";
        break;
    case FCLAW2D_FILE_ERR_FORMAT_V1:
        tstr = "Wrong file format";
        break;
    case FCLAW2D_FILE_ERR_SECTION_TYPE_V1:
        tstr = "Valid non-matching section type";
        break;
    case FCLAW2D_FILE_ERR_CONN_V1:
        tstr = "Invalid serialized connectivity data";
        break;
    case FCLAW2D_FILE_ERR_P4EST_V1:
        tstr = "Invalid " P4EST_STRING " data";
        break;
    case FCLAW2D_FILE_ERR_IN_DATA_V1:
        tstr = "Invalid input data";
        break;
    case FCLAW2D_FILE_ERR_COUNT_V1:
        tstr =
            "Read or write count error that is not classified as an other error";
        break;
    case FCLAW2D_FILE_ERR_DIM_V1:
        tstr = "The file has the wrong dimension";
        break;

        /* handle error codes as defined in libsc */
    case FCLAW2D_FILE_ERR_FILE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_FILE, string, resultlen);
    case FCLAW2D_FILE_ERR_NOT_SAME_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_NOT_SAME, string, resultlen);
    case FCLAW2D_FILE_ERR_AMODE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_AMODE, string, resultlen);
    case FCLAW2D_FILE_ERR_NO_SUCH_FILE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_NO_SUCH_FILE, string,
                                    resultlen);
    case FCLAW2D_FILE_ERR_FILE_EXIST_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_FILE_EXISTS, string,
                                    resultlen);
    case FCLAW2D_FILE_ERR_BAD_FILE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_BAD_FILE, string, resultlen);
    case FCLAW2D_FILE_ERR_ACCESS_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_ACCESS, string, resultlen);
    case FCLAW2D_FILE_ERR_NO_SPACE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_NO_SPACE, string, resultlen);
    case FCLAW2D_FILE_ERR_QUOTA_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_QUOTA, string, resultlen);
    case FCLAW2D_FILE_ERR_READ_ONLY_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_READ_ONLY, string, resultlen);
    case FCLAW2D_FILE_ERR_IN_USE_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_FILE_IN_USE, string,
                                    resultlen);
    case FCLAW2D_FILE_ERR_UNKNOWN_V1:
        return sc_MPI_Error_string (sc_MPI_ERR_UNKNOWN, string, resultlen);

    default:
        /* no valid fclaw2d file error code */
        SC_ABORT_NOT_REACHED ();
    }
    FCLAW_ASSERT (tstr != NULL);

    /* print into the output string */
    if ((retval = snprintf (string, sc_MPI_MAX_ERROR_STRING, "%s", tstr)) < 0)
    {
        /* unless something goes against the current standard of snprintf */
        return sc_MPI_ERR_NO_MEM;
    }
    if (retval >= sc_MPI_MAX_ERROR_STRING)
    {
        retval = sc_MPI_MAX_ERROR_STRING - 1;
    }
    *resultlen = retval;

    /* we have successfully placed a string in the output variables */
    return sc_MPI_SUCCESS;
}

/** Write a p4est to an opened file.
 * This function does not write the connectvity of the p4est.
 * If one want to write the connectivity, \ref fclaw2d_file_write_connectivity_v1
 * can be used.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            fclaw2d_file_open_create_v1.  It keeps track
 *                            of the data sets written one after another.
 * \param [in]    p4est       The p4est that is written to the file.
 * \param [in]    quad_string The string that is used as user string
 *                            for quadrant section.
 *                            An array of maximal \ref
 *                            FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes that
 *                            is written without the NUL-termination
 *                            after the section-dependent metadata and before
 *                            the actual data. If the char array is shorter the
 *                            written char array will be padded to the
 *                            right by spaces. The user_string is
 *                            written on rank 0 and therefore also only
 *                            required on rank 0. Can be NULL for other
 *                            ranks.
 * \param [in]    quad_data_string  The string that is used as user string
 *                            for quadrant data section if \a save_data is true
 *                            and \a p4est->data_size > 0.
 *                            An array of maximal \ref
 *                            FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes that
 *                            is written without the NUL-termination
 *                            after the section-dependent metadata and before
 *                            the actual data. If the char array is shorter the
 *                            written char array will be padded to the
 *                            right by spaces. The user_string is
 *                            written on rank 0 and therefore also only
 *                            required on rank 0. Can be NULL for other
 *                            ranks.
 * \param [in]    save_data   A Boolean to determine if the quadrant data is
 *                            stored in the file. The quadrant is only stored
 *                            \a p4est->data_size > 0.
 * \param [out]   errcode     An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_write_p4est_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                             p4est_t * p4est, const char *quad_string,
                             const char *quad_data_string, int save_data,
                             int *errcode)
{
    p4est_gloidx_t *pertree;
    sc_array_t arr;
    sc_array_t *quads, *quad_data;
    sc_array_t reshape;
    char write_data_bool[2];
    int write_data;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (p4est != NULL);

    /* initialize */
    pertree = NULL;

    /* allocate memory for pertree */
    pertree =
        FCLAW_ALLOC (p4est_gloidx_t, p4est->connectivity->num_trees + 1);

    /* get count per tree */
    p4est_comm_count_pertree (p4est, pertree);

    sc_array_init_data (&arr, pertree,
                        sizeof (p4est_gloidx_t) *
                        (p4est->connectivity->num_trees + 1), 1);
    /* write count per tree to the file */
    fc = fclaw2d_file_write_block_v1 (fc, arr.elem_size, &arr,
                                      P4EST_STRING " count per tree",
                                      errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc != NULL);
        FCLAW_FREE (pertree);
        return NULL;
    }
    FCLAW_FREE (pertree);

    /* check if we actually write data */
    write_data = save_data && (p4est->data_size > 0);

    quads = p4est_deflate_quadrants (p4est, write_data ? &quad_data : NULL);

    /* fclaw2d_file_write_field_v1 requires per rank local_num_quadrants many elements
     * and therefore we group the data per local quadrant by type casting.
     */
    sc_array_init_reshape (&reshape, quads,
                           FCLAW2D_FILE_COMPRESSED_QUAD_SIZE_V1,
                           p4est->local_num_quadrants);

  /** Write the current p4est to the file; we do not write the
   * connectivity to disk because the connectivity is assumed to
   * be known.
   */
    fc = fclaw2d_file_write_field_v1 (fc, reshape.elem_size, &reshape,
                                      quad_string, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        /* first write call failed */
        sc_array_destroy (quads);
        if (write_data) {
            sc_array_destroy (quad_data);
        }
        return NULL;
    }

    /* We always write a block that indicates if the quadrant data is written
     * or not.
     */
    snprintf (write_data_bool, 2, "%s", write_data ? "1" : "0");
    /* we do not write the terminating '\0' */
    sc_array_init_data (&arr, write_data_bool, 1, 1);
    fc = fclaw2d_file_write_block_v1 (fc, arr.elem_size, &arr,
                                      "Write quadrant data?", errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        /* first write call failed */
        sc_array_destroy (quads);
        if (write_data) {
            sc_array_destroy (quad_data);
        }
        return NULL;
    }

    if (write_data)
    {
        fc = fclaw2d_file_write_field_v1 (fc, quad_data->elem_size, quad_data,
                                          quad_data_string, errcode);
    }

    sc_array_destroy (quads);
    if (write_data)
    {
        sc_array_destroy (quad_data);
    }

    return fc;
}

/** Convert read data to a p4est.
 *
 * \param [in] mpicomm    MPI communicator of the p4est.
 * \param [in] mpisize    Number of MPI ranks.
 * \param [in] conn       connectivity used for the created
 *                        p4est.
 * \param [in] gfq        Global first quadrant array that
 *                        defines the partition of the
 *                        created p4est.
 * \param [in] pertree    The cumulative count per tree.
 * \param [in] quads      An array of compressed quadrants
 *                        that are used to create the new
 *                        p4est. This means an array
 *                        of (P4EST_DIM + 1) \ref p4est_qcoord_t
 *                        that contains the quadrant coordinates
 *                        succeeded by the quadrant level.
 * \param [in] quad_data  An array of quadrant data. This
 *                        array must have as many elements
 *                        as quadrants in the new p4est.
 * \param [out] errcode   An errcode that can be interpreted by \ref
 *                        fclaw2d_file_error_string_v1.
 * \return                A pointer to a newly allocated
 *                        p4est that consists of the given
 *                        quadrants and uses the given
 *                        connectivity.
 */
static p4est_t *
fclaw2d_file_data_to_p4est (sc_MPI_Comm mpicomm, int mpisize,
                            p4est_connectivity_t * conn,
                            const p4est_gloidx_t * gfq,
                            const p4est_gloidx_t * pertree,
                            sc_array_t * quads, sc_array_t * quad_data,
                            int *errcode)
{
    p4est_t *ptemp;
    sc_array_t quads_reshape;

    /* verify call convention and initialize error return */
    FCLAW_ASSERT (conn != NULL);
    FCLAW_ASSERT (gfq != NULL);
    FCLAW_ASSERT (quads != NULL &&
                  quads->elem_size == FCLAW2D_FILE_COMPRESSED_QUAD_SIZE_V1);
    FCLAW_ASSERT (errcode != NULL);
    *errcode = FCLAW2D_FILE_ERR_P4EST_V1;

    /* convert array interpretation for p4est_inflate_null */
    sc_array_init_reshape (&quads_reshape, quads, sizeof (p4est_qcoord_t),
                           (P4EST_DIM + 1) * quads->elem_count);

    ptemp =
        p4est_inflate_null (mpicomm, conn, gfq, pertree, &quads_reshape,
                            quad_data, NULL);
    if (ptemp != NULL)
    {
        *errcode = FCLAW2D_FILE_ERR_SUCCESS_V1;
    }
    return ptemp;
}

/** Read a p4est to an opened file using the MPI communicator of \a fc.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            fclaw2d_file_open_read_v1 (_ext).  It keeps track
 *                            of the data sets read one after another.
 * \param [in]    conn        A connectivity that is used to create
 *                            the \a p4est.
 * \param [in]    data_size   The data size of the p4est that will
 *                            be created by this function. The data size
 *                            must be zero if no quadrant data was stored.
 * \param [out]   p4est       The p4est that is created from the file.
 * \param [in,out] quad_string The user string of the quadrant section.
*                             At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                            The user string is read on rank 0 and internally
 *                            broadcasted to all ranks.
 * \param [in,out] quad_data_string  The user string of the quadrant data section.
                              At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                            The user string is read on rank 0 and internally
 *                            broadcasted to all ranks. The string stays
 *                            untouched if there is not quadrant data.
 * \param [out]   errcode     An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_p4est_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                            p4est_connectivity_t * conn, size_t data_size,
                            p4est_t ** p4est, char *quad_string,
                            char *quad_data_string, int *errcode)
{
    int mpisize, mpiret;
    int written_data;
    char write_bool[2];
    char user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1];
    p4est_topidx_t jt;
    p4est_gloidx_t jq;
    p4est_gloidx_t *gfq, *pertree;
    sc_array_t arr, quadrants, quad_data, pertree_arr;
    p4est_qcoord_t *comp_quad;

    /* verify call convention */
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (conn != NULL);
    FCLAW_ASSERT (p4est_connectivity_is_valid (conn));

    /* initialize error return context */
    FCLAW_ASSERT (p4est != NULL);
    *p4est = NULL;
    FCLAW_ASSERT (errcode != NULL);
    *errcode = FCLAW2D_FILE_ERR_UNKNOWN_V1;
    gfq = NULL;
    sc_array_init_size (&pertree_arr,
                        (conn->num_trees + 1) * sizeof (p4est_gloidx_t), 1);
    sc_array_init (&quadrants, FCLAW2D_FILE_COMPRESSED_QUAD_SIZE_V1);
    /* dummy initialization */
    sc_array_init (&quad_data, 1);

    /* temporary information */
    mpiret = sc_MPI_Comm_size (fc->mpicomm, &mpisize);
    SC_CHECK_MPI (mpiret);

  /** Read data to construct the underlying p4est.
   * One could also use p4est_{load,save} to read and write the p4est
   * and use the p4est_file functions only for quadrant data that is
   * stored externally.
   */

    /* read the count per tree */
    fc = fclaw2d_file_read_block_v1 (fc, pertree_arr.elem_size,
                                     &pertree_arr, quad_string, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        /* first read call failed */
        /* this error occurs in particular for a wrong tree number */
        FCLAW_ASSERT (fc == NULL);
        goto p4est_read_file_p4est_end;
    }

    pertree = (p4est_gloidx_t *) pertree_arr.array;
    /* check the read pertree array */
    if (pertree[0] != 0)
    {
        *errcode = FCLAW2D_FILE_ERR_P4EST_V1;
        sc_array_reset (&pertree_arr);
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1 " file_read_"
                                    P4EST_STRING, errcode);
    }
    for (jt = 0; jt < conn->num_trees; ++jt)
    {
        if (!(pertree[jt] <= pertree[jt + 1]))
        {
            *errcode = FCLAW2D_FILE_ERR_P4EST_V1;
            sc_array_reset (&pertree_arr);
            FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                        FCLAW2D_FILE_STRING_V1 " file_read_"
                                        P4EST_STRING, errcode);
        }
    }
    if (fc->global_num_quadrants != pertree[conn->num_trees])
    {
        *errcode = FCLAW2D_FILE_ERR_P4EST_V1;
        sc_array_reset (&pertree_arr);
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1 " file_read_"
                                    P4EST_STRING, errcode);
    }

    gfq = FCLAW_ALLOC (p4est_gloidx_t, mpisize + 1);
  /** Compute a uniform global first quadrant array to use a uniform
   * partition to read the data fields in parallel.
   */
    p4est_comm_global_first_quadrant (fc->global_num_quadrants, mpisize, gfq);

    FCLAW_ASSERT (gfq[mpisize] == pertree[conn->num_trees]);

    /* read the quadrants */
    fc = fclaw2d_file_read_field_ext_v1 (fc, gfq, quadrants.elem_size,
                                         &quadrants, quad_string, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        /* second read call failed */
        goto p4est_read_file_p4est_end;
    }

    /* check the read quadrants */
    for (jq = 0; jq < (p4est_gloidx_t) quadrants.elem_count; ++jq)
    {
        comp_quad =
            (p4est_qcoord_t *) sc_array_index (&quadrants, (size_t) jq);
        if (!p4est_coordinates_is_valid (comp_quad, comp_quad[P4EST_DIM]))
        {
            *errcode = FCLAW2D_FILE_ERR_P4EST_V1;
            /* clean up local variables and open file context */
            FCLAW_FREE (gfq);
            sc_array_reset (&pertree_arr);
            sc_array_reset (&quadrants);
            FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                        FCLAW2D_FILE_STRING_V1 " file_read_"
                                        P4EST_STRING, errcode);
        }
    }

    /* check if the quadrant data was written to the file */
    sc_array_init_data (&arr, write_bool, 1, 1);
    fc = fclaw2d_file_read_block_v1 (fc, arr.elem_count, &arr, user_string,
                                     errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        /* third read call failed */
        goto p4est_read_file_p4est_end;
    }
    write_bool[1] = '\0';
    written_data = 0;
    if (!strcmp (write_bool, "1"))
    {
        written_data = 1;
    }
    else if (!strcmp (write_bool, "0"))
    {
        written_data = 0;
    }
    else
    {
        /* invalid p4est data */
        *errcode = FCLAW2D_FILE_ERR_P4EST_V1;
        goto p4est_read_file_p4est_end;
    }

    if (written_data)
    {
        sc_array_init (&quad_data, data_size);
        /* read the quadrant data */
        fc = fclaw2d_file_read_field_ext_v1 (fc, gfq, quad_data.elem_size,
                                             &quad_data, quad_data_string,
                                             errcode);
        if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
        {
            FCLAW_ASSERT (fc == NULL);
            /* potential fourth read call failed */
            goto p4est_read_file_p4est_end;
        }
    }
    else {
        FCLAW_ASSERT (data_size == 0);
    }

    /* create the p4est from the read data */
    *p4est =
        fclaw2d_file_data_to_p4est (fc->mpicomm, mpisize, conn, gfq,
                                    (p4est_gloidx_t *) pertree_arr.array,
                                    &quadrants, written_data ? &quad_data : NULL,
                                    errcode);
    FCLAW_ASSERT ((p4est == NULL) ==
                  (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1));
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        /* clean up local variables and open file context */
        FCLAW_FREE (gfq);
        sc_array_reset (&pertree_arr);
        sc_array_reset (&quadrants);
        sc_array_reset (&quad_data);
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1 " file_read_"
                                    P4EST_STRING, errcode);
    }

    /* clean up und return */
  p4est_read_file_p4est_end:
    FCLAW_FREE (gfq);
    sc_array_reset (&pertree_arr);
    sc_array_reset (&quadrants);
    sc_array_reset (&quad_data);
    return fc;
}

/** Write a connectivity to an opened file.
 * This function writes two block sections to the opened file.
 * The first block contains the size of the serialized connectivity data
 * and the second data block contains serialized connectivity.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            fclaw2d_file_open_create_v1.  It keeps track
 *                            of the data sets written one after another.
 * \param [in]    conn        The connectivity that is written to the file.
 * \param [in]    conn_string The user string that written for the
 *                            connectivity data block section.
 *                            An array of maximal \ref
 *                            FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes that
 *                            is written without the NUL-termination
 *                            after the section-dependent metadata and before
 *                            the actual data. If the char array is shorter the
 *                            written char array will be padded to the
 *                            right by spaces. The user_string is
 *                            written on rank 0 and therefore also only
 *                            required on rank 0. Can be NULL for other
 *                            ranks.
 * \param [out]   errcode     An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_write_connectivity_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                                    p4est_connectivity_t * conn,
                                    const char *conn_string, int *errcode)
{
    uint64_t conn_size = 0;
    sc_array_t *conn_buffer, conn_size_arr, reshape;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (conn != NULL);
    FCLAW_ASSERT (conn_string != NULL);
    FCLAW_ASSERT (p4est_connectivity_is_valid (conn));

    /* \ref p4est_connectivity_deflate aborts on errors */
    conn_buffer = p4est_connectivity_deflate (conn, P4EST_CONN_ENCODE_NONE);

    conn_size = (uint64_t) (conn_buffer->elem_size * conn_buffer->elem_count);
    sc_array_init_data (&conn_size_arr, &conn_size, sizeof (uint64_t), 1);
    fc = fclaw2d_file_write_block_v1 (fc, sizeof (size_t),
                                      &conn_size_arr,
                                      P4EST_STRING " connectivity size",
                                      errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        sc_array_destroy (conn_buffer);
        return NULL;
    }

    /* reshape the array to fit for \ref fclaw2d_file_write_block_v1 */
    sc_array_init_reshape (&reshape, conn_buffer,
                           conn_buffer->elem_count * conn_buffer->elem_size,
                           1);
    fc = fclaw2d_file_write_block_v1 (fc, reshape.elem_size, &reshape,
                                      conn_string, errcode);

    /* clean up */
    sc_array_destroy (conn_buffer);

    return fc;
}

/** Read a connectivity from an opened file.
 * This function reads two block sections from the opened file.
 * The first block contains the size of the serialized connectivity data
 * and the second data block contains serialized connectivity.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            fclaw2d_file_open_read_v1 (_ext).  It keeps track
 *                            of the data sets written one after another.
 * \param [out]   conn        The connectivity that is read from the file.
 * \param [in,out] conn_string The user string that read for the
 *                            connectivity data block section.
 *                            At least \ref FCLAW2D_FILE_USER_STRING_BYTES_V1 bytes.
 *                            The user string is read on rank 0 and internally
 *                            broadcasted to all ranks.
 * \param [out]   errcode     An errcode that can be interpreted by \ref
 *                            fclaw2d_file_error_string_v1.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
static fclaw2d_file_context_p4est_v1_t *
fclaw2d_file_read_connectivity_v1 (fclaw2d_file_context_p4est_v1_t * fc,
                                   p4est_connectivity_t ** conn,
                                   char *conn_string, int *errcode)
{
    uint64_t read_conn_size;
    size_t conn_size;
    sc_array_t conn_size_arr;
    sc_array_t conn_arr;
    sc_array_t reshape;

    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (conn != NULL);
    FCLAW_ASSERT (conn_string != NULL);

    sc_array_init_data (&conn_size_arr, &read_conn_size, sizeof (uint64_t),
                        1);
    /* get the connectivity size */
    fc = fclaw2d_file_read_block_v1 (fc, sizeof (uint64_t), &conn_size_arr,
                                     conn_string, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        return NULL;
    }

    conn_size = (size_t) read_conn_size;

    sc_array_init_size (&conn_arr, conn_size, 1);
    /* read the connectivity */
    fc = fclaw2d_file_read_block_v1 (fc, conn_size, &conn_arr, conn_string,
                                     errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS_V1)
    {
        FCLAW_ASSERT (fc == NULL);
        return NULL;
    }

    /* reshape the connectivity data for \ref p4est_connectivity_inflate */
    sc_array_init_reshape (&reshape, &conn_arr, conn_arr.elem_size,
                           sizeof (char));

    /* create the connectivity from the read data */
    *conn = p4est_connectivity_inflate (&reshape);

    /* \ref p4est_connectivity_inflate returns NULL for an invalid
     * connectivity. Therefore, we do not explicitly check for the
     * validity of the returned connectivity.
     */
    if (*conn == NULL)
    {
        /* \ref p4est_connectivity_inflate failed due to wrong format */
        /* close, dealloc file and set specific error code */
        *errcode = FCLAW2D_FILE_ERR_CONN_V1;
        FCLAW2D_FILE_CHECK_NULL_V1 (*errcode, fc,
                                    FCLAW2D_FILE_STRING_V1
                                    " file_read_connectivity", errcode);
    }

    /* clean up */
    sc_array_reset (&conn_arr);

    return fc;
}

/** Close a file opened for parallel write/read and free the context.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in,out] fc       Context previously created by \ref
 *                          fclaw2d_file_open_create_v1 or \ref
 *                          fclaw2d_file_open_read_v1 (_ext).  Is freed.
 * \param [out] errcode     An errcode that can be interpreted by \ref
 *                          fclaw2d_file_error_string_v1.
 * \return                  0 for a successful call and -1 in case of
 *                          an error. See also errcode argument.
 */
static int
fclaw2d_file_close_v1 (fclaw2d_file_context_p4est_v1_t * fc, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int mpiret;

    mpiret = sc_io_close (&fc->file);
    FCLAW2D_FILE_CHECK_INT_V1 (mpiret, "Close file", errcode);

    if (fc->gfq_owned)
    {
        FCLAW_FREE (fc->global_first_quadrant);
    }
    FCLAW_FREE (fc);

    fclaw2d_file_error_code_v1 (*errcode, errcode);
    return 0;
}

/* End of legacy file functions */

typedef struct fclaw2d_file_context
{
    fclaw2d_domain_t *domain;
    fclaw2d_file_context_p4est_v1_t *fc;
    char basename[FCLAW2D_FILE_NAME_BYTES];
}
fclaw2d_file_context_t;

/**
 * This function translates the fclaw2d_file_v1 error codes to fclaw_file
 * error codes.
 *
 * \param [in]  errcode_v1     A fclaw2d_file_v1 error code; see \ref
 *                             fclaw2d_file_error_v1
 * \param [out]  errcode       On output a fclaw2d_file errorcode; see \ref
 *                             fclaw2d_file_error.
 * \return                     0 in case of success and -1 otherwise.
 */
static int
fclaw2d_file_translate_error_code_v1 (int errcode_v1, int *errcode)
{
    FCLAW_ASSERT (errcode != NULL);

    switch (errcode_v1)
    {
    case FCLAW2D_FILE_ERR_SUCCESS_V1:
        *errcode = FCLAW2D_FILE_ERR_SUCCESS;
        return 0;
    case FCLAW2D_FILE_ERR_FILE_V1:
        *errcode = FCLAW2D_FILE_ERR_FILE;
        return 0;
    case FCLAW2D_FILE_ERR_NOT_SAME_V1:
        *errcode = FCLAW2D_FILE_ERR_NOT_SAME;
        return 0;
    case FCLAW2D_FILE_ERR_AMODE_V1:
        *errcode = FCLAW2D_FILE_ERR_AMODE;
        return 0;
    case FCLAW2D_FILE_ERR_NO_SUCH_FILE_V1:
        *errcode = FCLAW2D_FILE_ERR_NO_SUCH_FILE;
        return 0;
    case FCLAW2D_FILE_ERR_FILE_EXIST_V1:
        *errcode = FCLAW2D_FILE_ERR_FILE_EXIST;
        return 0;
    case FCLAW2D_FILE_ERR_BAD_FILE_V1:
        *errcode = FCLAW2D_FILE_ERR_BAD_FILE;
        return 0;
    case FCLAW2D_FILE_ERR_ACCESS_V1:
        *errcode = FCLAW2D_FILE_ERR_ACCESS;
        return 0;
    case FCLAW2D_FILE_ERR_NO_SPACE_V1:
        *errcode = FCLAW2D_FILE_ERR_NO_SPACE;
        return 0;
    case FCLAW2D_FILE_ERR_QUOTA_V1:
        *errcode = FCLAW2D_FILE_ERR_QUOTA;
        return 0;
    case FCLAW2D_FILE_ERR_READ_ONLY_V1:
        *errcode = FCLAW2D_FILE_ERR_READ_ONLY;
        return 0;
    case FCLAW2D_FILE_ERR_IN_USE_V1:
        *errcode = FCLAW2D_FILE_ERR_IN_USE;
        return 0;
    case FCLAW2D_FILE_ERR_IO_V1:
        *errcode = FCLAW2D_FILE_ERR_IO;
        return 0;
    case FCLAW2D_FILE_ERR_FORMAT_V1:
        *errcode = FCLAW2D_FILE_ERR_FORMAT;
        return 0;
    case FCLAW2D_FILE_ERR_SECTION_TYPE_V1:
        *errcode = FCLAW2D_FILE_ERR_SECTION_TYPE;
        return 0;
    case FCLAW2D_FILE_ERR_CONN_V1:
        *errcode = FCLAW2D_FILE_ERR_CONN;
        return 0;
    case FCLAW2D_FILE_ERR_P4EST_V1:
        *errcode = FCLAW2D_FILE_ERR_P4EST;
        return 0;
    case FCLAW2D_FILE_ERR_IN_DATA_V1:
        *errcode = FCLAW2D_FILE_ERR_IN_DATA;
        return 0;
    case FCLAW2D_FILE_ERR_COUNT_V1:
        *errcode = FCLAW2D_FILE_ERR_COUNT;
        return 0;
    case FCLAW2D_FILE_ERR_DIM_V1:
        *errcode = FCLAW2D_FILE_ERR_DIM;
        return 0;
    case FCLAW2D_FILE_ERR_UNKNOWN_V1:
        *errcode = FCLAW2D_FILE_ERR_UNKNOWN;
        return 0;
    case FCLAW2D_FILE_ERR_LASTCODE_V1:
        *errcode = FCLAW2D_FILE_ERR_LASTCODE;
        return 0;
    default:
        *errcode = FCLAW2D_FILE_ERR_UNKNOWN;
        return 0;
    }
    return 0;
}

/**
 * This function translates the fclaw2d_file error codes to fclaw_file_v1
 * error classes.
 *
 * \param [in]  errcode         A fclaw2d_file error code; see \ref
 *                              fclaw2d_file_error
 * \param [out]  errcode_v1     On output a fclaw2d_file_v1 errorcode; see \ref
 *                              fclaw2d_file_error_v1.
 * \return                      0 in case of success and -1 otherwise.
 */
static int
fclaw2d_file_translate_error_code_to_v1 (int errcode, int *errcode_v1)
{
    FCLAW_ASSERT (errcode_v1 != NULL);

    switch (errcode)
    {
    case FCLAW2D_FILE_ERR_SUCCESS:
        *errcode_v1 = FCLAW2D_FILE_ERR_SUCCESS_V1;
        return 0;
    case FCLAW2D_FILE_ERR_FILE:
        *errcode_v1 = FCLAW2D_FILE_ERR_FILE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_NOT_SAME:
        *errcode_v1 = FCLAW2D_FILE_ERR_NOT_SAME_V1;
        return 0;
    case FCLAW2D_FILE_ERR_AMODE:
        *errcode_v1 = FCLAW2D_FILE_ERR_AMODE;
        return 0;
    case FCLAW2D_FILE_ERR_NO_SUCH_FILE:
        *errcode_v1 = FCLAW2D_FILE_ERR_NO_SUCH_FILE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_FILE_EXIST:
        *errcode_v1 = FCLAW2D_FILE_ERR_FILE_EXIST_V1;
        return 0;
    case FCLAW2D_FILE_ERR_BAD_FILE:
        *errcode_v1 = FCLAW2D_FILE_ERR_BAD_FILE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_ACCESS:
        *errcode_v1 = FCLAW2D_FILE_ERR_ACCESS_V1;
        return 0;
    case FCLAW2D_FILE_ERR_NO_SPACE:
        *errcode_v1 = FCLAW2D_FILE_ERR_NO_SPACE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_QUOTA:
        *errcode_v1 = FCLAW2D_FILE_ERR_QUOTA_V1;
        return 0;
    case FCLAW2D_FILE_ERR_READ_ONLY:
        *errcode_v1 = FCLAW2D_FILE_ERR_READ_ONLY_V1;
        return 0;
    case FCLAW2D_FILE_ERR_IN_USE:
        *errcode_v1 = FCLAW2D_FILE_ERR_IN_USE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_IO:
        *errcode_v1 = FCLAW2D_FILE_ERR_IO_V1;
        return 0;
    case FCLAW2D_FILE_ERR_FORMAT:
        *errcode_v1 = FCLAW2D_FILE_ERR_FORMAT_V1;
        return 0;
    case FCLAW2D_FILE_ERR_SECTION_TYPE:
        *errcode_v1 = FCLAW2D_FILE_ERR_SECTION_TYPE_V1;
        return 0;
    case FCLAW2D_FILE_ERR_CONN:
        *errcode_v1 = FCLAW2D_FILE_ERR_CONN_V1;
        return 0;
    case FCLAW2D_FILE_ERR_P4EST:
        *errcode_v1 = FCLAW2D_FILE_ERR_P4EST_V1;
        return 0;
    case FCLAW2D_FILE_ERR_IN_DATA:
        *errcode_v1 = FCLAW2D_FILE_ERR_IN_DATA_V1;
        return 0;
    case FCLAW2D_FILE_ERR_COUNT:
        *errcode_v1 = FCLAW2D_FILE_ERR_COUNT_V1;
        return 0;
    case FCLAW2D_FILE_ERR_DIM:
        *errcode_v1 = FCLAW2D_FILE_ERR_DIM_V1;
        return 0;
    case FCLAW2D_FILE_ERR_UNKNOWN:
        *errcode_v1 = FCLAW2D_FILE_ERR_UNKNOWN_V1;
        return 0;
    case FCLAW2D_FILE_ERR_LASTCODE:
        *errcode_v1 = FCLAW2D_FILE_ERR_LASTCODE_V1;
        return 0;
    default:
        *errcode_v1 = FCLAW2D_FILE_ERR_UNKNOWN_V1;
        return 0;
    }
    return 0;
}

fclaw2d_file_context_t *
fclaw2d_file_open_write (const char *filename,
                         const char *user_string, int write_partition,
                         fclaw2d_domain_t * domain, int *errcode)
{
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (domain->pp != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;
    size_t file_len;
    p4est_wrap_t *wrap;
    p4est_t *p4est;
    fclaw2d_file_context_p4est_v1_t *fc;
    fclaw2d_file_context_t *fclaw_fc;
    char buf[FCLAW2D_FILE_NAME_BYTES];

    file_len = strlen (filename) + strlen ("." FCLAW2D_FILE_EXT) + 1;
    if (file_len > FCLAW2D_FILE_NAME_BYTES)
    {
        /* filename too long */
        *errcode = FCLAW2D_FILE_ERR_BAD_FILE;
        return NULL;
    }

    /* get p4est_wrap_t from domain */
    wrap = (p4est_wrap_t *) domain->pp;
    p4est = wrap->p4est;
    FCLAW_ASSERT (p4est_is_valid (p4est));

    /* create the file */
    sc_strcopy (buf, file_len, filename);
    strcat (buf, "." FCLAW2D_FILE_EXT);
    fc = fclaw2d_file_open_create_v1 (p4est, buf, user_string,
                                      &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        return NULL;
    }

    /* write the connectiviy */
    fc = fclaw2d_file_write_connectivity_v1 (fc, p4est->connectivity,
                                             "p4est connectivity",
                                             &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        return NULL;
    }

    /* write the p4est */
    fc = fclaw2d_file_write_p4est_v1 (fc, p4est, "p4est quadrants",
                                      "p4est quadrant data", 0,
                                      &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        return NULL;
    }

    /* allocate and set flcaw file context */
    fclaw_fc = FCLAW_ALLOC (fclaw2d_file_context_t, 1);
    fclaw_fc->fc = fc;
    fclaw_fc->domain = domain;
    sc_strcopy (fclaw_fc->basename, strlen (filename) + 1, filename);

    return fclaw_fc;
}

fclaw2d_file_context_t *
fclaw2d_file_write_block (fclaw2d_file_context_t *
                          fc, const char *user_string,
                          size_t block_size,
                          sc_array_t * block_data, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (block_data != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;

    fc->fc = fclaw2d_file_write_block_v1 (fc->fc, block_size, block_data, user_string,
                                          &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc->fc == NULL);
        /* The p4est file context was closed and deallocated. */
        /* deallocate fclaw2d file context */
        FCLAW_FREE (fc);
        return NULL;
    }

    return fc;
}

fclaw2d_file_context_t *
fclaw2d_file_write_array (fclaw2d_file_context_t *
                          fc, const char *user_string,
                          size_t patch_size,
                          sc_array_t * patch_data, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (patch_data != NULL);
    FCLAW_ASSERT (patch_data->elem_size == sizeof (sc_array_t));
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;
    char *data_dest, *data_src;
    size_t si;
    sc_array_t contig_arr, *patch_arr;

    /* copy the data to a contigous array  */
    sc_array_init_size (&contig_arr, patch_size, patch_data->elem_count);
    for (si = 0; si < patch_data->elem_count; ++si) {
        patch_arr = (sc_array_t *) sc_array_index (patch_data, si);
        data_src = (char *) sc_array_index (patch_arr, 0);
        data_dest = (char *) sc_array_index (&contig_arr, si);

        memcpy (data_dest, data_src, patch_size);
    }

    /* we use the partiton of the underlying p4est */
    fc->fc = fclaw2d_file_write_field_v1 (fc->fc, patch_size, &contig_arr,
                                          user_string, &errcode_internal);
    sc_array_reset (&contig_arr);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS) {
        FCLAW_ASSERT (fc->fc == NULL);
        /* The p4est file context was closed and deallocated. */
        /* deallocate fclaw2d file context */
        FCLAW_FREE (fc);
        return NULL;
    }

    return fc;
}

fclaw2d_file_context_t *
fclaw2d_file_open_read (const char *filename, char *user_string,
                        sc_MPI_Comm mpicomm, int read_partition,
                        fclaw2d_domain_t ** domain, int *errcode)
{
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;
    size_t file_len;
    char buf[FCLAW2D_FILE_NAME_BYTES];
    char read_user_string[FCLAW2D_FILE_USER_STRING_BYTES_V1 + 1];
    p4est_gloidx_t global_num_quadrants;
    fclaw2d_file_context_p4est_v1_t *p4est_fc;
    fclaw2d_file_context_t *fclaw_fc;
    p4est_connectivity_t *conn;
    p4est_t              *p4est;

    file_len = strlen (filename) + strlen ("." FCLAW2D_FILE_EXT) + 1;
    if (file_len > FCLAW2D_FILE_NAME_BYTES)
    {
        /* filename is too long */
        *errcode = FCLAW2D_FILE_ERR_BAD_FILE;
        return NULL;
    }

    /* WARNING: Currently, we do not handle wrong endianness. */
    /* open the given file */
    sc_strcopy (buf, file_len, filename);
    strcat (buf, "." FCLAW2D_FILE_EXT);
    p4est_fc = fclaw2d_file_open_read_ext_v1 (mpicomm, buf, user_string,
                                              &global_num_quadrants,
                                              &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (p4est_fc == NULL);
        return NULL;
    }

    /* read the p4est connectivity */
    p4est_fc = fclaw2d_file_read_connectivity_v1 (p4est_fc, &conn,
                                                  read_user_string,
                                                  &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (p4est_fc == NULL);
        return NULL;
    }

    /* read the p4est */
    p4est_fc = fclaw2d_file_read_p4est_v1 (p4est_fc, conn, 0, &p4est,
                                           read_user_string, read_user_string,
                                           &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (p4est_fc == NULL);
        return NULL;
    }

    /* create domain from p4est */
    *domain = fclaw2d_domain_new_p4est (p4est);

    /* allocate and set fclaw file context */
    fclaw_fc = FCLAW_ALLOC (fclaw2d_file_context_t, 1);
    fclaw_fc->fc = p4est_fc;
    fclaw_fc->domain = *domain;
    sc_strcopy (fclaw_fc->basename, strlen (filename) + 1, filename);

    return fclaw_fc;
}

fclaw2d_file_context_t *
fclaw2d_file_read_block (fclaw2d_file_context_t *
                         fc, char *user_string,
                         size_t block_size,
                         sc_array_t * block_data, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (block_data != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;

    fc->fc = fclaw2d_file_read_block_v1 (fc->fc, block_size, block_data, user_string,
                                          &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc->fc == NULL);
        /* The p4est file context was closed and deallocated. */
        /* deallocate fclaw2d file context */
        FCLAW_FREE (fc);
        return NULL;
    }

    return fc;
}

fclaw2d_file_context_t *
fclaw2d_file_read_array (fclaw2d_file_context_t *
                         fc, char *user_string,
                         size_t patch_size,
                         sc_array_t * patch_data, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (patch_data != NULL);
    FCLAW_ASSERT (patch_data->elem_size == sizeof (sc_array_t));
    FCLAW_ASSERT (errcode != NULL);

    size_t si;
    int errcode_internal;
    char *data_src, *data_dest;
    sc_array_t contig_arr, *current_arr;

    sc_array_init (&contig_arr, patch_size);
    /* we use the partition of the underlying p4est */
    fc->fc = fclaw2d_file_read_field_v1 (fc->fc, patch_size, &contig_arr,
                                          user_string, &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS) {
        FCLAW_ASSERT (fc->fc == NULL);
        /* The p4est file context was closed and deallocated. */
        /* deallocate fclaw2d file context */
        FCLAW_FREE (fc);
        return NULL;
    }

    /* copy data to discontigous array */
    sc_array_resize (patch_data, contig_arr.elem_count);
    for (si = 0; si < contig_arr.elem_count; ++si) {
        data_src = (char *) sc_array_index (&contig_arr, si);

        current_arr = (sc_array_t *) sc_array_index (patch_data, si);
        sc_array_init_size (current_arr, patch_size, 1);
        data_dest = (char *) sc_array_index (current_arr, 0);

        memcpy (data_dest, data_src, patch_size);
    }
    sc_array_reset (&contig_arr);

    return fc;
}

int
fclaw2d_file_close (fclaw2d_file_context_t * fc, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int retval, errcode_internal;

    retval = fclaw2d_file_close_v1 (fc->fc, &errcode_internal);
    fclaw2d_file_translate_error_code_v1 (errcode_internal, errcode);
    if (*errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        FCLAW_EXECUTE_ASSERT_TRUE (retval != 0);
        return -1;
    }

    FCLAW_FREE (fc);

    return 0;
}

int
fclaw2d_file_error_string (int errcode, char *string, int *resultlen)
{
    FCLAW_ASSERT (string != NULL);
    FCLAW_ASSERT (resultlen != NULL);
    FCLAW_ASSERT (FCLAW2D_FILE_ERR_SUCCESS <= errcode
                  && errcode < FCLAW2D_FILE_ERR_LASTCODE);

    int ret;
    int errcode_v1;
    const char *tstr = NULL;

    /* The version 1 originated error codes come first, FCLAW2D_FILE_ERR_UNKNOWN
     * is the last error code that orgrinates from the first version of the
     * file format.
     */
    if (FCLAW2D_FILE_ERR_SUCCESS <= errcode && errcode <= FCLAW2D_FILE_ERR_UNKNOWN) {
        /* use error string as in file format version 1 implementation */
        fclaw2d_file_translate_error_code_to_v1 (errcode, &errcode_v1);
        return fclaw2d_file_error_string_v1 (errcode_v1, string, resultlen);
    }

    /* handle remaining error codes */
    switch (errcode)
    {
    case FCLAW2D_FILE_ERR_NOT_IMPLEMENTED:
        tstr = "The functionality must be still implemented";
        break;

    default:
        /* no valid fclaw2d file error code */
        SC_ABORT_NOT_REACHED ();
    }

    if ((ret = snprintf (string, sc_MPI_MAX_ERROR_STRING, "%s", tstr)) < 0) {
           return sc_MPI_ERR_NO_MEM;
    }
    if (ret >= sc_MPI_MAX_ERROR_STRING) {
            ret = sc_MPI_MAX_ERROR_STRING - 1;
    }
    *resultlen = ret;

    return sc_MPI_SUCCESS;
}
