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

/** This is a private header that is not intented to be used by the user.
 *
 * This file exports macros and error codes of the version data file format,
 * which was moved from p4est to ForestClaw. The only purpose of this file is to
 * make the respective macros available to fclaw2d_file.h.
 */

#ifndef FCLAW3D_FILE_V1_H
#define FCLAW3D_FILE_V1_H

#include<sc_mpi.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

#define FCLAW3D_FILE_USER_STRING_BYTES_V1 48 /**< number of user string bytes */
#define FCLAW3D_FILE_MAX_GLOBAL_QUAD_V1 9999999999999999 /**< maximal number of global quadrants */
#define FCLAW3D_FILE_MAX_BLOCK_SIZE_V1 9999999999999 /**< maximal number of block bytes */
#define FCLAW3D_FILE_MAX_FIELD_ENTRY_SIZE_V1 9999999999999 /**< maximal number of bytes per field entry */

/** Error values for fclaw3d_file_v1 functions.
 */
typedef enum fclaw3d_file_error_v1
{
    FCLAW3D_FILE_ERR_SUCCESS_V1 = sc_MPI_ERR_LASTCODE,
                                                     /**< file function completed with success */
    FCLAW3D_FILE_ERR_FILE_V1,
                            /**< invalid file handle */
    FCLAW3D_FILE_ERR_NOT_SAME_V1,
                                /**< collective arg not identical */
    FCLAW3D_FILE_ERR_AMODE_V1,
                             /**< access mode error */
    FCLAW3D_FILE_ERR_NO_SUCH_FILE_V1,
                                    /**< file does not exist */
    FCLAW3D_FILE_ERR_FILE_EXIST_V1,
                                  /**< file exists already */
    FCLAW3D_FILE_ERR_BAD_FILE_V1,
                                /**< invalid file name */
    FCLAW3D_FILE_ERR_ACCESS_V1,
                              /**< permission denied */
    FCLAW3D_FILE_ERR_NO_SPACE_V1,
                                /**< not enough space */
    FCLAW3D_FILE_ERR_QUOTA_V1,
                             /**< quota exceeded */
    FCLAW3D_FILE_ERR_READ_ONLY_V1,
                                 /**< read only file (system) */
    FCLAW3D_FILE_ERR_IN_USE_V1,
                              /**< file currently open by other process */
    FCLAW3D_FILE_ERR_IO_V1,
                          /**< other I/O error */
    FCLAW3D_FILE_ERR_FORMAT_V1,/**< read file has a wrong format */
    FCLAW3D_FILE_ERR_SECTION_TYPE_V1,
                                    /**< a valid non-matching section type */
    FCLAW3D_FILE_ERR_CONN_V1,
                            /**< invalid serialized connectivity data */
    FCLAW3D_FILE_ERR_P4EST_V1,
                             /**< invalid p4est data */
    FCLAW3D_FILE_ERR_IN_DATA_V1,
                               /**< input data of file function is invalid */
    FCLAW3D_FILE_ERR_COUNT_V1, /**< read or write count error that was not
                                 classified as a format error */
    FCLAW3D_FILE_ERR_UNKNOWN_V1,
                               /**< unknown error */
    FCLAW3D_FILE_ERR_LASTCODE_V1
                               /**< to define own error codes for
                                  a higher level application
                                  that is using fclaw3d_file_v1
                                  functions */
}
fclaw3d_file_error_v1_t;

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW3D_FILE_H */
