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
#include <p4est_io.h>
#include <p4est_algorithms.h>
#include <p4est_wrap.h>
#else
#include <fclaw3d_file.h>
#include <p8est_io.h>
#include <p8est_algorithms.h>
#include <p8est_wrap.h>
#endif

typedef struct fclaw2d_file_context
{
    fclaw2d_domain_t *domain;
    p4est_file_context_t *fc;
}
fclaw2d_file_context_t;

fclaw2d_file_context_t *
fclaw2d_file_open_write (const char *filename,
                         const char *user_string,
                         fclaw2d_domain_t * domain, int *errcode)
{
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (domain->pp != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int errcode_internal;
    p4est_wrap_t *wrap;
    p4est_t *p4est;
    p4est_file_context_t *fc;
    fclaw2d_file_context_t *fclaw_fc;
    sc_array_t parameters;
    uint64_t parameters_buffer[2];

    /* get p4est_wrap_t from domain */
    wrap = (p4est_wrap_t *) domain->pp;
    p4est = wrap->p4est;
    FCLAW_ASSERT (p4est_is_valid (p4est));

    /* create the file */
    fc = p4est_file_open_create (p4est, filename, user_string,
                                 &errcode_internal);
    if (errcode_internal != P4EST_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        /* TODO: Introduce more FCLAW_FILE error codes */
        *errcode = errcode_internal;
        return NULL;
    }

    /* write the connectiviy */
    fc = p4est_file_write_connectivity (fc, p4est->connectivity,
                                        "p4est connectivity",
                                        &errcode_internal);
    if (errcode_internal != P4EST_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        /* TODO: Introduce more FCLAW_FILE error codes */
        *errcode = errcode_internal;
        return NULL;
    }

    /* write the p4est */
    fc = p4est_file_write_p4est (fc, p4est, "p4est quadrants",
                                 "p4est quadrant data", &errcode_internal);
    if (errcode_internal != P4EST_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        /* TODO: Introduce more FCLAW_FILE error codes */
        *errcode = errcode_internal;
        return NULL;
    }

    /* pack parameters for \ref p4est_wrap_new_p4est */
    parameters_buffer[0] = (uint64_t) wrap->hollow;
    parameters_buffer[1] = (uint64_t) wrap->btype;
    sc_array_init_data (&parameters, (void *) parameters_buffer, 2 * 64, 1);

    /* write paramters for \ref p4est_wrap_new_p4est */
    fc = p4est_file_write_block (fc, 2 * 64, &parameters,
                                 "p4est wrap hollow and btype",
                                 &errcode_internal);
    if (errcode_internal != P4EST_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (fc == NULL);
        /* TODO: Introduce more FCLAW_FILE error codes */
        *errcode = errcode_internal;
        return NULL;
    }

    /* allocate and set flcaw file context */
    fclaw_fc = FCLAW_ALLOC (fclaw2d_file_context_t, 1);
    fclaw_fc->fc = fc;
    fclaw_fc->domain = domain;

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

    /* The functionality must be still implemented. */
    *errcode = FCLAW2D_FILE_ERR_NOT_IMPLEMENTED;

    return NULL;
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
    FCLAW_ASSERT (errcode != NULL);

    /* The functionality must be still implemented. */
    *errcode = FCLAW2D_FILE_ERR_NOT_IMPLEMENTED;

    return NULL;
}

fclaw2d_file_context_t *
fclaw2d_file_open_read (sc_MPI_Comm mpicomm,
                        const char *filename,
                        char *user_string,
                        fclaw2d_domain_t ** domain, int *errcode)
{
    FCLAW_ASSERT (filename != NULL);
    FCLAW_ASSERT (user_string != NULL);
    FCLAW_ASSERT (domain != NULL);
    FCLAW_ASSERT (errcode != NULL);

    /* The functionality must be still implemented. */
    *errcode = FCLAW2D_FILE_ERR_NOT_IMPLEMENTED;

    return NULL;
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

    /* The functionality must be still implemented. */
    *errcode = FCLAW2D_FILE_ERR_NOT_IMPLEMENTED;

    return NULL;
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
    FCLAW_ASSERT (errcode != NULL);

    /* The functionality must be still implemented. */
    *errcode = FCLAW2D_FILE_ERR_NOT_IMPLEMENTED;

    return NULL;
}

int
fclaw2d_file_close (fclaw2d_file_context_t * fc, int *errcode)
{
    FCLAW_ASSERT (fc != NULL);
    FCLAW_ASSERT (errcode != NULL);

    int retval, errcode_internal;

    retval = p4est_file_close (fc->fc, &errcode_internal);
    if (errcode_internal != P4EST_FILE_ERR_SUCCESS)
    {
        FCLAW_ASSERT (retval != 0);
        /* TODO: Introduce more FCLAW_FILE error codes */
        *errcode = errcode_internal;
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

    int ret;

    if (errcode == FCLAW2D_FILE_ERR_NOT_IMPLEMENTED)
    {
        if ((ret = snprintf (string, sc_MPI_MAX_ERROR_STRING, "%s",
                             "The functionality must be still implemented.\n"))
            < 0)
        {
            return sc_MPI_ERR_NO_MEM;
        }
        if (ret >= sc_MPI_MAX_ERROR_STRING)
        {
            ret = sc_MPI_MAX_ERROR_STRING - 1;
        }
        *resultlen = ret;
    }
    else
    {
        /* Other error codes must be still implemented. */
        return -1;
    }

    return 0;
}
