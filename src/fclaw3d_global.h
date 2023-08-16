/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW3D_GLOBAL_H
#define FCLAW3D_GLOBAL_H

#include <forestclaw3d.h>  /* Needed to declare callbacks (below) */
#include <fclaw3d_map.h>   /* Needed to store the map context */
#include <fclaw3d_convenience.h> /* Needed for the file context */

#include <fclaw_timer.h>   /* Needed to create statically allocated array of timers */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* these are dimension-specific functions */

void fclaw3d_iterate_patch_cb
    (fclaw3d_domain_t * domain, fclaw3d_patch_t * patch,
     int blockno, int patchno, void *user);

void fclaw3d_iterate_family_cb
    (fclaw3d_domain_t * domain, fclaw3d_patch_t * patch,
     int blockno, int patchno, void *user);

/* much of the following will move into fclaw_global.h */

typedef struct fclaw3d_global fclaw3d_global_t;
typedef struct fclaw3d_global_iterate fclaw3d_global_iterate_t;

struct fclaw3d_global
{
    int count_amr_advance;
    int count_ghost_exchange;
    int count_amr_regrid;
    int count_amr_new_domain;
    int count_single_step;
    int count_elliptic_grids;
    int count_multiproc_corner;
    int count_grids_per_proc;
    int count_grids_remote_boundary;
    int count_grids_local_boundary;
    fclaw2d_timer_t timers[FCLAW2D_TIMER_COUNT];

    /* Time at start of each subcycled time step */
    double curr_time;
    double curr_dt;

    sc_MPI_Comm mpicomm;
    int mpisize;              /**< Size of communicator. */
    int mpirank;              /**< Rank of this process in \b mpicomm. */

    /** Solver packages for internal use. */
    struct fclaw_package_container *pkg_container;

    struct fclaw_pointer_map *vtables;    /**< Vtables */
    struct fclaw_pointer_map *options;    /**< options */

    struct fclaw3d_map_context* cont;
    struct fclaw3d_domain *domain;

#if 0
    /* CB: is this a good place for the accumulator?
           Would it be possible to add and retrieve it as an anonymous
           object that does not need to be known to this file? */

    struct fclaw3d_diagnostics_accumulator *acc;
#endif

    void *user;
};

struct fclaw3d_global_iterate
{
    fclaw3d_global_t* glob;
    void* user;
};

/** Allocate a new global structure. */
fclaw3d_global_t* fclaw3d_global_new (void);

fclaw3d_global_t* fclaw3d_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank);

void fclaw3d_global_destroy (fclaw3d_global_t * glob);

void fclaw3d_global_store_domain (fclaw3d_global_t* glob,
                                  struct fclaw3d_domain* domain);

void fclaw3d_global_store_map (fclaw3d_global_t* glob,
                               fclaw3d_map_context_t * map);

/**
 * @brief Pack global structure into buffer
 * 
 * @param glob the global structure
 * @param buffer the buffer to write to
 * @return size_t number of bytes written
 */
size_t fclaw3d_global_pack(const fclaw3d_global_t * glob, char* buffer);

/**
 * @brief Get the number of bytes needed to pack the global structure
 * 
 * @param glob the structure
 * @return size_t the number of bytes needed to store structure
 */
size_t fclaw3d_global_packsize(const fclaw3d_global_t * glob);

/**
 * @brief Unpack global structure from buffer
 * 
 * @param buffer the buffer to read from
 * @param glob newly create global structure
 * @return size_t number of bytes read
 */
size_t fclaw3d_global_unpack(char* buffer, fclaw3d_global_t** glob);

/** Question: Do we want to fix the glob for each writing workflow */

/** Write the options of global to the opened file.
 *
 * This is a collective function.
 * This function also writes the pack size of the global options to the
 * file.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc         Context previously created by \ref
 *                             fclaw3d_file_open_create.  It keeps track
 *                             of the data sets written one after another.
 * \param [in]     user_string A user string that is written to the file.
 *                             Only \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                             bytes without NUL-termination are
 *                             written to the file. If the user gives less
 *                             bytes the user_string in the file header is padded
 *                             by spaces.
 * \param [in]     glob        The global structure that is used to write the
 *                             options to the file. \b glob->domain must
 *                             coincide with the domain that was passed to open
 *                             the file.
 * \param [out]    errcode     An errcode that can be interpreted by
 *                             \ref fclaw3d_file_error_string.
 * \return                     Return a pointer to input context or NULL in case
 *                             of errors that does not abort the program.
 *                             In case of error the file is tried to close
 *                             and \b fc is freed.
 */
fclaw3d_file_context_t * fclaw3d_file_write_global_opt (fclaw3d_file_context_t *fc,
                                                        const char *user_string,
                                                        fclaw3d_global_t *glob,
                                                        int *errcode);

/** Read the options of global to the an opened file.
 *
 * This is a collective function.
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
 *                            and a trailing NUL-terminat
 * \param [out] glob          Allocated global structure with set options
 *                            according to the read options from file.
 * \param [out] errcode       An errcode that can be interpreted by
 *                            \ref fclaw3d_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and \b fc is freed.
 */
fclaw3d_file_context_t * fclaw3d_file_read_global_opt (fclaw3d_file_context_t *fc,
                                                       const char *user_string,
                                                       fclaw3d_global_t *glob,
                                                       int *errcode);

/** Write the global struct to an opened file such that it can be used to restart.
 *
 * This is a collective function.
 * The data is written in so-called file sections according to a prescribed
 * convention for storing data in ForestClaw.
 * TODO: Explicitly specify and document this convention.
 * TODO: Currently we do not store the diagnostcs and accumulator.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in, out] fc          Context previously created by \ref
 *                              fclaw3d_file_open_create.  It keeps track
 *                              of the data sets written one after another.
 * \param [in]       global     The global structure that is written to the file.
 *                              \b global->domain must coincide with the domain
 *                              that was used to open the file.
 * \param [in]       user_string  A user string that is written to the file.
 *                              Only \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                              bytes without NUL-termination are
 *                              written to the file. If the user gives less
 *                              bytes the user_string in the file header is padded
 *                              by spaces.
 * \param [out]       errcode   An errcode that can be interpreted by
 *                              \ref fclaw3d_file_error_string.
 * \return                      Return a pointer to input context or NULL in case
 *                              of errors that does not abort the program.
 *                              In case of error the file is tried to close
 *                              and \b fc is freed.
 */
fclaw3d_file_context_t * fclaw3d_file_write_global (fclaw3d_file_context_t *fc,
                                                    fclaw3d_global_t *global,
                                                    const char *user_string,
                                                    int *errcode);

/** Read a global struct to an opened file such that it can be used to restart.
 *
 * This is a collective function.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 * Without MPI I/O the function may abort on file system dependent
 * errors.
 *
 * \param [in]  fc            Context previously created by \ref
 *                            fclaw3d_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [in]  filename      The path to the file that is opened.
 * \param [out] user_string   At least \ref FCLAW3D_FILE_USER_STRING_BYTES
 *                            bytes. The user string is written
 *                            to the passed array including padding spaces
 *                            and a trailing NUL-termination.
 * \param [in] replace_fn     Qudrant replace callback that is passed to
 *                            \ref fclaw3d_file_read_domain. This callback can be
 *                            NULL.
 * \param [in] attributes     Attributes of the read domain. Currently the
 *                            attributes of the domain are not written to the
 *                            file by \ref fclaw3d_file_read_domain and
 *                            therefore must be passed to this function. Can be
 *                            NULL.
 * \param [in] wrap_user_pointer  A pointer to anonymous user data that is
 *                            passed to \ref fclaw3d_file_read_domain to create
 *                            the p8est_wrap structure that is used as element
 *                            of the read domain.
 * \param [out] global        Newly allocated global that is read from the file.
 * \param [out] errcode       An errcode that can be interpreted by
 *                            \ref fclaw3d_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
fclaw3d_file_context_t * fclaw3d_file_read_global (fclaw3d_file_context_t *fc,
                                                 const char *filename,
                                                 char *user_string,
                                                 p8est_replace_t replace_fn,
                                                 sc_keyvalue_t *attributes,
                                                 void *wrap_user_pointer,
                                                 fclaw3d_global_t **global,
                                                 int *errcode);

void fclaw3d_global_iterate_level (fclaw3d_global_t * glob, int level,
                                   fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_patches (fclaw3d_global_t * glob,
                                     fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_families (fclaw3d_global_t * glob,
                                      fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_adapted (fclaw3d_global_t * glob,
                                     struct fclaw3d_domain* new_domain,
                                     fclaw3d_match_callback_t mcb, void *user);

void fclaw3d_global_iterate_level_mthread (fclaw3d_global_t * glob, int level,
                                           fclaw3d_patch_callback_t pcb, void *user);

void fclaw3d_global_iterate_partitioned (fclaw3d_global_t * glob,
                                         struct fclaw3d_domain * new_domain,
                                         fclaw3d_transfer_callback_t tcb,
                                         void *user);
/**
 * @brief Store an options structure in the glob
 * 
 * @param glob the global context
 * @param key the key to store the options under
 * @param options the options structure
 */
void fclaw3d_global_options_store (fclaw3d_global_t* glob, const char* key, void* options);

/**
 * @brief Get an options structure from the glob
 * 
 * @param glob the global context
 * @param key the key to retrieve the options from
 * @return void* the options
 */
void* fclaw3d_global_get_options (fclaw3d_global_t* glob, const char* key);

/**
 * @brief Store a glob variable in static memory
 *
 * @param glob the glob variable
 */
void fclaw3d_global_set_global (fclaw3d_global_t* glob);

/**
 * @brief Set the static glob variable to NULL
 */
void fclaw3d_global_unset_global (void);

/**
 * @brief Get the static glob variable
 *
 * @return fclaw2d_global_t* the glob variable
 */
fclaw3d_global_t* fclaw3d_global_get_global (void);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_set_global_context(fclaw3d_global_t *glob);

/**
 * @brief
 *
 * @param glob
 */
void fclaw3d_clear_global_context(fclaw3d_global_t *glob);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW3D_GLOBAL_H */
