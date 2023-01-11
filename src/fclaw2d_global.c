/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_filesystem.h>
#include <fclaw_global.h>

#include <fclaw_package.h>
#include <fclaw_timer.h>
#include <fclaw_pointer_map.h>
#include <fclaw_packing.h>

#ifndef P4_TO_P8
#include <fclaw2d_defs.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>

#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_map.h>
#else
#include <fclaw3d_defs.h>
#include <fclaw3d_global.h>

#include <fclaw3d_domain.h>
/* figure out dimension-independent diagnostics */
#include <fclaw3d_map.h>
#endif

void
fclaw2d_iterate_patch_cb
  (fclaw2d_domain_t *domain, fclaw2d_patch_t *patch,
   int blockno, int patchno, void *user)
{
  fclaw_global_iterate_t *gi = (fclaw_global_iterate_t *) user;

  FCLAW_ASSERT (gi->gpcb != NULL);
  gi->gpcb (gi->glob, (fclaw_patch_t *) patch->user, blockno, patchno, gi->user);
}

void
fclaw2d_iterate_family_cb
  (fclaw2d_domain_t *domain, fclaw2d_patch_t *patch,
   int blockno, int patchno, void *user)
{
  fclaw_global_iterate_t *gi = (fclaw_global_iterate_t *) user;
  fclaw_patch_t *family[FCLAW2D_NUMSIBLINGS];
  int i;

  for (i = 0; i < FCLAW2D_NUMSIBLINGS; ++i) {
    family[i] = (fclaw_patch_t *) patch[i].user;
  }

  FCLAW_ASSERT (gi->gfcb != NULL);
  gi->gfcb (gi->glob, family, blockno, patchno, gi->user);
}

/* much of this will eventually move into fclaw_global.c */

fclaw2d_global_t* fclaw2d_global_new (void)
{
    fclaw2d_global_t *glob;

    glob = FCLAW_ALLOC (fclaw2d_global_t, 1);

    /* these variables need to be set after calling this function */
    glob->mpicomm = sc_MPI_COMM_NULL;
    glob->mpisize = 0;
    glob->mpirank = -1;

    glob->pkg_container = fclaw_package_container_new ();
    glob->vtables = fclaw_pointer_map_new ();
    glob->options = fclaw_pointer_map_new ();

    glob->count_amr_advance = 0;
    glob->count_ghost_exchange = 0;
    glob->count_amr_regrid = 0;
    glob->count_amr_new_domain = 0;
    glob->count_multiproc_corner = 0;
    glob->count_grids_per_proc = 0;
    glob->count_grids_remote_boundary = 0;
    glob->count_grids_local_boundary = 0;
    glob->count_single_step = 0;
    glob->count_elliptic_grids = 0;
    glob->curr_time = 0;
    glob->cont = NULL;

#ifndef P4_TO_P8
    /* think about how this can work independent of dimension */
    glob->acc = FCLAW_ALLOC (fclaw2d_diagnostics_accumulator_t, 1);
#endif /* P4_TO_P8 */

    return glob;
}

fclaw2d_global_t* fclaw2d_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank)
{
    fclaw2d_global_t *glob = fclaw2d_global_new ();

    /*
     * Set the communicator.
     * With the current code, overridden by fclaw2d_global_store_domain.
     * Maybe we should streamline this in the future.
     */
    glob->mpicomm = mpicomm;
    glob->mpisize = mpisize;
    glob->mpirank = mpirank;

    return glob;
}

#ifndef P4_TO_P8

static void check_vt(fclaw_userdata_vtable_t* vt, const char* name)
{
    char msg[1024];
    sprintf(msg,"Unregistered vtable for options %s",name);
    SC_CHECK_ABORT ((vt != NULL), msg);
}

static void 
pack_iterator_callback(const char* key, void* value, void* user)
{
    char** buffer_ptr = (char **) user;

    *buffer_ptr += fclaw_pack_string(key, *buffer_ptr);

    fclaw_userdata_vtable_t* vt = fclaw_app_options_get_vtable(key);
    check_vt(vt,key);

    *buffer_ptr += vt->pack(value,*buffer_ptr);
}

size_t 
fclaw2d_global_pack(const fclaw2d_global_t * glob, char* buffer)
{
    const char* buffer_start = buffer;

    buffer += fclaw_pack_double(glob->curr_time, buffer);
    buffer += fclaw_pack_double(glob->curr_dt, buffer);

    buffer += fclaw_pack_size_t(buffer,fclaw_pointer_map_size(glob->options));

    fclaw_pointer_map_iterate(glob->options, pack_iterator_callback, &buffer);

    return (buffer-buffer_start);
}

static void 
packsize_iterator_callback(const char* key, void* value, void* user)
{
    size_t* options_size = (size_t*) user;
    fclaw_userdata_vtable_t* vt = fclaw_app_options_get_vtable(key);
    check_vt(vt,key);
    (*options_size) += fclaw_packsize_string(key) + vt->size(value);
}

size_t 
fclaw2d_global_packsize(const fclaw2d_global_t * glob)
{
    size_t options_size = sizeof(size_t);
    fclaw_pointer_map_iterate(glob->options, packsize_iterator_callback, &options_size);
    return 2*sizeof(double) + options_size;
}

size_t 
fclaw2d_global_unpack(char* buffer, fclaw2d_global_t ** glob_ptr)
{
    char* buffer_start = buffer;

	fclaw2d_global_t* glob = fclaw2d_global_new();
    *glob_ptr = glob;

    buffer += fclaw_unpack_double(buffer,&glob->curr_time);
    buffer += fclaw_unpack_double(buffer,&glob->curr_dt);

    size_t num_option_structs;
    buffer += fclaw_unpack_size_t(buffer,&num_option_structs);
    glob->options = fclaw_pointer_map_new();
    for(size_t i = 0; i< num_option_structs; i++)
    {
        char * key;
        buffer += fclaw_unpack_string(buffer,&key);
        fclaw_userdata_vtable_t* vt = fclaw_app_options_get_vtable(key);
        check_vt(vt,key);
        void * options;
        buffer += vt->unpack(buffer,&options);
        fclaw_pointer_map_insert(glob->options, key, options, vt->destroy);
        FCLAW_FREE(key);
    }

    return buffer-buffer_start;
}

#endif

void
fclaw2d_global_store_domain (fclaw2d_global_t* glob, fclaw2d_domain_t* domain)
{
    glob->domain = domain;

    /* this is redundant if global has been created with a communicator */
    if (glob->mpisize > 0) {
        /* double-check for extra paranoia */
        FCLAW_ASSERT (glob->mpisize == domain->mpisize);
        FCLAW_ASSERT (glob->mpirank == domain->mpirank);
    }
    glob->mpicomm = domain->mpicomm;
    glob->mpisize = domain->mpisize;
    glob->mpirank = domain->mpirank;

    /*
     * This is an assignment that might get removed in the future.
     * There is a separate function, fclow2d_global_store_map (see below),
     * wich accomplishes this without accessing the domain attributes.
     */
    glob->cont = (fclaw2d_map_context_t*)
           fclaw2d_domain_attribute_access (glob->domain, "fclaw_map_context", NULL);
}

void
fclaw2d_global_store_map (fclaw2d_global_t* glob,
                          fclaw2d_map_context_t * map)
{
    glob->cont = map;
}

void
fclaw2d_global_destroy (fclaw2d_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);

    if(glob->pkg_container != NULL) fclaw_package_container_destroy ((fclaw_package_container_t *)glob->pkg_container);
    if(glob->vtables != NULL) fclaw_pointer_map_destroy (glob->vtables);
    if(glob->options != NULL) fclaw_pointer_map_destroy (glob->options);

#ifndef P4_TO_P8
    FCLAW_FREE (glob->acc);
#endif
    FCLAW_FREE (glob);
}

void fclaw2d_global_iterate_level (fclaw2d_global_t * glob, int level,
                                   fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_level (glob->domain, level, pcb, &g);
}

void fclaw2d_global_iterate_patches (fclaw2d_global_t * glob,
                                     fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_patches (glob->domain, pcb, &g);
}

void fclaw2d_global_iterate_families (fclaw2d_global_t * glob,
                                      fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_families (glob->domain, pcb, &g);
}

void fclaw2d_global_iterate_adapted (fclaw2d_global_t * glob, fclaw2d_domain_t* new_domain,
                                     fclaw2d_match_callback_t mcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_adapted (glob->domain, new_domain,mcb,&g);
}

void fclaw2d_global_iterate_level_mthread (fclaw2d_global_t * glob, int level,
                                           fclaw2d_patch_callback_t pcb, void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_level_mthread (glob->domain, level,pcb,&g);
}

void fclaw2d_global_iterate_partitioned (fclaw2d_global_t * glob,
                                         fclaw2d_domain_t * new_domain,
                                         fclaw2d_transfer_callback_t tcb,
                                         void *user)
{
    fclaw2d_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw2d_domain_iterate_partitioned (glob->domain,new_domain,tcb,&g);
}

static fclaw2d_global_t* fclaw2d_global_glob = NULL;

void fclaw2d_global_set_global (fclaw2d_global_t* glob)
{
    FCLAW_ASSERT (fclaw2d_global_glob == NULL);
    fclaw2d_global_glob = glob;
}

void fclaw2d_global_unset_global (void)
{
    FCLAW_ASSERT (fclaw2d_global_glob != NULL);
    fclaw2d_global_glob = NULL;
}

fclaw2d_global_t* fclaw2d_global_get_global (void)
{
    FCLAW_ASSERT(fclaw2d_global_glob != NULL);
    return fclaw2d_global_glob;
}

// Only 2d for now need fclaw2d_options
#ifndef P4_TO_P8

static char* old_path = NULL;

void fclaw2d_set_global_context(fclaw2d_global_t *glob)
{
    fclaw_options_t* opts = fclaw2d_get_options(glob);
    fclaw_set_logging_prefix(opts->logging_prefix);

    // Change run directory
    if(opts->run_directory != NULL){
        FCLAW_ASSERT(old_path == NULL);
        old_path = fclaw_cwd();
        fclaw_cd(opts->run_directory);
    }
}

void fclaw2d_clear_global_context(fclaw2d_global_t *glob)
{
    fclaw_set_logging_prefix(NULL);

    // Return to previous cwd
    if(old_path != NULL){
        fclaw_cd(old_path);
        FCLAW_FREE(old_path);
        old_path = NULL;
    }
}

#endif