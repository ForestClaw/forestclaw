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

#include <fclaw_filesystem.h>

#include <fclaw_package.h>
#include <fclaw_timer.h>
#include <fclaw_pointer_map.h>
#include <fclaw_packing.h>

#include <fclaw2d_defs.h>
#include <fclaw_global.h>
#include <fclaw_options.h>

#include <fclaw_domain.h>
#include <fclaw_diagnostics.h>
#include <fclaw2d_map.h>

/* much of this will eventually move into fclaw_global.c */

fclaw_global_t* fclaw_global_new (void)
{
    fclaw_global_t *glob;

    glob = FCLAW_ALLOC (fclaw_global_t, 1);

    /* these variables need to be set after calling this function */
    glob->mpicomm = sc_MPI_COMM_NULL;
    glob->mpisize = 0;
    glob->mpirank = -1;

    glob->domain = NULL;

    glob->vtables = fclaw_pointer_map_new ();
    glob->options = fclaw_pointer_map_new ();
    glob->attributes = fclaw_pointer_map_new ();

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

    return glob;
}

fclaw_global_t* fclaw_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank)
{
    fclaw_global_t *glob = fclaw_global_new ();

    /*
     * Set the communicator.
     * With the current code, overridden by fclaw_global_store_domain.
     * Maybe we should streamline this in the future.
     */
    glob->mpicomm = mpicomm;
    glob->mpisize = mpisize;
    glob->mpirank = mpirank;

    return glob;
}

static void check_vt(fclaw_packing_vtable_t* vt, const char* name)
{
    char msg[1024];
    sprintf(msg,"Unregistered options packing vtable for \"%s\"",name);
    SC_CHECK_ABORT ((vt != NULL), msg);
}

static void 
pack_iterator_callback(const char* key, void* value, void* user)
{
    char** buffer_ptr = (char **) user;

    *buffer_ptr += fclaw_pack_string(key, *buffer_ptr);

    fclaw_packing_vtable_t* vt = fclaw_app_get_options_packing_vtable(key);
    check_vt(vt,key);

    // advance buffer pointer
    *buffer_ptr += vt->pack(value,*buffer_ptr);
}

size_t 
fclaw_global_pack(const fclaw_global_t * glob, char* buffer)
{
    if(glob->domain->refine_dim == 3)
    {
        fclaw_abortf("fclaw_global_packsize not implemented for 3d\n");
    }
    const char* buffer_start = buffer;

    buffer += fclaw_pack_double(glob->curr_time, buffer);
    buffer += fclaw_pack_double(glob->curr_dt, buffer);

    buffer += fclaw_pack_size_t(fclaw_pointer_map_size(glob->options), buffer);

    fclaw_pointer_map_iterate(glob->options, pack_iterator_callback, &buffer);

    return (buffer-buffer_start);
}

static void 
packsize_iterator_callback(const char* key, void* value, void* user)
{
    size_t* options_size = (size_t*) user;
    fclaw_packing_vtable_t* vt = fclaw_app_get_options_packing_vtable(key);
    check_vt(vt,key);
    (*options_size) += fclaw_packsize_string(key) + vt->size(value);
}

size_t 
fclaw_global_packsize(const fclaw_global_t * glob)
{
    if(glob->domain->refine_dim == 3)
    {
        fclaw_abortf("fclaw_global_packsize not implemented for 3d\n");
    }
    size_t options_size = sizeof(size_t);
    fclaw_pointer_map_iterate(glob->options, packsize_iterator_callback, &options_size);
    return 2*sizeof(double) + options_size;
}

size_t 
fclaw_global_unpack(char* buffer, fclaw_global_t ** glob_ptr)
{
    char* buffer_start = buffer;

	fclaw_global_t* glob = fclaw_global_new();
    *glob_ptr = glob;

    buffer += fclaw_unpack_double(buffer,&glob->curr_time);
    buffer += fclaw_unpack_double(buffer,&glob->curr_dt);

    size_t num_option_structs;
    buffer += fclaw_unpack_size_t(buffer,&num_option_structs);

    for(size_t i = 0; i< num_option_structs; i++)
    {
        char * key;
        buffer += fclaw_unpack_string(buffer,&key);
        fclaw_packing_vtable_t* vt = fclaw_app_get_options_packing_vtable(key);
        check_vt(vt,key);
        void * options;
        buffer += vt->unpack(buffer,&options);
        fclaw_pointer_map_insert(glob->options, key, options, vt->destroy);
        FCLAW_FREE(key);
    }

    return buffer-buffer_start;
}

void
fclaw_global_store_domain (fclaw_global_t* glob, fclaw_domain_t* domain)
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

    /* done for backwards compatibility */
    fclaw2d_map_context_t* map = (fclaw2d_map_context_t*)
           fclaw_domain_attribute_access (glob->domain, "fclaw_map_context", NULL);

    fclaw2d_map_store(glob, map);
}


void
fclaw_global_destroy (fclaw_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);

    fclaw_pointer_map_destroy (glob->vtables);
    fclaw_pointer_map_destroy (glob->options);
    fclaw_pointer_map_destroy (glob->attributes);

    FCLAW_FREE (glob);
}

void fclaw_global_iterate_level (fclaw_global_t * glob, int level,
                                   fclaw_patch_callback_t pcb, void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_level (glob->domain, level, pcb, &g);
}

void fclaw_global_iterate_patches (fclaw_global_t * glob,
                                     fclaw_patch_callback_t pcb, void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_patches (glob->domain, pcb, &g);
}

void fclaw_global_iterate_families (fclaw_global_t * glob,
                                      fclaw_patch_callback_t pcb, void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_families (glob->domain, pcb, &g);
}

void fclaw_global_iterate_adapted (fclaw_global_t * glob, fclaw_domain_t* new_domain,
                                     fclaw_match_callback_t mcb, void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_adapted (glob->domain, new_domain,mcb,&g);
}

void fclaw_global_iterate_level_mthread (fclaw_global_t * glob, int level,
                                           fclaw_patch_callback_t pcb, void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_level_mthread (glob->domain, level,pcb,&g);
}

void fclaw_global_iterate_partitioned (fclaw_global_t * glob,
                                         fclaw_domain_t * new_domain,
                                         fclaw_transfer_callback_t tcb,
                                         void *user)
{
    fclaw_global_iterate_t g;
    g.glob = glob;
    g.user = user;
    fclaw_domain_iterate_partitioned (glob->domain,new_domain,tcb,&g);
}

void fclaw_global_options_store (fclaw_global_t* glob, const char* key, void* options)
{
    
    if(fclaw_pointer_map_get(glob->options,key) != NULL)
    {
        fclaw_abortf("Options %s already stored in glob\n", key);
    }
    fclaw_pointer_map_insert(glob->options, key, options, NULL);
}

void* fclaw_global_get_options (fclaw_global_t* glob, const char* key)
{
    
    void* options = fclaw_pointer_map_get(glob->options, key);
    if(options == NULL)
    {
        fclaw_abortf("Options %s not found in glob\n", key);
    }
    return options;   
}

void fclaw_global_attribute_store (fclaw_global_t* glob, 
                                     const char* key, 
                                     void* options,
                                     fclaw_pointer_map_value_destroy_t destroy)
{
    
    fclaw_pointer_map_insert(glob->attributes, key, options, destroy);
}

void* fclaw_global_get_attribute (fclaw_global_t* glob, const char* key)
{
    
    return fclaw_pointer_map_get(glob->attributes, key);
}

static fclaw_global_t* fclaw_global_glob = NULL;

void fclaw_global_set_static (fclaw_global_t* glob)
{
    fclaw_global_glob = glob;
}

void fclaw_global_clear_static (void)
{
    fclaw_global_glob = NULL;
}

fclaw_global_t* fclaw_global_get_static_global (void)
{
    return fclaw_global_glob;
}

// Only 2d for now need fclaw2d_options
#ifndef P4_TO_P8

static char* old_path = NULL;

void fclaw_set_global_context(fclaw_global_t *glob)
{
    fclaw_options_t* opts = fclaw_get_options(glob);
    fclaw_set_logging_prefix(opts->logging_prefix);

    // Change run directory
    if(opts->run_directory != NULL){
        FCLAW_ASSERT(old_path == NULL);
        old_path = fclaw_cwd();
        fclaw_cd(opts->run_directory);
    }
}

void fclaw_clear_global_context(fclaw_global_t *glob)
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