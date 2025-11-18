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
#include <fclaw_map.h>
#include <fclaw_packing.h>

static
fclaw_global_t* global_new (void)
{
    fclaw_global_t *glob;

    glob = FCLAW_ALLOC (fclaw_global_t, 1);

    /* these variables need to be set after calling this function */
    glob->domain = NULL;
    glob->mpicomm = sc_MPI_COMM_NULL;
    glob->mpisize = 0;
    glob->mpirank = -1;

    glob->pkg_container = fclaw_package_container_new ();
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
    glob->curr_dt = 0;
    glob->cont = NULL;
    glob->user = NULL;

    return glob;
}

fclaw_global_t* fclaw_global_new (fclaw_app_t * app)
{
    fclaw_global_t *glob = global_new ();

    glob->mpicomm = fclaw_app_get_mpi_size_rank(app, &glob->mpisize, &glob->mpirank);

    sc_options_t *sc_options = fclaw_app_get_options(app);
    fclaw_global_attribute_store(glob,
                                 "fclaw_options",
                                 sc_options,
                                 NULL, NULL);

    void *fclaw_opt_sections = fclaw_app_get_attribute(app, "fclaw_opt_sections", NULL);
    fclaw_global_attribute_store(glob,
                                 "fclaw_opt_sections",
                                 fclaw_opt_sections,
                                 NULL, NULL);

    return glob;
}

fclaw_global_t* fclaw_global_new_comm (sc_MPI_Comm mpicomm,
                                           int mpisize, int mpirank)
{
    fclaw_global_t *glob = global_new ();

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
}

void
fclaw_map_store (fclaw_global_t* glob,
                          fclaw_map_context_t * map)
{
    glob->cont = map;
}

fclaw_map_context_t*
fclaw_map_get(fclaw_global_t* glob)
{
    return glob->cont;
}

void
fclaw_global_destroy (fclaw_global_t * glob)
{
    FCLAW_ASSERT (glob != NULL);

    if(glob->pkg_container != NULL) fclaw_package_container_destroy ((fclaw_package_container_t *)glob->pkg_container);
    if(glob->vtables != NULL) fclaw_pointer_map_destroy (glob->vtables);
    if(glob->options != NULL) fclaw_pointer_map_destroy (glob->options);
    if(glob->attributes != NULL) fclaw_pointer_map_destroy (glob->attributes);

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
    
    fclaw_pointer_map_insert(glob->options, key, options, NULL);
}

void* fclaw_global_get_options (fclaw_global_t* glob, const char* key)
{
    
    void* options = fclaw_pointer_map_get(glob->options, key);
    FCLAW_ASSERT(options != NULL);
    return options;   
}

/*
 * @brief struct for store attributes and metadata
 */
typedef
struct attribute_entry
{
    /* the key to the packing vtable, NULL if not used */
    char* packing_vtable_key;
    /* the attribute */
    void* attribute;
    /* the callback to destroy the attribute, NULL if not used*/
    fclaw_pointer_map_value_destroy_t destroy;
} attribute_entry_t;

/* callback to destroy attribute entry */
static void 
attribute_entry_destroy(void* value)
{
    attribute_entry_t* entry = (attribute_entry_t*) value;
    if(entry->destroy != NULL)
    {
        entry->destroy(entry->attribute);
    }
    FCLAW_FREE(entry->packing_vtable_key);
    FCLAW_FREE(entry);
}

void 
fclaw_global_attribute_store (fclaw_global_t * glob, 
                              const char * key, 
                              void* attribute,
                              const char * packing_vtable_key, 
                              fclaw_pointer_map_value_destroy_t destroy)
{
    attribute_entry_t *entry = FCLAW_ALLOC(attribute_entry_t,1);
    if(packing_vtable_key != NULL)
    {
        entry->packing_vtable_key = FCLAW_ALLOC(char,strlen(packing_vtable_key)+1);
        strcpy(entry->packing_vtable_key,packing_vtable_key);
    }
    else
    {
        entry->packing_vtable_key = NULL;
    }
    entry->attribute = attribute;
    entry->destroy = destroy;
    fclaw_pointer_map_insert(glob->attributes, key, entry, attribute_entry_destroy);
}

void * 
fclaw_global_get_attribute (fclaw_global_t* glob, const char* key)
{
    attribute_entry_t *entry = 
        (attribute_entry_t*) fclaw_pointer_map_get(glob->attributes, key);
    void * attribute = NULL;
    if(entry != NULL)
    {
        attribute = entry->attribute;
    }
    return attribute;
}

/* ***************************************
 *  Packing and unpacking functions
 * ***************************************/

static void 
check_vt(fclaw_packing_vtable_t* vt, 
         const char* name, 
         const char* vtable_name)
{
    if(vt == NULL)
    {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Unregistered attribute packing vtable \"%s\" for attribute \"%s\"",vtable_name,name);
        SC_CHECK_ABORT ((vt != NULL), msg);
    }
}

static void 
num_to_pack_cb(const char* key, void* value, void* user)
{
    size_t *num_to_pack = (size_t*) user;
    attribute_entry_t* entry = (attribute_entry_t*) value;
    if(entry->packing_vtable_key != NULL)
    {
        (*num_to_pack)++;
    }
}
typedef
struct pack_iter
{
    fclaw_global_t* glob;
    char** buffer_ptr;
    size_t num_packed;
} pack_iter_t;

static void 
pack_attribute_cb(const char* key, void* value, void* user)
{
    pack_iter_t *iter = (pack_iter_t*) user;
    attribute_entry_t* entry = (attribute_entry_t*) value;
    if(entry->packing_vtable_key != NULL)
    {
        fclaw_packing_vtable_t* vt 
            = (fclaw_packing_vtable_t*) fclaw_global_get_vtable(iter->glob, entry->packing_vtable_key);
        check_vt(vt, key, entry->packing_vtable_key);
        // advance buffer pointer
        *iter->buffer_ptr += fclaw_pack_string(entry->packing_vtable_key, *iter->buffer_ptr);
        *iter->buffer_ptr += fclaw_pack_string(key, *iter->buffer_ptr);
        *iter->buffer_ptr += vt->pack(iter->glob, entry->attribute, *iter->buffer_ptr);
    }
}

size_t 
fclaw_global_pack(fclaw_global_t * glob, char* buffer)
{
    const char* buffer_start = buffer;

    buffer += fclaw_pack_double(glob->curr_time, buffer);
    buffer += fclaw_pack_double(glob->curr_dt, buffer);

    size_t num_to_pack  = 0;
    fclaw_pointer_map_iterate(glob->attributes, num_to_pack_cb, &num_to_pack);
    buffer += fclaw_pack_size_t(num_to_pack, buffer);

    pack_iter_t iter;
    iter.glob = glob;
    iter.buffer_ptr = &buffer;
    fclaw_pointer_map_iterate(glob->attributes, pack_attribute_cb, &iter);

    return (buffer-buffer_start);
}

typedef
struct packsize_iter
{
    fclaw_global_t* glob;
    size_t size;
} packsize_iter_t;

static void 
attribute_packsize_cb(const char* key, void* value, void* user)
{
    packsize_iter_t* iter = (packsize_iter_t*) user;
    attribute_entry_t* entry = (attribute_entry_t*) value;
    if(entry->packing_vtable_key != NULL)
    {
        fclaw_packing_vtable_t* vt = 
            (fclaw_packing_vtable_t*) fclaw_global_get_vtable(iter->glob, entry->packing_vtable_key);
        check_vt(vt, key, entry->packing_vtable_key);
        iter->size += 
            fclaw_packsize_string(entry->packing_vtable_key) 
            + fclaw_packsize_string(key) 
            + vt->size(iter->glob, entry->attribute);
    }
}

size_t 
fclaw_global_packsize(fclaw_global_t * glob)
{
    packsize_iter_t iter;
    iter.glob = glob;
    iter.size = sizeof(size_t);
    fclaw_pointer_map_iterate(glob->attributes, attribute_packsize_cb, &iter);
    return 2*sizeof(double) + iter.size;
}

size_t 
fclaw_global_unpack(char * buffer, fclaw_global_t * glob)
{
    char *buffer_start = buffer;

    buffer += fclaw_unpack_double(buffer,&glob->curr_time);
    buffer += fclaw_unpack_double(buffer,&glob->curr_dt);

    size_t num_attributes;
    buffer += fclaw_unpack_size_t(buffer,&num_attributes);

    for(size_t i = 0; i< num_attributes; i++)
    {
        char *packing_vtable_key;
        buffer += fclaw_unpack_string(buffer,&packing_vtable_key);

        char *attribute_key;
        buffer += fclaw_unpack_string(buffer,&attribute_key);

        fclaw_packing_vtable_t *vt = 
            (fclaw_packing_vtable_t *) fclaw_global_get_vtable(glob, packing_vtable_key);

        check_vt(vt, attribute_key, packing_vtable_key);

        attribute_entry_t *entry 
            = (attribute_entry_t *) fclaw_pointer_map_get(glob->attributes, attribute_key);

        if(entry == NULL)
        {
            entry = FCLAW_ALLOC(attribute_entry_t,1);
            entry->packing_vtable_key = packing_vtable_key;
            entry->destroy = vt->destroy;
            fclaw_pointer_map_insert(glob->attributes, attribute_key, entry, attribute_entry_destroy);
        }
        else
        {
            FCLAW_ASSERT(strcmp(entry->packing_vtable_key,packing_vtable_key) == 0);
            FCLAW_FREE(packing_vtable_key);
            entry->destroy(entry->attribute);
        }

        entry->attribute = vt->new_data(glob);

        buffer += vt->unpack(glob, buffer, entry->attribute);

        FCLAW_FREE(attribute_key);
    }

    return buffer-buffer_start;
}
void 
fclaw_global_vtable_store (fclaw_global_t * glob, 
                              const char * key, 
                              void * vtable,
                              fclaw_pointer_map_value_destroy_t destroy)
{
    fclaw_pointer_map_insert(glob->vtables, key, vtable, destroy);
}

void * 
fclaw_global_get_vtable (fclaw_global_t* glob, const char* key)
{
    return fclaw_pointer_map_get(glob->vtables, key);
}

static fclaw_global_t* fclaw2d_global_glob = NULL;

void fclaw_global_set_static (fclaw_global_t* glob)
{
    FCLAW_ASSERT (fclaw2d_global_glob == NULL);
    fclaw2d_global_glob = glob;
}

void fclaw_global_clear_static (void)
{
    FCLAW_ASSERT (fclaw2d_global_glob != NULL);
    fclaw2d_global_glob = NULL;
}

fclaw_global_t* fclaw_global_get_static_global (void)
{
    FCLAW_ASSERT(fclaw2d_global_glob != NULL);
    return fclaw2d_global_glob;
}

// Only 2d for now need fclaw2d_options
#ifndef P4_TO_P8

static char* old_path = NULL;

void fclaw_set_global_context(fclaw_global_t *glob)
{
    fclaw_options_t* opts = fclaw_get_options(glob);

    // Change run directory
    if(strcmp(opts->run_directory,"") != 0){
        FCLAW_ASSERT(old_path == NULL);
        fclaw_set_logging_prefix(opts->logging_prefix);
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