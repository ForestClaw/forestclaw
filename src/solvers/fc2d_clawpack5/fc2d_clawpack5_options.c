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

#include "fc2d_clawpack5_options.h"
#include "fclaw_base.h"

#include <fclaw2d_clawpatch_options.h>
#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_pointer_map.h>
#include <fclaw_packing.h>

static void*
clawpack5_register (fc2d_clawpack5_options_t* clawopt, sc_options_t * opt)
{
    fclaw_options_add_int_array (opt, 0, "order", &clawopt->order_string,
                               "2 2", &clawopt->order, 2,
                               "[clawpack5] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &clawopt->mcapa, -1,
                        "[clawpack5] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &clawopt->src_term, 0,
                         "[clawpack5] Source term option [F]");

    sc_options_add_bool (opt, 0, "use-fwaves", &clawopt->use_fwaves, 0,
                         "[clawpack5] Use fwaves flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &clawopt->mwaves, 1,
                        "[clawpack5] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &clawopt->mthlim_string, NULL,
                                 &clawopt->mthlim, clawopt->mwaves,
                                 "[clawpack5] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");
    
    fclaw_options_add_int_array (opt, 0, "mthbc", &clawopt->mthbc_string, "1 1 1 1",
                                 &clawopt->mthbc, 4,
                                 "[clawpack5] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "ascii-out", &clawopt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &clawopt->vtk_out, 0,
                           "Output VTK formatted data [F]");

    clawopt->is_registered = 1;
    clawopt->is_unpacked = 0;

    return NULL;
}

fclaw_exit_type_t
clawpack5_postprocess (fc2d_clawpack5_options_t * clawopt)
{
    fclaw_options_convert_int_array (clawopt->mthlim_string, &clawopt->mthlim,
                                     clawopt->mwaves);
    fclaw_options_convert_int_array (clawopt->order_string, &clawopt->order,
                                     2);
    fclaw_options_convert_int_array (clawopt->mthbc_string, &clawopt->mthbc,4);
    
    return FCLAW_NOEXIT;
}


fclaw_exit_type_t
clawpack5_check (fc2d_clawpack5_options_t * clawopt)
{
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;

    /* There is probably a better place to do this */    
    CLAWPACK5_SET_AMR_MODULE(&clawopt->mwaves, &clawopt->mcapa,
                   clawopt->mthlim, clawopt->method, 
                   &clawopt->use_fwaves);

    /* Should also check mthbc, mthlim, etc. */

    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}

void
clawpack5_destroy (fc2d_clawpack5_options_t * clawopt)
{
    fclaw_options_destroy_array (clawopt->mthbc);
    fclaw_options_destroy_array (clawopt->order);
    fclaw_options_destroy_array (clawopt->mthlim);

    //free strings if unpacked
    if (clawopt->is_unpacked)
    {
        FCLAW_FREE((void*) clawopt->order_string);
        FCLAW_FREE((void*) clawopt->mthlim_string);
        FCLAW_FREE((void*) clawopt->mthbc_string);
    }

    FCLAW_FREE (clawopt);
}

static void clawpack5_destroy_void(void* user)
{
    fc2d_clawpack5_options_t* clawopt = (fc2d_clawpack5_options_t*) user;
    clawpack5_destroy(clawopt);
}

static size_t 
options_packsize(void* user)
{
    fc2d_clawpack5_options_t* opts = (fc2d_clawpack5_options_t*) user;

    size_t size = sizeof(fc2d_clawpack5_options_t);
    size += fclaw_packsize_string(opts->order_string);
    size += 2*sizeof(int);  /* order */
    size += opts->mwaves*sizeof(int);  /* mthlim */
    size += fclaw_packsize_string(opts->mthlim_string);
    size += 4*sizeof(int);  /* mthbc */
    size += fclaw_packsize_string(opts->mthbc_string);

    return size;
}

static size_t 
options_pack(void* user, char* buffer)
{
    char* buffer_start = buffer;

    fc2d_clawpack5_options_t* opts = (fc2d_clawpack5_options_t*) user;

    //pack entire struct
    buffer += FCLAW_PACK(*opts, buffer);

    //append arrays to buffer
    buffer += fclaw_pack_string(opts->order_string,buffer);
    buffer += fclaw_pack_int(opts->order[0],buffer);
    buffer += fclaw_pack_int(opts->order[1],buffer);
    for(size_t i = 0; i < opts->mwaves; i++)
    {
        buffer += fclaw_pack_int(opts->mthlim[i],buffer);
    }
    buffer += fclaw_pack_string(opts->mthlim_string,buffer);
    for(size_t i = 0; i < 4; i++)
    {
        buffer += fclaw_pack_int(opts->mthbc[i],buffer);
    }
    buffer += fclaw_pack_string(opts->mthbc_string,buffer);

    return buffer-buffer_start;
}

static size_t 
options_unpack(char* buffer, void** user)
{
    char* buffer_start = buffer;

    fc2d_clawpack5_options_t** opts_ptr = (fc2d_clawpack5_options_t**) user;
    *opts_ptr = FCLAW_ALLOC(fc2d_clawpack5_options_t,1);
    fc2d_clawpack5_options_t* opts = *opts_ptr;

    buffer += FCLAW_UNPACK(buffer, opts);

    //unpack arrays
    buffer += fclaw_unpack_string(buffer,(char**) &opts->order_string);
    opts->order = FCLAW_ALLOC(int,2);
    buffer += fclaw_unpack_int(buffer,&opts->order[0]);
    buffer += fclaw_unpack_int(buffer,&opts->order[1]);
    opts->mthlim = FCLAW_ALLOC(int,opts->mwaves);
    for(size_t i = 0; i < opts->mwaves; i++)
    {
        buffer += fclaw_unpack_int(buffer,&opts->mthlim[i]);
    }
    buffer += fclaw_unpack_string(buffer,(char**) &opts->mthlim_string);
    opts->mthbc = FCLAW_ALLOC(int,4);
    for(size_t i = 0; i < 4; i++)
    {
        buffer += fclaw_unpack_int(buffer,&opts->mthbc[i]);
    }
    buffer += fclaw_unpack_string(buffer,(char**) &opts->mthbc_string);

    opts->is_unpacked = 1;
   
    return buffer-buffer_start;
}

static fclaw_packing_vtable_t packing_vt = 
{
	options_pack,
	options_unpack,
	options_packsize,
	clawpack5_destroy_void
};

const fclaw_packing_vtable_t* 
fc2d_clawpack5_options_get_packing_vtable()
{
    return &packing_vt;
}


/* -------------------------------------------------------------------
   Functions below are part of the options vtable;  no need to change
    these They call functions above
   ------------------------------------------------------------------- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fclaw_app_register_options_packing_vtable("fc2d_clawpack5", &packing_vt);

    fc2d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawopt = (fc2d_clawpack5_options_t*) package;

    return clawpack5_register(clawopt,opt);
}



static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_clawpack5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    return clawpack5_postprocess (clawopt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_clawpack5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    return clawpack5_check(clawopt);    /* Nothing can go wrong here! */
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_clawpack5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    clawpack5_destroy (clawopt);
}

static const fclaw_app_options_vtable_t clawpack5_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc2d_clawpack5_options_t*  fc2d_clawpack5_options_register (fclaw_app_t * app,
                                                            const char *section,
                                                            const char *configfile)
{
    fc2d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);

    clawopt = FCLAW_ALLOC (fc2d_clawpack5_options_t, 1);
    fclaw_app_options_register (app, section, configfile,
                                &clawpack5_options_vtable, clawopt);

    fclaw_app_set_attribute(app, section, clawopt);
    return clawopt;
}

fc2d_clawpack5_options_t* fc2d_clawpack5_get_options(fclaw_global_t *glob)
{
    return (fc2d_clawpack5_options_t*) fclaw_global_get_options(glob,"fc2d_clawpack5");
}

void fc2d_clawpack5_options_store (fclaw_global_t* glob, fc2d_clawpack5_options_t* clawopt)
{
    fclaw_global_options_store(glob,"fc2d_clawpack5",clawopt);
}
