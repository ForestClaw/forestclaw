/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#include "fc3d_clawpack46_options.h"

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw_options.h>
#include <fclaw_package.h>

static int s_clawpack_options_package_id = -1;

static void*
clawpack_register (fc3d_clawpack46_options_t* clawopt, sc_options_t * opt)
{
    fclaw_options_add_int_array (opt, 0, "order", &clawopt->order_string,
                               "2 2", &clawopt->order, 2,
                               "[clawpack] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &clawopt->mcapa, -1,
                        "[clawpack] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &clawopt->src_term, 0,
                         "[clawpack] Source term option [F]");

    sc_options_add_bool (opt, 0, "use-fwaves", &clawopt->use_fwaves, 0,
                         "[clawpack] Use fwave flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &clawopt->mwaves, 1,
                        "[clawpack] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &clawopt->mthlim_string, NULL,
                                 &clawopt->mthlim, clawopt->mwaves,
                                 "[clawpack] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");
    
    /* Array of NumFaces=4 values */
    fclaw_options_add_int_array (opt, 0, "mthbc", &clawopt->mthbc_string, "1 1 1 1",
                                 &clawopt->mthbc, 4,
                                 "[clawpack] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "ascii-out", &clawopt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &clawopt->vtk_out, 0,
                           "Output VTK formatted data [F]");


    clawopt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
clawpack_postprocess (fc3d_clawpack46_options_t * clawopt)
{
    fclaw_options_convert_int_array (clawopt->mthbc_string, &clawopt->mthbc,4);
    
    fclaw_options_convert_int_array (clawopt->mthlim_string, &clawopt->mthlim,
                                     clawopt->mwaves);
    fclaw_options_convert_int_array (clawopt->order_string, &clawopt->order,2);

    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
clawpack_check(fc3d_clawpack46_options_t *clawopt,
                 fclaw2d_clawpatch_options_t *clawpatch_opt)
{
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;
    clawopt->method[6] = clawpatch_opt->maux;

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack : bad maux/mcapa combination\n");
        return FCLAW_EXIT_ERROR;
    }

    /* Should also check mthbc, mthlim, etc. */
    return FCLAW_NOEXIT;
}

static
void clawpack_destroy (fc3d_clawpack46_options_t * clawopt)
{
    fclaw_options_destroy_array (clawopt->mthbc);
    fclaw_options_destroy_array (clawopt->order);
    fclaw_options_destroy_array (clawopt->mthlim);
}

/* ------------------------------------------------------
   Generic calls to options handling;  each calls 
   clawpack-specific options call back
   ------------------------------------------------------ */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc3d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawopt = (fc3d_clawpack46_options_t*) package;

    return clawpack_register(clawopt,opt);
}


static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc3d_clawpack46_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    return clawpack_postprocess (clawopt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack46_options_t *clawopt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc3d_clawpack46_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    clawpatch_opt = (fclaw2d_clawpatch_options_t *)
        fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return clawpack_check(clawopt,clawpatch_opt);    
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc3d_clawpack46_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    clawpack_destroy (clawopt);

    FCLAW_FREE (clawopt);
}

static const fclaw_app_options_vtable_t clawpack_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc3d_clawpack46_options_t*  fc3d_clawpack46_options_register (fclaw_app_t * app,
                                                              const char *configfile)
{
    fc3d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);

    clawopt = FCLAW_ALLOC (fc3d_clawpack46_options_t, 1);
    fclaw_app_options_register (app, "clawpack", configfile,
                                &clawpack_options_vtable, clawopt);
    
    fclaw_app_set_attribute(app,"clawpack",clawopt);
    return clawopt;
}

fc3d_clawpack46_options_t* fc3d_clawpack46_get_options(fclaw2d_global_t *glob)
{
    int id = s_clawpack_options_package_id;
    return (fc3d_clawpack46_options_t*) fclaw_package_get_options(glob,id);
}

void fc3d_clawpack46_options_store (fclaw2d_global_t* glob, fc3d_clawpack46_options_t* clawopt)
{
    int id = fclaw_package_container_add_pkg(glob,clawopt);
    s_clawpack_options_package_id = id;
}
