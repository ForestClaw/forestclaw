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

#include "fc2d_cudaclaw5_options.h"

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw_options.h>

static void*
cudaclaw5_register (fc2d_cudaclaw5_options_t* clawopt, sc_options_t * opt)
{
    fclaw_options_add_int_array (opt, 0, "order", &clawopt->order_string,
                               "2 2", &clawopt->order, 2,
                               "[cudaclaw5] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &clawopt->mcapa, -1,
                        "[cudaclaw5] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &clawopt->src_term, 0,
                         "[cudaclaw5] Source term option [F]");

    sc_options_add_bool (opt, 0, "use_fwaves", &clawopt->use_fwaves, 0,
                         "[cudaclaw5] Use fwaves flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &clawopt->mwaves, 1,
                        "[cudaclaw5] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &clawopt->mthlim_string, NULL,
                                 &clawopt->mthlim, clawopt->mwaves,
                                 "[cudaclaw5] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");
    
    fclaw_options_add_int_array (opt, 0, "mthbc", &clawopt->mthbc_string, "1 1 1 1",
                                 &clawopt->mthbc, 4,
                                 "[cudaclaw5] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "ascii-out", &clawopt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &clawopt->vtk_out, 0,
                           "Output VTK formatted data [F]");

    clawopt->is_registered = 1;
    return NULL;
}

fclaw_exit_type_t
cudaclaw5_postprocess (fc2d_cudaclaw5_options_t * clawopt)
{
    fclaw_options_convert_int_array (clawopt->mthlim_string, &clawopt->mthlim,
                                     clawopt->mwaves);
    fclaw_options_convert_int_array (clawopt->order_string, &clawopt->order,
                                     2);
    fclaw_options_convert_int_array (clawopt->mthbc_string, &clawopt->mthbc,4);

    return FCLAW_NOEXIT;
}


fclaw_exit_type_t
cudaclaw5_check (fc2d_cudaclaw5_options_t * clawopt, 
                 fclaw2d_clawpatch_options_t *clawpatch_opt)
{
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of cudaclaw */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;
    clawopt->method[6] = clawpatch_opt->maux;

    if (clawopt->use_fwaves)
    {
        fclaw_global_essentialf("cudaclaw5 : fwaves not yet implemented\n");
        return FCLAW_EXIT_QUIET;
    }

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("cudaclaw5 : mcapa > 0 but maux == 0.\n");
        return FCLAW_EXIT_ERROR;
    }

    /* There is probably a better place to do this */    
    CUDACLAW5_SET_AMR_MODULE(&clawopt->mwaves, &clawopt->mcapa,
                             clawopt->mthlim, clawopt->method);

    /* Should also check mthbc, mthlim, etc. */

    /* check hardwired dimensions in cudaclaw5_flux2 */
    int check = cudaclaw5_check_dims(clawpatch_opt->meqn, 
                                     clawpatch_opt->maux, 
                                     clawopt->mwaves);
    if (!check)
    {        
        fclaw_global_essentialf("cudaclaw5 : Problem with dimensions " \
                                "setting in cudaclaw5_flux2\n");
        return FCLAW_EXIT_ERROR;
    }

    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}

void
cudaclaw5_destroy (fc2d_cudaclaw5_options_t * clawopt)
{
    fclaw_options_destroy_array (clawopt->mthbc);
    fclaw_options_destroy_array (clawopt->order);
    fclaw_options_destroy_array (clawopt->mthlim);
}


/* -------------------------------------------------------------------
   Functions below are part of the options vtable;  no need to change
    these They call functions above
   ------------------------------------------------------------------- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc2d_cudaclaw5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawopt = (fc2d_cudaclaw5_options_t*) package;

    return cudaclaw5_register(clawopt,opt);
}



static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_cudaclaw5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_cudaclaw5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    return cudaclaw5_postprocess (clawopt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_cudaclaw5_options_t *clawopt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_cudaclaw5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    clawpatch_opt = (fclaw2d_clawpatch_options_t *)
        fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return cudaclaw5_check(clawopt,clawpatch_opt);    /* Nothing can go wrong here! */
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_cudaclaw5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    clawopt = (fc2d_cudaclaw5_options_t*) package;
    FCLAW_ASSERT (clawopt->is_registered);

    cudaclaw5_destroy (clawopt);
    FCLAW_FREE (clawopt);
}

static const fclaw_app_options_vtable_t cudaclaw5_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to cudaclaw options
   ---------------------------------------------------------- */
fc2d_cudaclaw5_options_t*  fc2d_cudaclaw5_options_register (fclaw_app_t * app,
                                                            const char *configfile)
{
    fc2d_cudaclaw5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);

    clawopt = FCLAW_ALLOC (fc2d_cudaclaw5_options_t, 1);
    fclaw_app_options_register (app, "cudaclaw5", configfile,
                                &cudaclaw5_options_vtable, clawopt);

    fclaw_app_set_attribute(app,"cudaclaw5",clawopt);
    return clawopt;
}

fc2d_cudaclaw5_options_t* fc2d_cudaclaw5_get_options(fclaw2d_global_t *glob)
{
    return (fc2d_cudaclaw5_options_t*) fclaw2d_global_get_options(glob, "cudaclaw5");
}

void fc2d_cudaclaw5_options_store (fclaw2d_global_t* glob, fc2d_cudaclaw5_options_t* clawopt)
{
    fclaw2d_global_options_store(glob, "cudaclaw5", clawopt);
}
