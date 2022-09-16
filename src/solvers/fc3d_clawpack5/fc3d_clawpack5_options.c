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

#include "fc3d_clawpack5_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fc3d_clawpack5_package
{
  int is_registered;
  fc3d_clawpack5_options_t clawopt;
}
fc3d_clawpack5_package_t; 

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc3d_clawpack5_package_t *clawpkg;
    fc3d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc3d_clawpack5_package_t*) package;
    clawopt = &clawpkg->clawopt;

    FCLAW_ASSERT (clawopt != NULL);

    fclaw_options_add_int_array (opt, 0, "order", &clawopt->order_string,
                               "2 2", &clawopt->order, 2,
                               "[clawpack5] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &clawopt->mcapa, -1,
                        "[clawpack5] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &clawopt->src_term, 0,
                         "[clawpack5] Source term option [F]");

    sc_options_add_bool (opt, 0, "use_fwaves", &clawopt->use_fwaves, 0,
                         "[clawpack5] Use fwaves flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &clawopt->mwaves, 1,
                        "[clawpack5] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &clawopt->mthlim_string, NULL,
                                 &clawopt->mthlim, clawopt->mwaves,
                                 "[clawpack5] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");
    
    fclaw_options_add_int_array (opt, 0, "mthbc", &clawopt->mthbc_string, "1 1 1 1",
                                 &clawopt->mthbc, FCLAW2D_NUMFACES,
                                 "[clawpack5] Physical boundary condition type [1 1 1 1]");

    clawpkg->is_registered = 1;
    return NULL;
}

fclaw_exit_type_t
fc3d_clawpack5_postprocess (fc3d_clawpack5_options_t * clawopt)
{
    fclaw_options_convert_int_array (clawopt->mthlim_string, &clawopt->mthlim,
                                     clawopt->mwaves);
    fclaw_options_convert_int_array (clawopt->order_string, &clawopt->order,
                                     2);
    fclaw_options_convert_int_array (clawopt->mthbc_string, &clawopt->mthbc,
                                     FCLAW2D_NUMFACES);
    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack5_package_t *clawpkg;
    fc3d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc3d_clawpack5_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);

    return fc3d_clawpack5_postprocess (clawopt);
}

static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack5_package_t *clawpkg;
    fc3d_clawpack5_options_t *clawopt;

    fclaw_clawpatch3_options_t *clawpatch_opt;
    clawpatch_opt = fclaw_app_get_attribute(app,"clawpatch",NULL);

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc3d_clawpack5_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);
#if 0
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;
    clawopt->method[6] = clawpatch_opt->maux;

    if (clawopt->use_fwaves)
    {
        fclaw_global_essentialf("clawpack5 : fwaves not yet implemented\n");
        return FCLAW_EXIT_QUIET;
    }

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack5 : bad maux/mcapa combination\n");
        return FCLAW_EXIT_ERROR;
    }
    /// Need to do it after convert option to array
    
    SET_AMR_MODULE(&clawopt->mwaves, &clawopt->mcapa,
                   clawopt->mthlim, clawopt->method);
    /* Should also check mthbc, mthlim, etc. */
#endif
    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}

void
fc3d_clawpack5_reset (fc3d_clawpack5_options_t * clawopt)
{
    fclaw_options_destroy_array (clawopt->mthbc);
    fclaw_options_destroy_array (clawopt->order);
    fclaw_options_destroy_array (clawopt->mthlim);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc3d_clawpack5_package_t *clawpkg;
    fc3d_clawpack5_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc3d_clawpack5_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);

    fc3d_clawpack5_reset (clawopt);
    FCLAW_FREE (clawpkg);
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
fc3d_clawpack5_options_t*  fc3d_clawpack5_options_register (fclaw_app_t * app,
                                                            const char *section,
                                                            const char *configfile)
{
    fc3d_clawpack5_package_t *clawpkg;

    FCLAW_ASSERT (app != NULL);

    clawpkg = FCLAW_ALLOC (fc3d_clawpack5_package_t, 1);
    fclaw_app_options_register (app, section, configfile,
                                &clawpack5_options_vtable, clawpkg);
    fclaw_app_set_attribute(app, section, clawpkg);
    return &clawpkg->clawopt;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
