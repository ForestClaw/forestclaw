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

#include "fc2d_clawpack46_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fc2d_clawpack46_package
{
  int is_registered;
  fc2d_clawpack46_options_t clawopt;
}
fc2d_clawpack46_package_t; 

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc2d_clawpack46_package_t *clawpkg;
    fc2d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc2d_clawpack46_package_t*) package;
    clawopt = &clawpkg->clawopt;

    FCLAW_ASSERT (clawopt != NULL);

    fclaw_options_add_int_array (opt, 0, "order", &clawopt->order_string,
                               "2 2", &clawopt->order, 2,
                               "[clawpack46] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &clawopt->mcapa, -1,
                        "[clawpack46] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &clawopt->src_term, 0,
                         "[clawpack46] Source term option [F]");

    sc_options_add_bool (opt, 0, "use_fwaves", &clawopt->use_fwaves, 0,
                         "[clawpack46] Use fwaves flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &clawopt->mwaves, 1,
                        "[clawpack46] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &clawopt->mthlim_string, NULL,
                                 &clawopt->mthlim, clawopt->mwaves,
                                 "[clawpack46] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");
    
    /* Array of NumFaces many values */
    fclaw_options_add_int_array (opt, 0, "mthbc", &clawopt->mthbc_string, "1 1 1 1",
                                 &clawopt->mthbc, fclaw2d_NumFaces,
                                 "[clawpack46] Physical boundary condition type [1 1 1 1]");

    clawpkg->is_registered = 1;
    return NULL;
}

fclaw_exit_type_t
fc2d_clawpack46_postprocess (fc2d_clawpack46_options_t * clawopt)
{
    fclaw_options_convert_int_array (clawopt->mthbc_string, &clawopt->mthbc,
                                     fclaw2d_NumFaces);    
    fclaw_options_convert_int_array (clawopt->mthlim_string, &clawopt->mthlim,
                                     clawopt->mwaves);
    fclaw_options_convert_int_array (clawopt->order_string, &clawopt->order,
                                     2);

    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack46_package_t *clawpkg;
    fc2d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc2d_clawpack46_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);

    return fc2d_clawpack46_postprocess (clawopt);
}

#if 0
fclaw_exit_type_t
fc2d_clawpack46_check (fc2d_clawpack46_options_t * clawopt)
{
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;
    clawopt->method[6] = clawopt->maux;

    if (clawopt->use_fwaves)
    {
        fclaw_global_essentialf("clawpack46 : fwaves not yet implemented\n");
        return FCLAW_EXIT_QUIET;
    }

    if (clawopt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack46 : bad maux/mcapa combination\n");
        return FCLAW_EXIT_ERROR;
    }

    /* Should also check mthbc, mthlim, etc. */

    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}
#endif

static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack46_package_t *clawpkg;
    fc2d_clawpack46_options_t *clawopt;

    fclaw2d_clawpatch_options_t *clawpatch_opt;
    clawpatch_opt = fclaw_app_get_attribute(app,"clawpatch",NULL);

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc2d_clawpack46_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);
    
    clawopt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    clawopt->method[1] = clawopt->order[0];
    clawopt->method[2] = clawopt->order[1];
    clawopt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    clawopt->method[4] = clawopt->src_term;
    clawopt->method[5] = clawopt->mcapa;
    clawopt->method[6] = clawpatch_opt->maux;

    if (clawopt->use_fwaves)
    {
        fclaw_global_essentialf("clawpack46 : fwaves not yet implemented\n");
        return FCLAW_EXIT_QUIET;
    }

    if (clawpatch_opt->maux == 0 && clawopt->mcapa > 0)
    {
        fclaw_global_essentialf("clawpack46 : bad maux/mcapa combination\n");
        return FCLAW_EXIT_ERROR;
    }

    /* Should also check mthbc, mthlim, etc. */

    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}

void
fc2d_clawpack46_reset (fc2d_clawpack46_options_t * clawopt)
{
    fclaw_options_destroy_array (clawopt->mthbc);
    fclaw_options_destroy_array (clawopt->order);
    fclaw_options_destroy_array (clawopt->mthlim);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_clawpack46_package_t *clawpkg;
    fc2d_clawpack46_options_t *clawopt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    clawpkg = (fc2d_clawpack46_package_t*) package;
    FCLAW_ASSERT (clawpkg->is_registered);

    clawopt = &clawpkg->clawopt;
    FCLAW_ASSERT (clawopt != NULL);

    fc2d_clawpack46_reset (clawopt);
    FCLAW_FREE (clawpkg);
}

static const fclaw_app_options_vtable_t clawpack46_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc2d_clawpack46_options_t*  fc2d_clawpack46_options_register (fclaw_app_t * app,
                                                              const char *configfile)
{
    fc2d_clawpack46_package_t *clawpkg;

    FCLAW_ASSERT (app != NULL);

    clawpkg = FCLAW_ALLOC (fc2d_clawpack46_package_t, 1);
    fclaw_app_options_register (app, "clawpack46", configfile,
                                &clawpack46_options_vtable, clawpkg);
    fclaw_app_set_attribute(app,"clawpack46",clawpkg);
    return &clawpkg->clawopt;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
