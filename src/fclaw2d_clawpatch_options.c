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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw_package.h>

static int s_clawpatch_package_id = -1;

static void *
clawpatch_options_register(fclaw_app_t * a, void *optpkg, sc_options_t * opt)
{
    fclaw2d_clawpatch_options_t *clawpatch_options = (fclaw2d_clawpatch_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (opt != NULL);

    /* allocated storage for this package's option values */
    FCLAW_ASSERT (clawpatch_options != NULL);

    sc_options_add_int (opt, 0, "mx", &clawpatch_options->mx, 8,
                        "Number of grid cells per patch in x [8]");

    sc_options_add_int (opt, 0, "my", &clawpatch_options->my, 8,
                        "Number of grid cells per patch in y [8]");

    sc_options_add_int (opt, 0, "maux", &clawpatch_options->maux, 0,
                        "Number of auxilliary variables [0]");

    sc_options_add_int (opt, 0, "mbc", &clawpatch_options->mbc, 2,
                        "Number of ghost cells [2]");

    /* we do not need to work with the return value */
    clawpatch_options->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
clawpatch_options_postprocess(fclaw_app_t * a, void *optpkg, void *registered)
{
    fclaw2d_clawpatch_options_t *clawpatch_options = (fclaw2d_clawpatch_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* postprocess this package */
    FCLAW_ASSERT (clawpatch_options != NULL);
    FCLAW_ASSERT (clawpatch_options->is_registered);

    /* Convert strings to arrays (no strings to process here) */

    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
clawpatch_options_check(fclaw_app_t * app, void *optpkg, void *registered)
{
    fclaw2d_clawpatch_options_t *clawpatch_options = (fclaw2d_clawpatch_options_t *) optpkg;

    if (clawpatch_options->mx != clawpatch_options->my)
    {
        fclaw_global_essentialf("Clawpatch error : mx != my\n");
        return FCLAW_EXIT_ERROR;    }

    if (2*clawpatch_options->mbc > clawpatch_options->mx)
    {
        fclaw_global_essentialf("Clawpatch error : 2*mbc > mx or 2*mbc > my\n");
        return FCLAW_EXIT_ERROR;
    }

    return FCLAW_NOEXIT;
}

static void
clawpatch_options_destroy (fclaw_app_t * a, void *optpkg, void *registered)
{
    fclaw2d_clawpatch_options_t *clawpatch_options = (fclaw2d_clawpatch_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* free this optpkg */
    FCLAW_ASSERT (clawpatch_options != NULL);
    FCLAW_ASSERT (clawpatch_options->is_registered);

    FCLAW_FREE (clawpatch_options);
}

static
const fclaw_app_options_vtable_t fclaw2d_clawpatch_options_vtable = {
    clawpatch_options_register,
    clawpatch_options_postprocess,
    clawpatch_options_check,
    clawpatch_options_destroy
};

// fclaw2d_clawpatch_options_t* fclaw2d_clawpatch_get_options(fclaw_app_t* app)
// {
//     fclaw2d_clawpatch_options_t*  clawpatch_options;
//     clawpatch_options = (fclaw2d_clawpatch_options_t*) fclaw_app_get_attribute(app,"clawpatch",NULL);
//     FCLAW_ASSERT(clawpatch_options != NULL);
//     return clawpatch_options;
// }

static const fclaw_package_vtable_t clawpatch_vtable_notused = {
    NULL, // patch_data_new
    NULL //patch_data_delete
};

fclaw2d_clawpatch_options_t *
fclaw2d_clawpatch_register(fclaw2d_global_t* glob, fclaw_app_t* app, const char* configfile)
{
    int id; 

    fclaw2d_clawpatch_options_t* clawpatch_options;

    FCLAW_ASSERT (app != NULL);

    /* allocate storage for fclaw_options */
    /* we will free it in the options_destroy callback */
    clawpatch_options = FCLAW_ALLOC(fclaw2d_clawpatch_options_t,1);

    /* Could also pass in a section header (set to NULL for now) */
    fclaw_app_options_register (app,"clawpatch",
                                configfile,
                                &fclaw2d_clawpatch_options_vtable,
                                clawpatch_options);

    fclaw_app_set_attribute(app,"clawpatch",clawpatch_options);

    FCLAW_ASSERT(s_clawpatch_package_id == -1);
    id = fclaw_package_container_add_pkg_new(glob,
                                             clawpatch_options,
                                             &clawpatch_vtable_notused);
    s_clawpatch_package_id = id;

    return clawpatch_options;
}

static
int fclaw2d_clawpatch_get_package_id (void)
{
    return s_clawpatch_package_id;
}

fclaw2d_clawpatch_options_t* fclaw2d_clawpatch_get_options(fclaw2d_global_t* glob)
{
    int id = fclaw2d_clawpatch_get_package_id();
    return (fclaw2d_clawpatch_options_t*) 
            fclaw_package_get_options_new(glob, id);
}
