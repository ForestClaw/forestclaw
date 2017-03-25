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
#include <fclaw2d_clawpatch3_options.h>
#include <fclaw_package.h>

static int s_clawpatch3_package_id = -1;

static void *
clawpatch3_options_register(fclaw_app_t * a, void *optpkg, sc_options_t * opt)
{
    fclaw2d_clawpatch3_options_t *clawpatch3_options = (fclaw2d_clawpatch3_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (opt != NULL);

    /* allocated storage for this package's option values */
    FCLAW_ASSERT (clawpatch3_options != NULL);

    sc_options_add_int (opt, 0, "mx", &clawpatch3_options->mx, 8,
                        "Number of grid cells per patch in x [8]");

    sc_options_add_int (opt, 0, "my", &clawpatch3_options->my, 8,
                        "Number of grid cells per patch in y [8]");

    sc_options_add_int (opt, 0, "mz", &clawpatch3_options->mz, 8,
                        "Number of grid cells per patch in z [8]");

    sc_options_add_int (opt, 0, "maux", &clawpatch3_options->maux, 0,
                        "Number of auxilliary variables [0]");

    sc_options_add_int (opt, 0, "mbc", &clawpatch3_options->mbc, 2,
                        "Number of ghost cells [2]");

    sc_options_add_int (opt, 0, "meqn", &clawpatch3_options->meqn, 1,
                        "Number of equations [2]");

    /* we do not need to work with the return value */
    clawpatch3_options->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
clawpatch3_options_postprocess(fclaw_app_t * a, void *optpkg, void *registered)
{
    fclaw2d_clawpatch3_options_t *clawpatch3_options = (fclaw2d_clawpatch3_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* postprocess this package */
    FCLAW_ASSERT (clawpatch3_options != NULL);
    FCLAW_ASSERT (clawpatch3_options->is_registered);

    /* Convert strings to arrays (no strings to process here) */

    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
clawpatch3_options_check(fclaw_app_t * app, void *optpkg, void *registered)
{
    fclaw2d_clawpatch3_options_t *clawpatch3_options = (fclaw2d_clawpatch3_options_t *) optpkg;

    if (clawpatch3_options->mx != clawpatch3_options->my)
    {
        fclaw_global_essentialf("Clawpatch error : mx != my\n");
        return FCLAW_EXIT_ERROR;    }

    if (2*clawpatch3_options->mbc > clawpatch3_options->mx)
    {
        fclaw_global_essentialf("Clawpatch error : 2*mbc > mx or 2*mbc > my\n");
        return FCLAW_EXIT_ERROR;
    }

    return FCLAW_NOEXIT;
}

static void
clawpatch3_options_destroy (fclaw_app_t * a, void *optpkg, void *registered)
{
    fclaw2d_clawpatch3_options_t *clawpatch3_options = (fclaw2d_clawpatch3_options_t *) optpkg;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* free this optpkg */
    FCLAW_ASSERT (clawpatch3_options != NULL);
    FCLAW_ASSERT (clawpatch3_options->is_registered);

    FCLAW_FREE (clawpatch3_options);
}

static
const fclaw_app_options_vtable_t fclaw2d_clawpatch3_options_vtable = {
    clawpatch3_options_register,
    clawpatch3_options_postprocess,
    clawpatch3_options_check,
    clawpatch3_options_destroy
};

fclaw2d_clawpatch3_options_t *
fclaw2d_clawpatch3_options_register(fclaw_app_t* app, const char* configfile)
{
    fclaw2d_clawpatch3_options_t* clawpatch3_options;

    FCLAW_ASSERT (app != NULL);

    /* allocate storage for fclaw_options */
    /* we will free it in the options_destroy callback */
    clawpatch3_options = FCLAW_ALLOC(fclaw2d_clawpatch3_options_t,1);

    /* Could also pass in a section header (set to NULL for now) */
    fclaw_app_options_register (app,"clawpatch",
                                configfile,
                                &fclaw2d_clawpatch3_options_vtable,
                                clawpatch3_options);
    fclaw_app_set_attribute(app,"clawpatch",clawpatch3_options);
    return clawpatch3_options;
}

void fclaw2d_clawpatch3_set_options (fclaw2d_global_t *glob, fclaw2d_clawpatch3_options_t* clawpatch3_options)
{
    int id;

    FCLAW_ASSERT(s_clawpatch3_package_id == -1);
    id = fclaw_package_container_add_pkg(glob,
                                         clawpatch3_options);
    s_clawpatch3_package_id = id;
}

static
int fclaw2d_clawpatch3_get_package_id (void)
{
    return s_clawpatch3_package_id;
}

fclaw2d_clawpatch3_options_t* fclaw2d_clawpatch3_get_options(fclaw2d_global_t* glob)
{
    int id = fclaw2d_clawpatch3_get_package_id();
    return (fclaw2d_clawpatch3_options_t*) 
            fclaw_package_get_options(glob, id);
}
