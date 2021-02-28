/*
Copyright (c) 2021 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "sgn_options.h"

#include <fclaw2d_clawpatch_options.h>

static int s_sgn_options_package_id = -1;

static void *
sgn_register (sgn_options_t* sgn, sc_options_t * opt)
{
    /* [sgn] User options */

    sc_options_add_double (opt, 0, "g",    &sgn->g, 1.0,  "[sgn] g [1.0]");
    sc_options_add_double (opt, 0, "a",    &sgn->a, 0.1,  "[sgn] a [0.1]");
    sc_options_add_double (opt, 0, "b",    &sgn->b,  12,  "[sgn] b [12]");
    sc_options_add_double (opt, 0, "h0",    &sgn->h0,  12,  "[sgn] h0 [1]");

    sc_options_add_double (opt, 0, "breaking",    &sgn->breaking,  1,  "[sgn] 1 [1]");
    sc_options_add_double (opt, 0, "alpha",    &sgn->alpha,  1.153,  "[sgn] alpha [1.153]");


    sc_options_add_double (opt, 0, "dry-tolerance", &sgn->dry_tolerance, 1e-4,  
                           "[sgn] dry_tolerance [1e-4]");

    sc_options_add_double (opt, 0, "sea-level", &sgn->sea_level, 0,  
                           "[sgn] sea-level [0]");

    sc_options_add_int (opt, 0, "claw-version", &sgn->claw_version, 4,
                           "Clawpack_version (4 or 5) [4]");

    sgn->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
sgn_check (sgn_options_t *sgn)
{
    return FCLAW_NOEXIT;
}

static void
sgn_destroy(sgn_options_t *sgn)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    sgn_options_t *sgn = (sgn_options_t*) package;


    return sgn_register(sgn,opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    sgn_options_t           *sgn;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    sgn = (sgn_options_t*) package;

    return sgn_check(sgn);
}


static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    sgn_options_t *sgn;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    sgn = (sgn_options_t*) package;
    FCLAW_ASSERT (sgn->is_registered);

    sgn_destroy (sgn);

    FCLAW_FREE (sgn);
}


static const fclaw_app_options_vtable_t options_vtable_sgn =
{
    options_register,
    NULL,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

sgn_options_t* sgn_options_register (fclaw_app_t * app,
                                          const char *configfile)
{
    sgn_options_t *sgn;
    FCLAW_ASSERT (app != NULL);

    sgn = FCLAW_ALLOC (sgn_options_t, 1);
    fclaw_app_options_register (app,"sgn", configfile, &options_vtable_sgn,
                                sgn);

    fclaw_app_set_attribute(app,"sgn",sgn);
    return sgn;
}

void sgn_options_store (fclaw2d_global_t* glob, sgn_options_t* sgn)
{
    FCLAW_ASSERT(s_sgn_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,sgn);
    s_sgn_options_package_id = id;
}

sgn_options_t* sgn_get_options(fclaw2d_global_t* glob)
{
    int id = s_sgn_options_package_id;
    return (sgn_options_t*) fclaw_package_get_options(glob, id);    
}

