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

#include "disk_user.h"

static void *
disk_register (user_options_t *user_opt, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user_opt->example, 0,
                        "[user] 0 = pillowdisk;  1 = pillowdisk5;  [0]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");

    sc_options_add_double (opt, 0, "alpha", &user_opt->alpha, 0.5,
                           "[user] Ratio used for squared- and pillow-disk [0.5]");

    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
disk_postprocess(user_options_t *user)
{
    /* nothing to post-process */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
disk_check (user_options_t *user_opt)
{

    if (user_opt->example < 0 || user_opt->example > 1) {
        fclaw_global_essentialf ("Option --user:example must be 0 or 1\n");
    }
    /* Print the summary after we have done a global option check */
    return FCLAW_EXIT_QUIET;
}

static void
disk_destroy(user_options_t *user)
{
    /* Nothing to destroy */
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user_opt = (user_options_t*) package;

    return disk_register(user_opt,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    user_options_t *user_opt = (user_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user_opt->is_registered);

    /* Convert strings to arrays */
    return disk_postprocess (user_opt);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user_options_t *user_opt = (user_options_t*) package;
    return disk_check(user_opt);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t *user_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user_opt = (user_options_t*) package;
    FCLAW_ASSERT (user_opt->is_registered);

    disk_destroy (user_opt);

    FCLAW_FREE (user_opt);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

user_options_t* disk_options_register (fclaw_app_t * app,
                                       const char *configfile)
{
    user_options_t *user_opt;
    FCLAW_ASSERT (app != NULL);

    user_opt = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user_opt);

    fclaw_app_set_attribute(app,"user",user_opt);

    return user_opt;
}

void disk_options_store (fclaw_global_t* glob, user_options_t* user_opt)
{
    fclaw_global_options_store(glob, "user", user_opt);
}

const user_options_t* disk_get_options(fclaw_global_t* glob)
{
    return (user_options_t*) fclaw_global_get_options(glob,"user");
}

#if 0
void disk_global_post_process(fclaw_options_t *fclaw_opt,
                              fclaw2d_clawpatch_options_t *clawpatch_opt,
                              user_options_t *user_opt)
{
    if (user_opt->example == 1)
        if (clawpatch_opt->mx*pow_int(2,fclaw_opt->minlevel) < 32)
        {
            fclaw_global_essentialf("The five patch mapping requires mx*2^minlevel >= 32\n");
            exit(0);
        }
}
#endif




