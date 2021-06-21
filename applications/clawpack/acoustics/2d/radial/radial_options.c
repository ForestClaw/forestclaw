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

#include "radial_user.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

static int s_user_options_package_id = -1;

static void *
radial_register (user_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int    (opt, 0, "example", &user->example, 0, "[user] 0 - nomap; " \
                                   "1 - 5-patch disk; 2-pillow-disk [0]");
    sc_options_add_double (opt, 0, "rho",     &user->rho, 1, "[user] rho [1]");
    sc_options_add_double (opt, 0, "bulk",    &user->bulk, 4, "[user] bulk modulus [4]");
    sc_options_add_double (opt, 0, "alpha",   &user->alpha, 0.5, "[user] alpha (for 5-patch map) [0.5]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
radial_postprocess (user_options_t *user)
{

    //user->cc = sqrt(user->bulk/user->rho);
    //user->zz = user->rho*user->cc;

    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
radial_check (user_options_t *user,
             fclaw_options_t *fclaw_opt,
             fclaw2d_clawpatch_options_t *clawpatch_opt)
{
    if (user->example < 0 || user->example > 2) {
        fclaw_global_essentialf ("Option --user:example must be 0 or 1\n");
        return FCLAW_EXIT_QUIET;
    }
    if (user->example > 0)
        if (clawpatch_opt->mx*pow_int(2,fclaw_opt->minlevel) < 32)
        {
            fclaw_global_essentialf("The five patch mapping requires mx*2^minlevel >= 32\n");
            return FCLAW_EXIT_QUIET;
        }
    return FCLAW_NOEXIT;
}

static void
radial_destroy(user_options_t *user)
{
    /* Nothing to destroy */
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (user_options_t*) package;

    return radial_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    user_options_t *user = (user_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return radial_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;
    fclaw_options_t *fclaw_opt = 
                 (fclaw_options_t*) fclaw_app_get_attribute(app,"Options",NULL);

    fclaw2d_clawpatch_options_t *clawpatch_opt = 
                 (fclaw2d_clawpatch_options_t*)  fclaw_app_get_attribute(app,"clawpatch",NULL);

    return radial_check(user,fclaw_opt, clawpatch_opt);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    radial_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};

/* ------------------------- ... and here ---------------------------- */

user_options_t* radial_options_register (fclaw_app_t * app,
                            const char *configfile)
{
    user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void radial_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

user_options_t* radial_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}

