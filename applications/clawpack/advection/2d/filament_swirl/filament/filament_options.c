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

#include "filament_user.h"

static void *
filament_register (filament_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 1,
                        "[user] 0 = nomap; 1 = brick; 2 = five patch square [0]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 0.5,
                           "[user] Ratio used for squared- and pillow-disk [0.5]");

    fclaw_options_add_double_array (opt, 0, "center", &user->center_string, "0 0",
                                    &user->center, 2, 
                                    "Center point for bilinear mapping  [(0,0)]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
filament_check (filament_options_t *user)
{
    if (user->example < 0 || user->example > 2) {
        fclaw_global_essentialf ("Option --user:example must be 0, 1, or 2\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;
}


static void
filament_destroy (filament_options_t *user)
{
    /* Nothing to destroy */
}


/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    filament_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (filament_options_t*) package;

    return filament_register(user,opt);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    filament_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (filament_options_t*) package;
    return filament_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    filament_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (filament_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    filament_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    options_check,
    options_destroy,
};


/* ------------- User options access functions --------------------- */

filament_options_t* filament_options_register (fclaw_app_t * app,
                                               const char *section,
                                               const char *configfile)
{
    filament_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (filament_options_t, 1);
    fclaw_app_options_register (app,section, configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,section,user);
    return user;
}

void filament_options_store (fclaw2d_global_t* glob, filament_options_t* user)
{
    fclaw2d_global_options_store(glob, "user", user);
}

const filament_options_t* filament_get_options(fclaw2d_global_t* glob)
{
    return (filament_options_t*) fclaw2d_global_get_options(glob, "user");
}
