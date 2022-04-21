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

#include "shockbubble_user.h"

static int s_user_options_package_id = -1;

static void *
shockbubble_register (user_options_t* user, sc_options_t * opt)
{

    /* [user] User options */
    sc_options_add_int (opt, 0, "idisc", &user->idisc,    2, "[user] idisc [2]");
    sc_options_add_double (opt, 0, "gamma", &user->gamma, 1.4, "[user] gamma [1.4]");

    sc_options_add_double (opt, 0, "x0",    &user->x0,    0.5, "[user] x0 [0.5]");
    sc_options_add_double (opt, 0, "y0",    &user->y0,    0.0, "[user] y0 [0.0]");
    sc_options_add_double (opt, 0, "r0",    &user->r0,    0.2, "[user] r0 [0.2]");
    sc_options_add_double (opt, 0, "rhoin", &user->rhoin, 0.1, "[user] rhoin [0.1]");
    sc_options_add_double (opt, 0, "pinf",  &user->pinf,  5.0, "[user] pinf [5.0x]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 4,
                           "Clawpack_version (4 or 5) [4]");

    sc_options_add_bool (opt, 0, "cuda", &user->cuda, 0,
                           "Use cudaclaw [F]");   

    user->is_registered = 1;
    return NULL;
}

static void
shockbubble_destroy(user_options_t *user)
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

    return shockbubble_register(user,opt);
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

    shockbubble_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    NULL,
    options_destroy
};

/* ------------- User options access functions --------------------- */

user_options_t* shockbubble_options_register (fclaw_app_t * app,
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

void shockbubble_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

user_options_t* shockbubble_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}


