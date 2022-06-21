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

#include "swirl_user.h"

#include <fclaw_pointer_map.h>

static void *
swirl_register (user_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_double (opt, 0, "period", &user->period, 4,
                           "Period of the flow field [4]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                           "Clawpack_version (4 or 5) [5]");

    sc_options_add_bool (opt, 0, "cuda", &user->cuda, 0,
                           "Use cudaclaw5 [F]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
swirl_postprocess(user_options_t *user)
{
    /* nothing to post-process yet ... */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
swirl_check (user_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
swirl_destroy(user_options_t *user)
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

    return swirl_register(user,opt);
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
    return swirl_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return swirl_check(user);
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

    swirl_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};

/* --------------------- Public interface access functions ---------------------------- */

user_options_t* swirl_options_register (fclaw_app_t * app,
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

void swirl_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"user") == NULL);
    fclaw_pointer_map_insert(glob->options, "user", user, NULL);
}

const user_options_t* swirl_get_options(fclaw2d_global_t* glob)
{
    user_options_t* user = (user_options_t*) 
                              fclaw_pointer_map_get(glob->options, "user");
    FCLAW_ASSERT(user != NULL);
    return user;   
}
