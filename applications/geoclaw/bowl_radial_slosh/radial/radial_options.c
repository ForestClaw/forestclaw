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

#include "radial_user.h"

static void *
radial_register (radial_user_options_t *user, sc_options_t * opt)
{
    user->is_registered = 1;
    return NULL;
}

static void
radial_destroy(radial_user_options_t *user)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    radial_user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (radial_user_options_t*) package;

    return radial_register(user,opt);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    radial_user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (radial_user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    radial_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    NULL,
    options_destroy
};

radial_user_options_t* radial_options_register (fclaw_app_t * app,
                                        const char *section,
                                        const char *configfile)
{
    radial_user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (radial_user_options_t, 1);
    fclaw_app_options_register (app, section, configfile, &options_vtable_user,
                                user);

    return user;
}

void radial_options_store (fclaw2d_global_t* glob, radial_user_options_t* user)
{
    fclaw2d_global_options_store(glob, "radial-user", user);
}

radial_user_options_t* radial_get_options(fclaw2d_global_t* glob)
{
    return (radial_user_options_t*) fclaw2d_global_get_options(glob, "radial-user");
}

