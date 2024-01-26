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

#include "swirlcons_user.h"

#include <fclaw_options.h>


#include <fclaw_pointer_map.h>

static void *
swirlcons_register (user_options_t *user, sc_options_t * opt)
{
    sc_options_add_int (opt, 0, "rp-solver", &user->rp_solver, 1,
                           "Conservative riemann solver (1-4) [T]");

    sc_options_add_int (opt, 0, "example", &user->example, 1,
                           "Velocity field choice (0-3) [1]");

    sc_options_add_int (opt, 0, "mapping", &user->mapping, 0,
                           "1 : u(x) > 0; 2: u(x) changes sign (1,2) [1]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 0.4,
                           "[user] Used in five-patch mapping [0.4]");

    sc_options_add_int (opt, 0, "initial-condition", &user->initial_condition, 0,
                        "[user] Initial condition : 0=non-smooth; 1=smooth [1]");

    sc_options_add_bool (opt, 0, "color-equation", &user->color_equation, 0,
                        "[user]  Solve color-equation using edge velocities [1]");

    sc_options_add_bool (opt, 0, "use-stream", &user->use_stream, 0,
                        "[user]  Use streamfunction [0]");


    fclaw_options_add_double_array (opt, 0, "center", &user->center_string, "0 0",
                                    &user->center, 2, 
                                    "Center point for bilinear mapping  [0,0]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
swirlcons_postprocess(user_options_t *user)
{
    fclaw_options_convert_double_array(user->center_string, &user->center,2);

    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
swirlcons_check (user_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
swirlcons_destroy(user_options_t *user)
{
    fclaw_options_destroy_array (user->center);
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

    return swirlcons_register(user,opt);
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
    return swirlcons_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return swirlcons_check(user);
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

    swirlcons_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

user_options_t* swirlcons_options_register (fclaw_app_t * app,
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

void swirlcons_options_store (fclaw_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"user") == NULL);
    fclaw_pointer_map_insert(glob->options, "user", user, NULL);
}

const user_options_t* swirlcons_get_options(fclaw_global_t* glob)
{
    user_options_t* user = (user_options_t*) 
                              fclaw_pointer_map_get(glob->options, "user");
    FCLAW_ASSERT(user != NULL);
    return user;   
}
