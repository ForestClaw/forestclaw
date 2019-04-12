/*
  Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "mgtest_user.h"

static int s_user_options_package_id = -1;

static void *
mgtest_register (mgtest_options_t *user, sc_options_t * opt)
{
    /* [user] User options */

    sc_options_add_int (opt, 0, "rhs-choice", &user->rhs_choice, 1,
                           "[user] RHS choice [1]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 20,
                           "alpha (used rhs_choice == 1) [0.5]");

    sc_options_add_double (opt, 0, "x0", &user->x0, 0.5,
                           "x-location of center (used in rhs_choice == 1) [0.5]");

    sc_options_add_double (opt, 0, "y0", &user->y0, 0.5,
                           "y-location of center (used in rhs_choice == 1) [0.5]");

    sc_options_add_double (opt, 0, "a", &user->a, 2,
                           "x-frequency (used in rhs_choice == 2) [2]");

    sc_options_add_double (opt, 0, "b", &user->b, 2,
                           "y-frequency of center (used in rhs_choice == 2) [2]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
mgtest_postprocess(mgtest_options_t *user)
{
    /* nothing to post-process yet ... */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
mgtest_check (mgtest_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
mgtest_destroy(mgtest_options_t *user)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    mgtest_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (mgtest_options_t*) package;

    return mgtest_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    mgtest_options_t *user = (mgtest_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return mgtest_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    mgtest_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (mgtest_options_t*) package;

    return mgtest_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    mgtest_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (mgtest_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    mgtest_destroy (user);

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

mgtest_options_t* mgtest_options_register (fclaw_app_t * app,
                                       const char *configfile)
{
    mgtest_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (mgtest_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void mgtest_options_store (fclaw2d_global_t* glob, mgtest_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const mgtest_options_t* mgtest_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (mgtest_options_t*) fclaw_package_get_options(glob, id);    
}
