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

#include "phasefield_user.h"

static int s_user_options_package_id = -1;

static void *
phasefield_register (phasefield_options_t *user, sc_options_t * opt)
{
    /* [user] User options */

    sc_options_add_int (opt, 0, "example", &user->example, 1,
                           "[user] Example [1]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 20,
                           "alpha (used example == 1) [0.5]");

#if 0
    read(10,*) S
    read(10,*) alpha
    read(10,*) m
    read(10,*) xi
    read(10,*) k
    read(10,*) gamma
    read(10,*) r0
#endif    
    sc_options_add_double (opt, 0, "S", &user->S, 0.5,
                           "Stefan number (dimensionless undercooling) [0.5]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 400,
                           "Phasefield parameter alpha [400]");

    sc_options_add_double (opt, 0, "m", &user->m, 0.035,
                           "Phasefield parameter m [0.035]");

    sc_options_add_double (opt, 0, "xi", &user->xi, 0.00625,
                           "Dimensionless interface thickness [0.00625]");

    sc_options_add_double (opt, 0, "k", &user->k, 4.0,
                           "Mode number [4]");

    sc_options_add_double (opt, 0, "gamma", &user->gamma, 0.015,
                           "Magnitude of the anisotropy [0.015]");

    sc_options_add_double (opt, 0, "r0", &user->r0, 0.1,
                           "Radius of initial crystal [0.1]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
phasefield_postprocess(phasefield_options_t *user)
{
    /* nothing to post-process yet ... */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
phasefield_check (phasefield_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
phasefield_destroy(phasefield_options_t *user)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    phasefield_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (phasefield_options_t*) package;

    return phasefield_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    phasefield_options_t *user = (phasefield_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return phasefield_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    phasefield_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (phasefield_options_t*) package;

    return phasefield_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    phasefield_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (phasefield_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    phasefield_destroy (user);

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

phasefield_options_t* phasefield_options_register (fclaw_app_t * app,
                                       const char *configfile)
{
    phasefield_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (phasefield_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void phasefield_options_store (fclaw2d_global_t* glob, phasefield_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const phasefield_options_t* phasefield_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (phasefield_options_t*) fclaw_package_get_options(glob, id);    
}
