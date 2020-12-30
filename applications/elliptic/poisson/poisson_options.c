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

#include "poisson_user.h"

static int s_user_options_package_id = -1;

static void *
poisson_register (poisson_options_t *user, sc_options_t * opt)
{
    /* [user] User options */

    sc_options_add_int (opt, 0, "example", &user->example, 1,
                           "[user] Example [1]");

    sc_options_add_int (opt, 0, "beta-choice", &user->beta_choice, 1,
                           "[user] Beta choice [0]");

    sc_options_add_double (opt, 0, "alpha", &user->alpha, 20,
                           "alpha (used example == 1) [0.5]");

    sc_options_add_double (opt, 0, "x0", &user->x0, 0.5,
                           "x-location of center (used in example == 1) [0.5]");

    sc_options_add_double (opt, 0, "y0", &user->y0, 0.5,
                           "y-location of center (used in example == 1) [0.5]");

    sc_options_add_double (opt, 0, "a", &user->a, 2,
                           "x-frequency (used in example == 2) [2]");

    sc_options_add_double (opt, 0, "b", &user->b, 2,
                           "y-frequency of center (used in example == 2) [2]");

    /* For disk */
    sc_options_add_double (opt, 0, "eps-disk", &user->eps_disk, 0.015625,
                           "eps_disk defines sharpness of the interface [0.015625]");

    /* Parameters for polar plots */
    sc_options_add_int (opt, 0, "m-polar", &user->m_polar, 1,
                           "m_polar number of polar flowers [1]");

    fclaw_options_add_double_array (opt, 0, "x0-polar", &user->x0_polar_string,"0",
        &user->x0_polar,user->m_polar,"(polar plots) x location of center [0]");

    fclaw_options_add_double_array (opt, 0, "y0-polar", &user->y0_polar_string,"0",
        &user->y0_polar,user->m_polar, "(polar plots) y location of center [0]");

    fclaw_options_add_double_array (opt, 0, "r0-polar", &user->r0_polar_string, "0.25",
        &user->r0_polar,user->m_polar, "(polar plots) r0-polar radius for polar disks [0.25]");

    fclaw_options_add_double_array (opt, 0, "r1-polar", &user->r1_polar_string, "0.35",
        &user->r1_polar,user->m_polar,"(polar plots) r1-polar length of polar petals [0.35]");

    fclaw_options_add_int_array (opt, 0, "n-polar", &user->n_polar_string, "4",
        &user->n_polar,user->m_polar,"(polar plots) n_polar number of polar petals [4]");

#if 0
    /* Set operator type (star, fivepoint) */
    sc_keyvalue_t *kv = user->kv_patch_operator = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv, "starpatch",  STARPATCH);
    sc_keyvalue_set_int (kv, "fivepoint",  FIVEPOINT);
    sc_keyvalue_set_int (kv, "varpoisson",  VARPOISSON);
    sc_options_add_keyvalue (opt, 0, "patch_operator", &user->patch_operator,
                             "fivepoint", kv, "[user] Set operator type [fivepoint]");
#endif                             

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
poisson_postprocess(poisson_options_t *user)
{
    fclaw_options_convert_double_array (user->x0_polar_string,&user->x0_polar,user->m_polar);
    fclaw_options_convert_double_array (user->y0_polar_string,&user->y0_polar,user->m_polar);
    fclaw_options_convert_double_array (user->r0_polar_string,&user->r0_polar,user->m_polar);
    fclaw_options_convert_double_array (user->r1_polar_string,&user->r1_polar,user->m_polar);
    fclaw_options_convert_int_array (user->n_polar_string, &user->n_polar, user->m_polar);

    /* nothing to post-process yet ... */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
poisson_check (poisson_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
poisson_destroy(poisson_options_t *user)
{
    fclaw_options_destroy_array (user->x0_polar);
    fclaw_options_destroy_array (user->y0_polar);
    fclaw_options_destroy_array (user->r0_polar);
    fclaw_options_destroy_array (user->r1_polar);
    fclaw_options_destroy_array (user->n_polar);

#if 0
    FCLAW_ASSERT (user->kv_patch_operator != NULL);
    sc_keyvalue_destroy (user->kv_patch_operator);
#endif    
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    poisson_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (poisson_options_t*) package;

    return poisson_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    poisson_options_t *user = (poisson_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return poisson_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    poisson_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (poisson_options_t*) package;

    return poisson_check(user);
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    poisson_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (poisson_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    poisson_destroy (user);

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

poisson_options_t* poisson_options_register (fclaw_app_t * app,
                                       const char *configfile)
{
    poisson_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (poisson_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void poisson_options_store (fclaw2d_global_t* glob, poisson_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const poisson_options_t* poisson_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (poisson_options_t*) fclaw_package_get_options(glob, id);    
}
