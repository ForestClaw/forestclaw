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

#include "torus_user.h"

#include <fclaw_pointer_map.h>

static void *
torus_register (user_options_t *user_opt, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user_opt->example, 0,
                        "[user] 0 = torus; 1 = twisted torus [0]");

    sc_options_add_int (opt, 0, "initial-condition", &user_opt->initial_condition, 0,
                        "[user] Initial condition : 0=non-smooth; 1=smooth [1]");

    sc_options_add_double (opt, 0, "alpha", &user_opt->alpha, 0.4,
                           "[user] Ratio r/R, r=outer radius, R=inner radius " \
                           "(used for torus) [0.4]");

    sc_options_add_double (opt, 0, "beta", &user_opt->beta, 0.0,
                           "[user] beta > 0 gives variable cross section [0]");

    sc_options_add_double (opt, 0, "init-radius", &user_opt->init_radius, 0.1,
                           "[user] Initial radius [0.1]");


    fclaw_options_add_double_array (opt, 0, "theta", 
                                    &user_opt->theta_string,"0 1",&user_opt->theta,2,
                                    "[user] theta range [0,1]");    

    fclaw_options_add_double_array (opt, 0, "phi", 
                                    &user_opt->phi_string,"0 1",&user_opt->phi,2,
                                    "[user] phi range [0,1]");    

    sc_options_add_int (opt, 0, "refine-pattern", &user_opt->refine_pattern, 0,
                           "[user] Refinement pattern [0]");

    sc_options_add_double (opt, 0, "revs-per-s", &user_opt->revs_per_s, 0.5,
                           "[user] Revolutions per second [0.5]");

    sc_options_add_double (opt, 0, "cart_speed", &user_opt->cart_speed, 0.5,
                           "[user] Cartesian speed [1]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");

    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
torus_postprocess(user_options_t *user_opt)
{
    fclaw_options_convert_double_array (user_opt->phi_string, &user_opt->phi, 2);
    fclaw_options_convert_double_array (user_opt->theta_string, &user_opt->theta, 2);
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
torus_check(user_options_t *user_opt)
{
    /* Nothing to check */
    return FCLAW_NOEXIT;

}

static void
torus_destroy(user_options_t *user_opt)
{
    FCLAW_FREE (user_opt->theta);
    FCLAW_FREE (user_opt->phi);
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

    return torus_register(user_opt,opt);
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
    return torus_postprocess (user_opt);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user_opt = (user_options_t*) package;

    return torus_check(user_opt);
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

    torus_destroy (user_opt);

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

user_options_t* torus_options_register (fclaw_app_t * app,
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

void torus_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"user") == NULL);
    fclaw_pointer_map_insert(glob->options, "user", user, NULL);
}

const user_options_t* torus_get_options(fclaw2d_global_t* glob)
{
    user_options_t* user = (user_options_t*) 
                              fclaw_pointer_map_get(glob->options, "user");
    FCLAW_ASSERT(user != NULL);
    return user;   
}

