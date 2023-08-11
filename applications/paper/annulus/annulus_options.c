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

#include "annulus_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw_pointer_map.h>

static void *
annulus_register(user_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user->example, 0,
                        "[user] 0=rigid body rot.; 1=other velocity [0]");

    sc_options_add_int (opt, 0, "refine-pattern", &user->refine_pattern, 0,
                        "[user] 0 = constant theta; 1 = constant_r [0]");

    sc_options_add_double (opt, 0, "revs-per-s", &user->revs_per_s, 0.5,
                           "[user] Revolutions per second [0.5]");

    sc_options_add_double (opt, 0, "cart_speed", &user->cart_speed, 1.,
                           "[user] Cartesian speed [1]");

    sc_options_add_double (opt, 0, "amplitude", &user->amplitude, 0.,
                           "[user] Amplitude [0.05]");

    sc_options_add_double (opt, 0, "freq", &user->freq, 0.,
                           "[user] freq [1]");

    sc_options_add_double (opt, 0, "init_radius", &user->init_radius, 0.125,
                           "[user] Initial radius used in initial conditions [0.125]");

    sc_options_add_double (opt, 0, "beta", &user->beta, 0.0,
                           "[user] Inner radius of annulus [0.4]");

    fclaw_options_add_double_array (opt, 0, "theta", 
                                    &user->theta_string,"0 1",&user->theta,2,
                                    "[user] theta range [0,1]");

    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 4,
                        "[user] Clawpack version (4 or 5) [4]");

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
annulus_postprocess(user_options_t *user)
{
    fclaw_options_convert_double_array (user->theta_string, &user->theta, 2);

    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
annulus_check(user_options_t *user)
{
    if (user->example < 0 || user->example > 5)
    {
        fclaw_global_essentialf("example must be in [0-4]\n");
        return FCLAW_EXIT_QUIET;        
    }
    return FCLAW_NOEXIT;

}


static void
annulus_destroy (user_options_t *user)
{
    FCLAW_FREE (user->theta);
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

    return annulus_register(user,opt);
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
    return annulus_postprocess (user_opt);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user_opt = (user_options_t*) package;

    return annulus_check(user_opt);
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

    annulus_destroy (user);

    FCLAW_FREE (user);
}


static const
fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ------------- User options access functions --------------------- */

user_options_t* annulus_options_register (fclaw_app_t * app,
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

void annulus_options_store (fclaw_global_t* glob, user_options_t* user)
{
    fclaw_global_options_store(glob, "user", user);
}

const user_options_t* annulus_get_options(fclaw_global_t* glob)
{
    return (user_options_t*) fclaw_pointer_map_get(glob->options,"user");
}

