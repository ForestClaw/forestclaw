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

#include "cylinder_user.h"

static int s_user_options_package_id = -1;

static void *
cylinder_register (user_options_t *user_opt, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "example", &user_opt->example, 0,
                        "[user] 0 = cylinder; 1 = twisted cylinder [0]");

    sc_options_add_int (opt, 0, "initial-condition", &user_opt->initial_condition, 0,
                        "[user] Initial condition : 0=non-smooth; 1=smooth [1]");

    sc_options_add_int (opt, 0, "refine-pattern", &user_opt->refine_pattern, 0,
                           "[user] Refinement pattern [0]");

    sc_options_add_int (opt, 0, "exact-metric", &user_opt->exact_metric, 0,
                           "[user] Use exact metric [0]");

    sc_options_add_double (opt, 0, "R", &user_opt->R, 1.0,
                           "[user] Radius of cylinder [1]");

    sc_options_add_double (opt, 0, "H", &user_opt->H, 6.283185307179586,
                           "[user] Height of cylinder [1]");

    sc_options_add_double (opt, 0, "r0", &user_opt->r0, 0.1,
                           "[user] Initial radius [0.1]");

    sc_options_add_double (opt, 0, "xc0", &user_opt->xc0, 0.5,
                           "[user] Initial theta (in [0,1]) [0.5]");

    sc_options_add_double (opt, 0, "yc0", &user_opt->yc0, 0.25,
                           "[user] Initial height in [0,1]) [0.5]");

    sc_options_add_double (opt, 0, "revs-per-s", &user_opt->revs_per_s, 0.5,
                           "[user] Revolutions per second [0.5]");

    sc_options_add_double (opt, 0, "v-speed", &user_opt->v_speed, 0.5,
                           "[user] Vertical speed [0.5]");

    sc_options_add_int (opt, 0, "mapping", &user_opt->mapping, 0,
                        "[user] Mapping : 0=cylinder; 1=latlong [0]");


    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 5,
                        "[user] Clawpack version (4 or 5) [5]");

    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
cylinder_postprocess(user_options_t *user_opt)
{
    /* Not used */
    double pi = M_PI;
    user_opt->H = 2*pi*user_opt->R;
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
cylinder_check(user_options_t *user_opt)
{
    /* Nothing to check */
    return FCLAW_NOEXIT;

}

static void
cylinder_destroy(user_options_t *user_opt)
{
    /* Not used yet */
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

    return cylinder_register(user_opt,opt);
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
    return cylinder_postprocess (user_opt);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user_opt = (user_options_t*) package;

    return cylinder_check(user_opt);
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

    cylinder_destroy (user_opt);

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

user_options_t* cylinder_options_register (fclaw_app_t * app,
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

void cylinder_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const user_options_t* cylinder_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}

