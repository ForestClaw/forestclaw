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

#include "replicated_user.h"

static int s_user_options_package_id = -1;

static void *
replicated_register (user_options_t *user_opt, sc_options_t * opt)
{
    sc_options_add_int (opt, 0, "example", &user_opt->example, 0,
                        "[user] 0 - multi-block; 1 - single block [0]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 4,
                        "[user] Clawpack version (4 or 5) [5]");

    sc_options_add_int (opt, 0, "replicate-factor", &user_opt->replicate_factor, 1,
                        "[user] Replication factor [-1]");

    sc_options_add_int (opt, 0, "minlevel-base", &user_opt->minlevel_base, 4,
                        "[user] Block 1x1 minimum refinement [4]");

    sc_options_add_int (opt, 0, "maxlevel-base", &user_opt->maxlevel_base, 7,
                        "[user] Block 1x1 maximum refinement [7]");


    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t 
replicated_postprocess(user_options_t *user, fclaw_options_t *fclaw_opt)
{

    /* Check to make sure replicate factor is power of 2. 
       We do it here, since options_check only happens after post-process. */
    int pos = 0;    
    while (user->replicate_factor != (1 << pos))
    {
        pos++;
        if (pos > 32)
        {
            fclaw_global_essentialf("Replicated_postprocess : Replicate " \
                                    "factor should be a low power of 2 (< 2**32)\n");
            return FCLAW_EXIT_QUIET;
        }
    }

    fclaw_opt->ax = 0;
    fclaw_opt->ay = 0;
    fclaw_opt->bx = user->replicate_factor;
    fclaw_opt->by = user->replicate_factor;

    fclaw_opt->minlevel = user->minlevel_base;
    fclaw_opt->maxlevel = user->maxlevel_base;



    if (user->example == 0)
    {
        /* 1x1 block;  use replicate factor to dimension domain */
        fclaw_opt->mi = 1;
        fclaw_opt->mj = 1;

        /* We need to increase the refinement level when we increase domain size. 
           Code below computes pos=log2(replicate_factor) */
        fclaw_opt->minlevel += pos;
        fclaw_opt->maxlevel += pos;
    }
    else
    {
        /* NxN block */
        fclaw_opt->mi = user->replicate_factor;
        fclaw_opt->mj = user->replicate_factor;
    }

    fclaw_opt->periodic_x = 1;
    fclaw_opt->periodic_y = 1;
    fclaw_opt->smooth_level = fclaw_opt->maxlevel-1;
    return FCLAW_NOEXIT;
}

static void
replicated_destroy (user_options_t *user)
{
    /* Nothing to destroy */
}


/* --------------- Generic option handling routines that call above routines ---------- */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (user_options_t*) package;

    return replicated_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;
    fclaw_options_t *fclaw_opt = 
                 (fclaw_options_t*) fclaw_app_get_attribute(app,"Options",NULL);

    return replicated_postprocess(user, fclaw_opt);
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

    replicated_destroy (user);

    FCLAW_FREE (user);
}


static const
fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    NULL,
    options_destroy
};

/* ---------------------------- Public interface -------------------------------------- */

user_options_t* replicated_options_register (fclaw_app_t * app,
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

void replicated_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

const user_options_t* replicated_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}

