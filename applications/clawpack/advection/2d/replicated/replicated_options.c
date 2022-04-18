/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

    sc_options_add_int (opt, 0, "replicate-factor", &user_opt->replicate_factor, 1,
                        "[user] Replication factor [-1]");

    sc_options_add_int (opt, 0, "minlevel-base", &user_opt->minlevel_base, 4,
                        "[user] Block 1x1 minimum refinement [4]");

    sc_options_add_int (opt, 0, "maxlevel-base", &user_opt->maxlevel_base, 7,
                        "[user] Block 1x1 maximum refinement [7]");

    /* [user] User options */
    sc_options_add_double (opt, 0, "uvel", &user_opt->uvel, 1,
                           "[user] Velocity in x direction [1]");

    sc_options_add_double (opt, 0, "vvel", &user_opt->vvel, 1,
                           "[user] Velocity in y direction [1]");

    sc_options_add_double (opt, 0, "revs-per-s", &user_opt->revs_per_s, 0.5,
                           "[user] Revolutions per second [0.5]");

    sc_options_add_int (opt, 0, "claw-version", &user_opt->claw_version, 4,
                        "[user] Clawpack version (4 or 5) [4]");

    user_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t 
replicated_postprocess(user_options_t *user, fclaw_options_t *fclaw_opt)
{

    /* Do post-processing before printing out summary */
    return FCLAW_EXIT_QUIET;
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

void replicated_global_post_process(fclaw_options_t *fclaw_opt, 
                                    fclaw2d_clawpatch_options_t *clawpatch_opt,
                                    user_options_t *user_opt)
{
    /* This routine will be called to do any more global  post-processing */
    FCLAW_ASSERT(user_opt->replicate_factor >= 0);

    fclaw_opt->ax = 0;
    fclaw_opt->ay = 0;
    fclaw_opt->bx = user_opt->replicate_factor;
    fclaw_opt->by = user_opt->replicate_factor;

    fclaw_opt->minlevel = user_opt->minlevel_base;
    fclaw_opt->maxlevel = user_opt->maxlevel_base;

    fclaw_opt->periodic_x = 1;
    fclaw_opt->periodic_y = 1;
    fclaw_opt->smooth_level = fclaw_opt->maxlevel;

    switch (user_opt->example)
    {
        case 0:
        {
            /* In this case, we mimic the multiblock behavior in a single block */
            /* Find p so that user_factor = 2**p.  Require p < 32 */
            int p = 0;    
            while (user_opt->replicate_factor != (1 << p))
            {
                p++;
                if (p > 32)
                {
                    fclaw_global_essentialf("Replicated_postprocess : Replicate " \
                                            "factor should be a low power of 2 (< 2**32)\n");
                    return exit(1);
                }
            }

            /* 1x1 block;  use replicate factor to dimension domain */
            fclaw_opt->mi = 1;
            fclaw_opt->mj = 1;

            /* We need to increase the refinement level when we increase 
               domain size.  */
            fclaw_opt->minlevel += p;
            fclaw_opt->maxlevel += p;
        }
        break;    
    
        default:
        {
            /* Multiblock case : NxN block */
            fclaw_opt->mi = user_opt->replicate_factor;
            fclaw_opt->mj = user_opt->replicate_factor;            
        }
        break;
    }    
}

