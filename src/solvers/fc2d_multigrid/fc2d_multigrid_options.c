/*
Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc2d_multigrid_options.h"

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw_options.h>
#include <fclaw_package.h>

static int s_multigrid_options_package_id = -1;

static void*
multigrid_register (fc2d_multigrid_options_t* mg_opt, sc_options_t * opt)
{

#if 0
    sc_options_add_int (opt, 0, "mfields", &mg_opt->mfields, 0,
                        "The number of fields in solution [1]");
#endif                        


    /* Array of NumFaces=4 values */
    fclaw_options_add_int_array (opt, 0, "boundary_conditions", 
                                 &mg_opt->bc_cond_string, "1 1 1 1",
                                 &mg_opt->boundary_conditions, 4,
                                 "[multigrid] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "ascii-out", &mg_opt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &mg_opt->vtk_out, 0,
                           "Output VTK formatted data [F]");

    sc_options_add_bool (opt, 0, "mg-prec", &mg_opt->mg_prec, 1,
                           "Use multigrid preconditioner [T]");

    sc_options_add_int (opt, 0, "max-it", &mg_opt->max_it, 10000,
                           "Max iterations for BiCGStab solver. [10000]");

    sc_options_add_double (opt, 0, "tol", &mg_opt->tol, 1e-12,
                           "Tolerance for BiCGStab solver. [1e-12]");

    sc_options_add_int (opt, 0, "max-levels", &mg_opt->max_levels, 0,
                           "The max number of levels in GMG cycle. 0 means no limit. [0]");

    sc_options_add_double (opt, 0, "patches-per-proc", &mg_opt->patches_per_proc, 0,
                           "Lowest level is guaranteed to have at least this number of "
                           "patches per processor. [0]");

    sc_options_add_int (opt, 0, "pre-sweeps", &mg_opt->pre_sweeps, 1,
                           "Number of sweeps on down cycle [1]");

    sc_options_add_int (opt, 0, "post-sweeps", &mg_opt->post_sweeps, 1,
                           "Number of sweeps on up cycle [1]");

    sc_options_add_int (opt, 0, "mid-sweeps", &mg_opt->mid_sweeps, 1,
                           "Number of sweeps inbetween up and down [1]");

    sc_options_add_int (opt, 0, "coarse-sweeps", &mg_opt->coarse_sweeps, 1,
                           "Number of sweeps on coarse level [1]");

    sc_options_add_string (opt, 0, "cycle-type", &mg_opt->cycle_type, "V",
                           "Cycle type [V]");

    sc_options_add_double (opt, 0, "patch-bcgs-tol", &mg_opt->patch_bcgs_tol, 1e-1,
                           "Tolerance for patch-based bcgs solver [1e-1]");

    sc_options_add_int (opt, 0, "patch-bcgs-max-it", &mg_opt->patch_bcgs_max_it, 1000,
                           "Max allowed iterations for patch-based bcgs solver [1000]");


    /* Used by starpatch only */
    sc_options_add_string (opt, 0, "patch-solver-type", &mg_opt->patch_solver_type, "BCGS",
                           "Patch solver type. Can either be BCGS or FFT [BCGS]");

    /* Set operator type (starpatch, fivepoint) */
    sc_keyvalue_t *kv_op = mg_opt->kv_patch_operator = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv_op, "starpatch",  STARPATCH);     /* Uses patch-solver-type */
    sc_keyvalue_set_int (kv_op, "fivepoint",  FIVEPOINT);     /* Uses FFT or BICG */
    sc_keyvalue_set_int (kv_op, "varpoisson",  VARPOISSON);   /* Uses BICG */
    sc_keyvalue_set_int (kv_op, "user_operator",  USER_OPERATOR);   /* Uses BICG */
    sc_options_add_keyvalue (opt, 0, "patch_operator", &mg_opt->patch_operator,
                             "fivepoint", kv_op, "Set patch operator type [fivepoint]");

    /* Set solver type (FFT, BICG) */
    sc_keyvalue_t *kv_s = mg_opt->kv_patch_solver = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv_s, "bicg", BICG);
    sc_keyvalue_set_int (kv_s, "fft",  FFT);     
    sc_keyvalue_set_int (kv_s, "user_solver",  USER_SOLVER);     
    sc_options_add_keyvalue (opt, 0, "patch_solver", &mg_opt->patch_solver,
                             "bicg", kv_s, "Set patch solver type [BICG]");

    mg_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
multigrid_postprocess (fc2d_multigrid_options_t * mg_opt)
{
    fclaw_options_convert_int_array (mg_opt->bc_cond_string, 
                                     &mg_opt->boundary_conditions,4);
    
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
multigrid_check(fc2d_multigrid_options_t *mg_opt,
                 fclaw2d_clawpatch_options_t *clawpatch_opt)
{
    return FCLAW_NOEXIT;
}

static
void multigrid_destroy (fc2d_multigrid_options_t * mg_opt)
{
    fclaw_options_destroy_array (mg_opt->boundary_conditions);

    FCLAW_ASSERT (mg_opt->kv_patch_operator != NULL);
    sc_keyvalue_destroy (mg_opt->kv_patch_operator);

    FCLAW_ASSERT (mg_opt->kv_patch_solver != NULL);
    sc_keyvalue_destroy (mg_opt->kv_patch_solver);
}

/* ------------------------------------------------------
   Generic calls to options handling;  each calls 
   clawpack-specific options call back
   ------------------------------------------------------ */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc2d_multigrid_options_t *mg_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    mg_opt = (fc2d_multigrid_options_t*) package;

    return multigrid_register(mg_opt,opt);
}


static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_multigrid_options_t *mg_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    mg_opt = (fc2d_multigrid_options_t*) package;
    FCLAW_ASSERT (mg_opt->is_registered);

    return multigrid_postprocess (mg_opt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_multigrid_options_t *mg_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    mg_opt = (fc2d_multigrid_options_t*) package;
    FCLAW_ASSERT (mg_opt->is_registered);

    clawpatch_opt = (fclaw2d_clawpatch_options_t *)
        fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return multigrid_check(mg_opt,clawpatch_opt);    
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_multigrid_options_t *mg_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    mg_opt = (fc2d_multigrid_options_t*) package;
    FCLAW_ASSERT (mg_opt->is_registered);

    multigrid_destroy (mg_opt);

    FCLAW_FREE (mg_opt);
}

static const fclaw_app_options_vtable_t multigrid_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc2d_multigrid_options_t*  fc2d_multigrid_options_register (fclaw_app_t * app,
                                                              const char *configfile)
{
    fc2d_multigrid_options_t *mg_opt;

    FCLAW_ASSERT (app != NULL);

    mg_opt = FCLAW_ALLOC (fc2d_multigrid_options_t, 1);
    fclaw_app_options_register (app, "multigrid", configfile,
                                &multigrid_options_vtable, mg_opt);
    
    fclaw_app_set_attribute(app,"multigrid",mg_opt);
    return mg_opt;
}

fc2d_multigrid_options_t* fc2d_multigrid_get_options(fclaw2d_global_t *glob)
{
    int id = s_multigrid_options_package_id;
    return (fc2d_multigrid_options_t*) fclaw_package_get_options(glob,id);
}

void fc2d_multigrid_options_store (fclaw2d_global_t* glob, fc2d_multigrid_options_t* mg_opt)
{
    int id = fclaw_package_container_add_pkg(glob,mg_opt);
    s_multigrid_options_package_id = id;
}
