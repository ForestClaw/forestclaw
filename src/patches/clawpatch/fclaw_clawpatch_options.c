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

#include <fclaw_clawpatch_options.h>

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_packing.h>
#include <sc_keyvalue.h>

static sc_keyvalue_t *
kv_refinement_criterea_new()
{
    sc_keyvalue_t *kv = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv, "value",        FCLAW_REFINE_CRITERIA_VALUE);
    sc_keyvalue_set_int (kv, "difference",   FCLAW_REFINE_CRITERIA_DIFFERENCE);
    sc_keyvalue_set_int (kv, "minmax",       FCLAW_REFINE_CRITERIA_MINMAX);
    sc_keyvalue_set_int (kv, "gradient",     FCLAW_REFINE_CRITERIA_GRADIENT);
    sc_keyvalue_set_int (kv, "user",         FCLAW_REFINE_CRITERIA_USER);

    return kv;
}
static void *
clawpatch_register(fclaw_clawpatch_options_t *clawpatch_options,
                   sc_options_t * opt)
{
    int original_patch_dim = clawpatch_options->patch_dim;
    if(clawpatch_options->patch_dim == 0)
    {
        sc_options_add_int (opt, 0, "patch-dim", &clawpatch_options->patch_dim, 2,
                            "Dimension of the patch (2 or 3) [2]");
    }
    sc_options_add_int (opt, 0, "mx", &clawpatch_options->mx, 8,
                        "Number of grid cells per patch in x [8]");

    sc_options_add_int (opt, 0, "my", &clawpatch_options->my, 8,
                        "Number of grid cells per patch in y [8]");

    if(original_patch_dim == 2)
    {
        clawpatch_options->mz = 1;
    }
    else
    {
        sc_options_add_int (opt, 0, "mz", &clawpatch_options->mz, 8,
                           "Number of grid cells per patch in z [1]");
    }

    sc_options_add_int (opt, 0, "maux", &clawpatch_options->maux, 0,
                        "Number of auxilliary variables [0]");

    sc_options_add_int (opt, 0, "mbc", &clawpatch_options->mbc, 2,
                        "Number of ghost cells [2]");

    sc_options_add_int (opt, 0, "meqn", &clawpatch_options->meqn, 1,
                        "Number of equations [1]");

    sc_options_add_int (opt, 0, "rhs-fields", &clawpatch_options->rhs_fields, 0,
                        "Number of fields in rhs [0]");

    /* ---------------------- advanced options -------------------------- */
    sc_options_add_int (opt, 0, "interp_stencil_width",
                        &clawpatch_options->interp_stencil_width,
                        3, "Interpolation stencil width [3]");

    sc_options_add_bool (opt, 0, "ghost_patch_pack_aux", 
                         &clawpatch_options->ghost_patch_pack_aux,0,
                         "Pack aux. variables for parallel comm. of ghost patches [F]");

    sc_options_add_bool (opt, 0, "save-aux", 
                         &clawpatch_options->save_aux,0,
                         "Save aux variables when re-taking a time step [F]");

    /* ---------------------- vtk options -------------------------- */
    sc_options_add_int(opt, 0, "vtk-patch-threshold", 
                       &clawpatch_options->vtk_patch_threshold, 0,
                       "Number of patches to buffer before each write in vtk output. 0 means buffer all patches before writing [0]");

    fclaw_options_add_int_array(opt, 0, "vtk-aux-out",
                                &clawpatch_options->vtk_aux_out_string, "",
                                &clawpatch_options->vtk_aux_out, 0,
                                "List of aux field indexes to output in vtk. Indes are base 1, so 1 is the first aux field. "
                                "If empty, no aux variables are output []");

    /* ---------------------- hdf5 options -------------------------- */
    sc_options_add_int(opt, 0, "hdf5-patch-threshold", 
                       &clawpatch_options->hdf5_patch_threshold, 0,
                       "Number of patches to buffer before each write in hdf5 output. 0 means buffer all patches before writing [0]");
    sc_options_add_int(opt, 0, "hdf5-compression-level", 
                       &clawpatch_options->hdf5_compression_level, 5,
                       "Compression level for hdf5 output. 0 is no compression, 9 is most compression. [5]");
    fclaw_options_add_int_array(opt, 0, "hdf5-aux-out",
                                &clawpatch_options->hdf5_aux_out_string, "",
                                &clawpatch_options->hdf5_aux_out, 0,
                                "List of aux field indexes to output in hdf5. Indes are base 1, so 1 is the first aux field. "
                                "If empty, no aux variables are output []");



    /* Set verbosity level for reporting timing */
    sc_keyvalue_t *kv = clawpatch_options->kv_refinement_criteria;
    sc_options_add_keyvalue (opt, 0, "refinement-criteria", 
                             &clawpatch_options->refinement_criteria, "minmax",
                             kv, "Refinement criteria [minmax]");

    sc_options_add_int (opt, 0, "threshold-variable",
                        &clawpatch_options->threshold_variable,
                        1, "Index of variable used for tagging in [1,meqn] [1]");


    clawpatch_options->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
clawpatch_postprocess(fclaw_clawpatch_options_t *clawpatch_opt)
{
    /* Convert strings to arrays */
    fclaw_options_convert_int_array(clawpatch_opt->vtk_aux_out_string,
                                    &clawpatch_opt->vtk_aux_out,
                                    clawpatch_opt->maux);

    fclaw_options_convert_int_array(clawpatch_opt->hdf5_aux_out_string,
                                    &clawpatch_opt->hdf5_aux_out,
                                    clawpatch_opt->maux);
    return FCLAW_NOEXIT;
}

static fclaw_exit_type_t
clawpatch_check(fclaw_clawpatch_options_t *clawpatch_opt)
{
    if (clawpatch_opt->patch_dim != 2 && clawpatch_opt->patch_dim != 3)
    {
        fclaw_global_essentialf("Clawpatch error : patch_dim must be 2 or 3\n");
        return FCLAW_EXIT_ERROR;
    }

    if (2*clawpatch_opt->mbc > clawpatch_opt->mx)
    {
        fclaw_global_essentialf("Clawpatch error : 2*mbc > mx or 2*mbc > my\n");
        return FCLAW_EXIT_ERROR;
    }

    if (clawpatch_opt->interp_stencil_width/2 > clawpatch_opt->mbc)
    {
        fclaw_global_essentialf("Interpolation width is too large for number of " \
                                "ghost cells (mbc) specified.  We should have " \
                                "(width)/2 <= mbc");
    }

    /* Don't check value, in case use wants to set something */
    if (clawpatch_opt->refinement_criteria < 0 || clawpatch_opt->refinement_criteria > 4)
    {
        fclaw_global_essentialf("Clawpatch error : Default refinement criteria " \
                                "must be one of 'value','difference','minmax', " \
                                "'gradient', or 'user'.\n");
        return FCLAW_EXIT_ERROR;            
    }

    if (clawpatch_opt->vtk_patch_threshold < 0)
    {
        fclaw_global_essentialf("Clawpatch error : vtk-patch-threshold must be " \
                                "non-negative.\n");
        return FCLAW_EXIT_ERROR;            
    }

    if (clawpatch_opt->hdf5_patch_threshold < 0)
    {
        fclaw_global_essentialf("Clawpatch error : hdf5-patch-threshold must be " \
                                "non-negative.\n");
        return FCLAW_EXIT_ERROR;            
    }

    if (clawpatch_opt->hdf5_compression_level < 0 || clawpatch_opt->hdf5_compression_level > 9)
    {
        fclaw_global_essentialf("Clawpatch error : hdf5-compression-level must be " \
                                "between 0 and 9.\n");
        return FCLAW_EXIT_ERROR;            
    }

    for (int i = 0; i < clawpatch_opt->maux; i++)
    {
        if (clawpatch_opt->vtk_aux_out[i] < 0 || clawpatch_opt->vtk_aux_out[i] > clawpatch_opt->maux)
        {
            fclaw_global_essentialf("Clawpatch error : vtk-aux-out must be in the range [1,maux].\n");
            return FCLAW_EXIT_ERROR;
        }

        if (clawpatch_opt->hdf5_aux_out[i] < 0 || clawpatch_opt->hdf5_aux_out[i] > clawpatch_opt->maux)
        {
            fclaw_global_essentialf("Clawpatch error : vtk-aux-out must be in the range [1,maux].\n");
            return FCLAW_EXIT_ERROR;
        }
    }
    return FCLAW_NOEXIT;
}

fclaw_clawpatch_options_t*
fclaw_clawpatch_options_new (int dim)
{
    fclaw_clawpatch_options_t* clawpatch_options = FCLAW_ALLOC_ZERO(fclaw_clawpatch_options_t,1);

    clawpatch_options->patch_dim = dim;

    clawpatch_options->kv_refinement_criteria = kv_refinement_criterea_new();

    return clawpatch_options;
}

void
fclaw_clawpatch_options_destroy (fclaw_clawpatch_options_t *clawpatch_opt)
{
    FCLAW_FREE(clawpatch_opt->vtk_aux_out);
    FCLAW_FREE(clawpatch_opt->hdf5_aux_out);

    if(clawpatch_opt->kv_refinement_criteria != NULL)
    {
        sc_keyvalue_destroy (clawpatch_opt->kv_refinement_criteria);
    }

    FCLAW_FREE(clawpatch_opt);
}

/* ------------------------------------------------------------------------
  Generic functions - these call the functions above
  ------------------------------------------------------------------------ */

static void *
options_register(fclaw_app_t * a, void *optpkg, sc_options_t * opt)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (opt != NULL);

    fclaw_clawpatch_options_t *clawpatch_opt = 
                               (fclaw_clawpatch_options_t *) optpkg;

    return clawpatch_register(clawpatch_opt,opt);
}

static fclaw_exit_type_t
options_postprocess(fclaw_app_t * a, void *optpkg, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (optpkg != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_clawpatch_options_t *clawpatch_opt = 
                               (fclaw_clawpatch_options_t *) optpkg;

    FCLAW_ASSERT (clawpatch_opt->is_registered);

    return clawpatch_postprocess(clawpatch_opt);
}

static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_clawpatch_options_t *clawpatch_opt = 
                              (fclaw_clawpatch_options_t *) package;

    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return clawpatch_check (clawpatch_opt);
}

static void
options_destroy (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    fclaw_clawpatch_options_t *clawpatch_opt = 
                               (fclaw_clawpatch_options_t*) package;

    FCLAW_ASSERT (clawpatch_opt->is_registered);

    /* Destroy option arrays created in post-process */
    fclaw_clawpatch_options_destroy (clawpatch_opt);
}


static
const fclaw_app_options_vtable_t fclaw_clawpatch_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};


/* ---------------------------------------------------------
   Public interface to clawpatch options
   --------------------------------------------------------- */

static
fclaw_clawpatch_options_t *
fclaw_clawpatch_options_register(int dim, fclaw_app_t* app, const char* name, const char* configfile)
{

    FCLAW_ASSERT (app != NULL);

    /* allocate storage for fclaw_options */
    fclaw_clawpatch_options_t* clawpatch_options = fclaw_clawpatch_options_new(dim);

    fclaw_app_options_register (app,
                                name,
                                configfile,
                                &fclaw_clawpatch_options_vtable,
                                clawpatch_options);

    fclaw_app_set_attribute(app, name, clawpatch_options);
    return clawpatch_options;
}

fclaw_clawpatch_options_t *
fclaw_clawpatch_2d_options_register(fclaw_app_t* app, const char* name, const char* configfile)
{
    return fclaw_clawpatch_options_register(2,app,name,configfile);
}

fclaw_clawpatch_options_t *
fclaw_clawpatch_3d_options_register(fclaw_app_t* app, const char* name, const char* configfile)
{
    return fclaw_clawpatch_options_register(3,app,name,configfile);
}

fclaw_clawpatch_options_t *
fclaw_clawpatch_dim_ind_options_register(fclaw_app_t* app, const char* name, const char* configfile)
{
    return fclaw_clawpatch_options_register(0,app,name,configfile);
}

void 
fclaw_clawpatch_options_store (fclaw_global_t *glob, 
                               fclaw_clawpatch_options_t* clawpatch_options)
{
    fclaw_global_options_store(glob, "fclaw_clawpatch", clawpatch_options);
}

fclaw_clawpatch_options_t* 
fclaw_clawpatch_get_options(fclaw_global_t* glob)
{
    return (fclaw_clawpatch_options_t*) 
            fclaw_global_get_options(glob, "fclaw_clawpatch");
}
