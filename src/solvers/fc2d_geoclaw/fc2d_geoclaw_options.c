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

#include "fc2d_geoclaw_options.h"
#include <fclaw_options.h>
#include <fclaw_package.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static int s_geoclaw_options_package_id = -1;

static void*
geoclaw_register (fc2d_geoclaw_options_t *geo_opt, sc_options_t * opt)
{
    fclaw_options_add_int_array (opt, 0, "order", &geo_opt->order_string,
                               "2 2", &geo_opt->order, 2,
                               "[geoclaw] Normal and transverse orders [2 2]");

    sc_options_add_int (opt, 0, "mcapa", &geo_opt->mcapa, -1,
                        "[geoclaw] Location of capacity function in aux array [-1]");

    sc_options_add_bool (opt, 0, "src_term", &geo_opt->src_term, 0,
                         "[geoclaw] Source term option [F]");

    sc_options_add_bool (opt, 0, "use_fwaves", &geo_opt->use_fwaves, 0,
                         "[geoclaw] Use fwaves flux-form [F]");


    sc_options_add_int (opt, 0, "mwaves", &geo_opt->mwaves, 1,
                        "[geoclaw] Number of waves [1]");

    fclaw_options_add_int_array (opt, 0, "mthlim", &geo_opt->mthlim_string, NULL,
                                 &geo_opt->mthlim, geo_opt->mwaves,
                                 "[geoclaw] Waves limiters (one entry per wave; " \
                                 "values 0-4) [NULL]");

    /* Array of NumFaces many values */
    fclaw_options_add_int_array (opt, 0, "mthbc", &geo_opt->mthbc_string, "1 1 1 1",
                                 &geo_opt->mthbc, 4,
                                 "[clawpack46] Physical boundary condition type [1 1 1 1]");

    /* Coarsen criteria */
    sc_options_add_double (opt, 0, "dry_tolerance_c", &geo_opt->dry_tolerance_c, 1.0,
                           "[geoclaw] Coarsen criteria: Dry tolerance [1.0]");

    sc_options_add_double (opt, 0, "wave_tolerance_c", &geo_opt->wave_tolerance_c, 1.0,
                           "[geoclaw] Coarsen criteria: Wave tolerance [1.0]");

    sc_options_add_int (opt, 0, "speed_tolerance_entries_c",
                        &geo_opt->speed_tolerance_entries_c, 1,
                        "[geoclaw] Coarsen criteria: Number of speed tolerance entries [1]");

    fclaw_options_add_double_array (opt, 0, "speed_tolerance_c",
                                    &geo_opt->speed_tolerance_c_string, NULL,
                                    &geo_opt->speed_tolerance_c, 
                                    geo_opt->speed_tolerance_entries_c,
                                    "[geoclaw] Coarsen criteria: speed tolerance [NULL]");

    sc_options_add_int (opt, 0, "mbathy", &geo_opt->mbathy, 1,
                        "[geoclaw] Location of bathymetry in aux array [1]");

    sc_options_add_bool (opt, 0, "ascii-out", &geo_opt->ascii_out,1,
                         "Output ascii files for post-processing [T]");

    geo_opt->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t 
geoclaw_check (fc2d_geoclaw_options_t *geo_opt)
{
    geo_opt->method[0] = 0;  /* Time stepping is controlled outside of clawpack */

    geo_opt->method[1] = geo_opt->order[0];
    geo_opt->method[2] = geo_opt->order[1];
    geo_opt->method[3] = 0;  /* No verbosity allowed in fortran subroutines */
    geo_opt->method[4] = geo_opt->src_term;
    geo_opt->method[5] = geo_opt->mcapa;

    return FCLAW_NOEXIT;    /* Nothing can go wrong here! */
}


static fclaw_exit_type_t
geoclaw_postprocess (fc2d_geoclaw_options_t * geo_opt)
{
    fclaw_options_convert_int_array (geo_opt->mthbc_string, &geo_opt->mthbc,4);

    fclaw_options_convert_int_array (geo_opt->mthlim_string, &geo_opt->mthlim,
                                     geo_opt->mwaves);
    fclaw_options_convert_int_array (geo_opt->order_string, &geo_opt->order,
                                     2);
    fclaw_options_convert_double_array (geo_opt->speed_tolerance_c_string,
                                        &geo_opt->speed_tolerance_c,
                                        geo_opt->speed_tolerance_entries_c);

    return FCLAW_NOEXIT;
}

static void
geoclaw_destroy (fc2d_geoclaw_options_t * geo_opt)
{
    fclaw_options_destroy_array (geo_opt->mthbc);
    fclaw_options_destroy_array (geo_opt->order);
    fclaw_options_destroy_array (geo_opt->mthlim);
    fclaw_options_destroy_array (geo_opt->speed_tolerance_c);
}

/* ------------------------------------------------------------------------
  Generic functions - these call the functions above
  ------------------------------------------------------------------------ */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc2d_geoclaw_options_t *geo_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    geo_opt = (fc2d_geoclaw_options_t*) package;

    return geoclaw_register(geo_opt, opt);
}

static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_geoclaw_options_t *geo_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    geo_opt = (fc2d_geoclaw_options_t*) package;
    FCLAW_ASSERT(geo_opt->is_registered != 0);


    return geoclaw_check(geo_opt);
}    

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_geoclaw_options_t *geo_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    geo_opt = (fc2d_geoclaw_options_t*) package;
    FCLAW_ASSERT (geo_opt->is_registered);

    return geoclaw_postprocess (geo_opt);
}


static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_geoclaw_options_t *geo_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == 0);

    geo_opt = (fc2d_geoclaw_options_t*) package;
    FCLAW_ASSERT (geo_opt->is_registered);

    geoclaw_destroy (geo_opt);
    FCLAW_FREE (geo_opt);
}

static const fclaw_app_options_vtable_t geoclaw_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc2d_geoclaw_options_t*  
fc2d_geoclaw_options_register (fclaw_app_t * app,
                               const char *section,
                               const char *configfile)
{
    fc2d_geoclaw_options_t *geo_opt;

    FCLAW_ASSERT (app != NULL);

    geo_opt = FCLAW_ALLOC (fc2d_geoclaw_options_t, 1);
    fclaw_app_options_register (app, section, configfile,
                                &geoclaw_options_vtable, geo_opt);

    fclaw_app_set_attribute(app, section, geo_opt);

    return geo_opt;
}

fc2d_geoclaw_options_t* fc2d_geoclaw_get_options(fclaw2d_global_t *glob)
{
    int id = s_geoclaw_options_package_id;
    return (fc2d_geoclaw_options_t*)  fclaw_package_get_options(glob, id);
}

void fc2d_geoclaw_options_store (fclaw2d_global_t* glob, 
                               fc2d_geoclaw_options_t* geo_opt)
{
    int id; 

    /* Don't register a package more than once */
    FCLAW_ASSERT(s_geoclaw_options_package_id == -1);
    id = fclaw_package_container_add_pkg(glob,geo_opt);
    s_geoclaw_options_package_id = id;
}



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
