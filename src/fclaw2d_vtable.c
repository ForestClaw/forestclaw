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

#include <fclaw2d_vtable.h>
#include <fclaw2d_clawpatch.h>

static fclaw2d_vtable_t vt;

fclaw2d_vtable_t* fclaw2d_vt()
{
    return &vt;
}

/* Initialize any settings that can be set here */
void fclaw2d_init_vtable()
{

    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

    /* ------------------------------------------------------------
      Functions below here depend on q and could be solver specific
      ------------------------------------------------------------- */

    /* These may be redefined by the user */
    /* Problem setup */
    fclaw_vt->problem_setup             = NULL;

    /* Defaults for regridding */
    fclaw_vt->after_regrid              = NULL;

    /* Diagnostics */
    fclaw_vt->patch_init_diagnostics       = NULL;
    fclaw_vt->patch_compute_diagnostics    = NULL;
    fclaw_vt->patch_gather_diagnostics     = NULL;
    fclaw_vt->patch_reset_diagnostics      = NULL;
    fclaw_vt->patch_finalize_diagnostics   = NULL;

    fclaw_vt->solver_init_diagnostics      = NULL;
    fclaw_vt->solver_compute_diagnostics    = NULL;
    fclaw_vt->solver_gather_diagnostics    = NULL;
    fclaw_vt->solver_reset_diagnostics     = NULL;
    fclaw_vt->solver_finalize_diagnostics  = NULL;

    fclaw_vt->user_init_diagnostics       = NULL;
    fclaw_vt->user_compute_diagnostics    = NULL;
    fclaw_vt->user_gather_diagnostics     = NULL;
    fclaw_vt->user_reset_diagnostics      = NULL;
    fclaw_vt->user_finalize_diagnostics   = NULL;

    /* ------------------------------------------------------------
      Metric functions - only loosely depend on solvers
      ------------------------------------------------------------- */
    fclaw_vt->metric_setup_mesh          = &fclaw2d_metric_setup_mesh;
    fclaw_vt->fort_setup_mesh            = &FCLAW2D_FORT_SETUP_MESH;

    fclaw_vt->metric_compute_area        = &fclaw2d_metric_compute_area;
    fclaw_vt->metric_area_set_ghost      = &fclaw2d_metric_area_set_ghost;

    fclaw_vt->metric_compute_normals     = &fclaw2d_metric_compute_normals;
    fclaw_vt->fort_compute_normals       = &FCLAW2D_FORT_COMPUTE_NORMALS;
    fclaw_vt->fort_compute_tangents      = &FCLAW2D_FORT_COMPUTE_TANGENTS;
    fclaw_vt->fort_compute_surf_normals  = &FCLAW2D_FORT_COMPUTE_SURF_NORMALS;
}

void fclaw2d_set_vtable()
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();
    if (fclaw_vt->metric_compute_area == &fclaw2d_metric_compute_area)
    {
        fclaw_vt->metric_area_set_ghost = &fclaw2d_metric_area_set_ghost;
    }
    else if (fclaw_vt->metric_compute_area == &fclaw2d_metric_compute_area_exact)
    {
        fclaw_vt->metric_area_set_ghost = &fclaw2d_metric_area_set_ghost_exact;
    }
}

#if 0
fclaw2d_vtable_t fclaw2d_get_vtable(fclaw2d_domain_t* domain)
{
    fclaw2d_vtable_t *vt;
    vt = (fclaw2d_vtable_t*) fclaw2d_domain_attribute_access(domain,"vtable",NULL);
    FCLAW_ASSERT(vt != NULL);
    return *vt;
}
#endif
