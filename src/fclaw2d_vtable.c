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

#include "fclaw2d_vtable.h"
#include "fclaw2d_regrid_default.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* Initialize any settings that can be set here */
void fclaw2d_init_vtable(fclaw2d_vtable_t *vt)
{
    vt->problem_setup = NULL;
    vt->patch_setup = NULL;
    vt->patch_physical_bc = NULL;
    vt->patch_single_step_update = NULL;
    vt->run_diagnostics = NULL;
    vt->patch_copy2samesize      = &fclaw2d_patch_copy2samesize;

    vt->patch_average2coarse     = &fclaw2d_patch_average2coarse;
    vt->fort_average2coarse      = &FCLAW2D_FORT_AVERAGE2COARSE;

    vt->patch_interpolate2fine   = &fclaw2d_patch_interpolate2fine;
    vt->fort_interpolate2fine    = &FCLAW2D_FORT_INTERPOLATE2FINE;

    vt->patch_tag4refinement     = &fclaw2d_patch_tag4refinement;
    vt->fort_tag4refinement      = &FCLAW2D_FORT_TAG4REFINEMENT;

    vt->patch_tag4coarsening     = &fclaw2d_patch_tag4coarsening;
    vt->fort_tag4coarsening      = &FCLAW2D_FORT_TAG4COARSENING;

    vt->write_header             = &fclaw2d_output_header_ascii;
    vt->fort_write_header        = &FCLAW2D_FORT_WRITE_HEADER;

    vt->patch_write_file         = &fclaw2d_output_patch_ascii;
    vt->fort_write_file          = &FCLAW2D_FORT_WRITE_FILE;
}

void fclaw2d_set_vtable(fclaw2d_domain_t* domain, fclaw2d_vtable_t *vt)
{
    fclaw2d_domain_attribute_add (domain,"vtable",vt);
}

fclaw2d_vtable_t fclaw2d_get_vtable(fclaw2d_domain_t* domain)
{
    fclaw2d_vtable_t *vt;
    vt = (fclaw2d_vtable_t*) fclaw2d_domain_attribute_access(domain,"vtable",NULL);
    FCLAW_ASSERT(vt != NULL);
    return *vt;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
