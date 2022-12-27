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

#include "filament_user.h"
#include "filament_fort.h"

static
void filament_patch_setup(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *patch,
                          int blockno,
                          int patchno)
{
    const filament_options_t* user = filament_get_options(glob);
    advection_patch_setup_manifold(glob,patch,blockno,patchno,
                                   user->claw_version); 
}

void filament_link_solvers(fclaw2d_global_t *glob)
{
    const fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

    if (fclaw_opt->manifold)
    {
        fclaw2d_patch_vtable_t *patch_vt = fclaw2d_patch_vt(glob);
        patch_vt->setup = filament_patch_setup;        
    }

    const filament_options_t* user = filament_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);

        clawpack46_vt->fort_setprob   = &FILAMENT_SETPROB;
        clawpack46_vt->fort_qinit     = &FILAMENT_CLAWPACK46_QINIT;

        if (fclaw_opt->manifold)
        {
            clawpack46_vt->fort_rpn2  = &CLAWPACK46_RPN2ADV_MANIFOLD;
            clawpack46_vt->fort_rpt2  = &CLAWPACK46_RPT2ADV_MANIFOLD;
        }
        else
        {
            /* Used in non-manifold case */
            clawpack46_vt->fort_setaux    = &CLAWPACK46_SETAUX;  
            clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
            clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;
        }
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *claw5_vt = fc2d_clawpack5_vt(glob);

        claw5_vt->fort_setprob   = &FILAMENT_SETPROB;
        claw5_vt->fort_qinit     = &FILAMENT_CLAWPACK5_QINIT;
        
        if (fclaw_opt->manifold)
        {
            claw5_vt->fort_rpn2   = CLAWPACK5_RPN2ADV_MANIFOLD;
            claw5_vt->fort_rpt2   = CLAWPACK5_RPT2ADV_MANIFOLD;
        }
        else
        {
            claw5_vt->fort_setaux    = CLAWPACK5_SETAUX;   /* Used in non-manifold case */
            claw5_vt->fort_rpn2      = CLAWPACK5_RPN2ADV;
            claw5_vt->fort_rpt2      = CLAWPACK5_RPT2ADV;
        }
    }
}