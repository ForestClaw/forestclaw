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

#include "burgers_user.h"

#include <fclaw2d_include_all.h>

/* Two versions of Clawpack */
#include <fc2d_clawpack46.h>
#include <fc2d_clawpack46_options.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>

#include <clawpack46_user_fort.h>


void burgers_link_solvers(fclaw2d_global_t *glob)
{
    const user_options_t* user = burgers_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt();        

        clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->fort_rpn2_cons = &RPN2CONS_UPDATE;

        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt();
        clawpatch_vt->fort_tag4coarsening = &TAG4COARSENING;
        clawpatch_vt->fort_tag4refinement = &TAG4REFINEMENT;

        fc2d_clawpack46_options_t *claw46_opt = fc2d_clawpack46_get_options(glob);
        if (claw46_opt->order[0] == 3)
        {
            clawpack46_vt->flux2          = &BURGERS_FLUX2;
            clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2BU;
            clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2BU;
            clawpack46_vt->fort_b4step2   = &CLAWPACK46_B4STEP2;
            clawpatch_vt->fort_interpolate_face   = &PERIODIC_FORT_INTERPOLATE_FACE;
            clawpatch_vt->fort_interpolate_corner = &PERIODIC_FORT_INTERPOLATE_CORNER;
            clawpatch_vt->fort_interpolate2fine   = &PERIODIC_FORT_INTERPOLATE2FINE;
        }
        else
        {
            clawpack46_vt->fort_rpn2 = &CLAWPACK46_RPN2;
            clawpack46_vt->fort_rpt2 = &CLAWPACK46_RPT2;
        }
    }
}






