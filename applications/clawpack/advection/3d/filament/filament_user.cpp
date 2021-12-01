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

#include "filament_user.h"


void filament_link_solvers(fclaw2d_global_t *glob)
{
    const user_options_t* user = filament_get_options(glob);
    if (user->claw_version == 4)
    {        
        if (!user->use_claw3d)
        {
            fclaw_global_essentialf("filament_user.cpp: 2d example not implemented in this directory.");
        }
        else
        {
            fc3d_clawpack46_vtable_t *clawpack46_vt = fc3d_clawpack46_vt();

            clawpack46_vt->fort_setprob   = &SETPROB;
            clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;

            /* Used in non-manifold case */
            clawpack46_vt->fort_setaux     = &CLAWPACK46_SETAUX;  
            clawpack46_vt->fort_rpn3       = &CLAWPACK46_RPN3;
            clawpack46_vt->fort_rpt3       = &CLAWPACK46_RPT3;
            clawpack46_vt->fort_rptt3      = &CLAWPACK46_RPTT3;            
        }
    }
    else if (user->claw_version == 5)
    {
        printf("filament_user.cpp : Example not implemented for Claw version 5.\n");
        exit(0);
    }
}