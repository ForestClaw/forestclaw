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

/* User defined extruded mesh mapping */
static void
filament_map_3dx(fclaw2d_map_context_t * cont, int blockno,
               double xc, double yc, double zc,
               double *xp, double *yp, double *zp)
{
    /* Call 2d mapping to get surface.  2d mapping is not scaled in the 
       extruded case. */
    cont->mapc2m(cont,blockno,xc,yc,xp,yp,zp);

    /* Extrude in z direction to constant height maxelev. */
    double maxelev = cont->user_double_3dx[0]; 
    *zp = maxelev*zc;

    /* Whether it makes sense to scale/shift this mapping is up to the user */
    scale_map(cont,xp,yp,zp);
    shift_map(cont,xp,yp,zp);
}


void filament_map_extrude(fclaw2d_map_context_t* cont,
                          const double maxelev)

{
    /* May be needed to get more general mappings */
    cont->mapc2m_3dx = filament_map_3dx;

    cont->user_double_3dx[0] = maxelev;

    cont->is_extruded = 1;

    return;
}



