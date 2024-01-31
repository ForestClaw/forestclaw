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


#include "latlong_user.h"

/* For 2d mappings */
#include "../all/advection_user.h"

#include <fclaw_map.h>


/* User defined extruded mesh mapping */
static void
latlong_map_3dx(fclaw_map_context_t * cont, int blockno,
               double xc, double yc, double zc,
               double *xp, double *yp, double *zp)
{
    /* Call 2d mapping to get surface.  2d mapping is not scaled in the 
       extruded case. */
    double xp1, yp1, zp1;
    cont->mapc2m(cont,blockno,xc,yc,&xp1,&yp1,&zp1);

    // returns value in [-pi/2, pi/2]
    double phi = asin(zp1);     

    // returns value in [-pi, pi]
    double theta = atan2(yp1,xp1);      
    if (theta < 0)
        theta += 2*M_PI;

    /* Get scaling in the radial direction */
    scale_map(cont,&xp1,&yp1,&zp1);
    double rp = sqrt(pow(xp1,2) + pow(yp1,2) + pow(zp1,2));

    /* Extrude in z direction to constant height maxelev. */
    double maxelev = cont->user_double_3d[0];  
    double R = rp + maxelev*zc;  
    *xp = R*cos(phi)*cos(theta);
    *yp = R*cos(phi)*sin(theta);
    *zp = R*sin(phi);
}


void latlong_map_extrude(fclaw_map_context_t* cont,
                         const double maxelev)

{
    /* May be needed to get more general mappings */
    cont->mapc2m_3d = latlong_map_3dx;

    cont->user_double_3d[0] = maxelev;

    cont->is_extruded = 1;

    return;
}



