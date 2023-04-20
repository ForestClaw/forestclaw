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


#include "overpressure_user.h"

/* For 2d mappings */
#include "../all/euler3d_user.h"

#include <fclaw2d_map.h>


/* User defined extruded mesh mapping */
static void
overpressure_map_3dx(fclaw2d_map_context_t * cont, int blockno,
               double xc, double yc, double zc,
               double *xp, double *yp, double *zp)
{
    /* Call 2d mapping to get surface.  2d mapping is not scaled in the 
       extruded case. */
    double xp1, yp1, zp1;
    cont->mapc2m(cont,blockno,xc,yc,&xp1,&yp1,&zp1);

    /* In extruded case, no transformations are applied to the 2d mapping */
    int mapping = cont->user_int_3dx[0];

    double maxelev    = cont->user_double_3dx[0];  
    double minz       = cont->user_double_3dx[1];  
    double maxz       = cont->user_double_3dx[2];  
    double midz       = cont->user_double_3dx[3];  
    double scale_bump = cont->user_double_3dx[4];

    if (mapping == 1)
    {
        scale_map(cont,&xp1,&yp1,&zp1);
        shift_map(cont,&xp1,&yp1,&zp1);

        /* Cartesian mesh */
        *xp = xp1;
        *yp = yp1;

        /* Stretch zc into [minz+bump, maxz] */
        double rp2 = xp1*xp1 + yp1*yp1;
        if (minz < midz && midz < maxz)
        {
            double f = (midz-minz)/(maxz-minz);
            if (zc < f)
            {                
                double zlow = minz + scale_bump*exp(-30*rp2);
                *zp = zlow  + (midz-zlow)*(zc/f);
            }
            else
            {
                double zlow = midz;
                *zp = zlow  + (maxz-zlow)*((zc-f)/(1-f));
            }
        }
        else
        {
            /* midz <= minz or midz >= maxz;  assume zc in [0,1] */
            double zlow = minz + scale_bump*exp(-30*rp2);
            *zp = zlow  + (maxz-zlow)*zc;
        }

    }
    else 
    {
        /* Rotate spherical maps */
        rotate_map(cont,&xp1,&yp1,&zp1);

        /* Latlong or sphere meshes */
        double phi = asin(zp1);     

        // returns value in [-pi, pi]
        double theta = atan2(yp1,xp1);      

        /* Get scaling in the radial direction */
        scale_map(cont,&xp1,&yp1,&zp1);
        double rp = sqrt(pow(xp1,2) + pow(yp1,2) + pow(zp1,2));

        /* Extrude in z direction to constant height maxelev. 
           This assumes that zc in [0,1] (default). 
        */
        double R = rp + maxelev*zc;  
        *xp = R*cos(phi)*cos(theta);
        *yp = R*cos(phi)*sin(theta);
        *zp = R*sin(phi);
    }
}


void overpressure_map_extrude(fclaw2d_map_context_t* cont,
                         const double maxelev,
                         const int mapping,
                         const double minz,
                         const double maxz,
                         const double midz,
                         const double scale_bump)

{
    /* May be needed to get more general mappings */
    cont->mapc2m_3dx = overpressure_map_3dx;

    /* Store parameters for use in routine above */
    cont->user_double_3dx[0] = maxelev;
    cont->user_double_3dx[1] = minz;
    cont->user_double_3dx[2] = maxz;
    cont->user_double_3dx[3] = midz;
    cont->user_double_3dx[4] = scale_bump;

    cont->user_int_3dx[0] = mapping;

    /* This is checked in 2d mappings.  If `is_extruded=1`, then the 2d mapping
       will not be scaled, shifted or rotated
    */
    cont->is_extruded = 1;

    return;
}



