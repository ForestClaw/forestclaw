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
#include <fclaw2d_metric_default_fort.h>
#include <fclaw2d_metric.h>

#include <fclaw2d_global.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_clawpatch.h>

void fclaw2d_metric_average_area(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *fine_patches,
                                 fclaw2d_patch_t *coarse_patch,
                                 int blockno, int coarse_patchno,
                                 int fine0_patchno)

{
    int mx,my, mbc;
    double xlower,ylower,dx,dy;

    double *areacoarse, *areafine;
    int igrid;
    fclaw2d_patch_t *fine_patch;

    fclaw2d_clawpatch_grid_data(glob,coarse_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    areacoarse = fclaw2d_clawpatch_get_area(glob,coarse_patch);

    for(igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        areafine = fclaw2d_clawpatch_get_area(glob,fine_patch);

        FCLAW2D_FORT_AVERAGE_AREA(&mx,&my,&mbc,areacoarse,areafine,&igrid);
    }

    /* Use either exact or approximate method */
    fclaw2d_vt()->metric_area_set_ghost(glob,coarse_patch,blockno,coarse_patchno);
}
