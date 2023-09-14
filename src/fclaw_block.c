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

#include <fclaw_block.h>
#include <fclaw_global.h>
#include <fclaw2d_defs.h>
#include <fclaw3d_defs.h>

#include <fclaw_patch.h>
#include <fclaw2d_patch.h>
#include <fclaw3d_patch.h>
#include <forestclaw2d.h>
#include <forestclaw3d.h>


void fclaw_block_get_block_boundary(fclaw_global_t * glob,
                                      fclaw_patch_t * patch,
                                      int *intersects_block)
{
    if(glob->domain->refine_dim == 2)
    {
        for (int iside = 0; iside < FCLAW2D_NUMFACES; iside++)
        {
            int iface_flags = fclaw2d_patch_block_face_flags[iside];
            fclaw2d_patch_t* patch_2d = fclaw_patch_get_2d_patch(patch);
            int is_block_face = (patch_2d->flags & iface_flags) != 0;

            /* True for physical and block boundaries across a face */
            intersects_block[iside] = is_block_face;
        }
    }
    else
    {
        for (int iside = 0; iside < FCLAW3D_NUMFACES; iside++)
        {
            int iface_flags = fclaw3d_patch_block_face_flags[iside];
            fclaw3d_patch_t* patch_3d = fclaw_patch_get_3d_patch(patch);
            int is_block_face = (patch_3d->flags & iface_flags) != 0;

            /* True for physical and block boundaries across a face */
            intersects_block[iside] = is_block_face;
        }
    }
}
