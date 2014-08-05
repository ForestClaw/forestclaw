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

#include "amr_includes.H"


void
get_corner_type(fclaw2d_domain_t* domain,
                int icorner,
                fclaw_bool intersects_bdry[],
                fclaw_bool intersects_block[],
                fclaw_bool *interior_corner,
                fclaw_bool *is_block_corner,
                int *block_iface)
{
    // p4est has tons of lookup table like this, can be exposed similarly
    int corner_faces[SpaceDim];
    fclaw2d_domain_corner_faces(domain, icorner, corner_faces);

    // Both faces are at a physical boundary
    fclaw_bool is_phys_corner =
        intersects_bdry[corner_faces[0]] && intersects_bdry[corner_faces[1]];

    // Corner lies in interior of physical boundary edge.
    fclaw_bool corner_on_phys_face = !is_phys_corner &&
             (intersects_bdry[corner_faces[0]] || intersects_bdry[corner_faces[1]]);

    /* Either a corner is at a block boundary (but not a physical boundary),
       or internal to a block.  L-shaped domains are excluded for now
       (i.e. no reentrant corners). */
    *interior_corner = !corner_on_phys_face && !is_phys_corner;

    // Both faces are at a block boundary
    *is_block_corner =
        intersects_block[corner_faces[0]] && intersects_block[corner_faces[1]];

    *block_iface = -1;
    if (!*is_block_corner)
    {
        if (intersects_block[corner_faces[0]])
        {
            // Corner is on a block face.
            *block_iface = corner_faces[0];
        }
        else if (intersects_block[corner_faces[1]])
        {
            // Corner is on a block face.
            *block_iface = corner_faces[1];
        }
    }
}
