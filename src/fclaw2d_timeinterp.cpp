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

#include <fclaw2d_timeinterp.h>
#include <fclaw2d_clawpatch.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


int fclaw2d_timeinterp_has_finegrid_neighbors(fclaw2d_domain_t * domain,
                                              fclaw2d_patch_t *this_patch)
{

    for (int iface = 0; iface < 4; iface++)
    {
        fclaw2d_patch_relation_t nt;
        nt = fclaw2d_patch_get_face_type(this_patch,iface);
        if (nt == FCLAW2D_PATCH_HALFSIZE)
        {
            return 1;
        }
    }

    for (int icorner = 0; icorner < 4; icorner++)
    {
        fclaw2d_patch_relation_t nt;
        int has_corner = !fclaw2d_patch_corner_is_missing(this_patch,icorner);
        if (has_corner)
        {
            nt = fclaw2d_patch_get_corner_type(this_patch,icorner);
            if (nt == FCLAW2D_PATCH_HALFSIZE)
            {
                return 1;
            }
        }
    }
    return 0;


#if 0
    for (int iface = 0; iface < 4; iface++)
    {
        int rproc[2];
        int rblockno;
        int rpatchno[2];
        int rfaceno;

        fclaw2d_patch_relation_t neighbor_type =
            fclaw2d_patch_face_neighbors(domain,
                                         blockno,
                                         patchno,
                                         iface,
                                         rproc,
                                         &rblockno,
                                         rpatchno,
                                         &rfaceno);

        if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            return 1;
        }
    }


    for (int icorner = 0; icorner < 4; icorner++)
    {
        int rproc_corner;
        int cornerpatchno;
        int cornerblockno;
        int rcornerno;
        fclaw2d_patch_relation_t neighbor_type;

        int has_corner_neighbor =
            fclaw2d_patch_corner_neighbors(domain,
                                           blockno,
                                           patchno,
                                           icorner,
                                           &rproc_corner,
                                           &cornerblockno,
                                           &cornerpatchno,
                                           &rcornerno,
                                           &neighbor_type);

        if (has_corner_neighbor & (neighbor_type == FCLAW2D_PATCH_HALFSIZE))
        {
            return 1;
        }
    }
    return 0;
#endif
}



static
void cb_setup_time_interp(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int blockno,
                          int patchno,
                          void *user)
{
#if 0
    /* These ideas are unfortunately really slow! */
    if (fclaw2d_timeinterp_has_finegrid_neighbors(domain,blockno,patchno))
    {
        /* Only create time interpolated for patches with
           fine grid neighbors. */
        double &alpha = *((double*) user);
        fclaw2d_clawpatch_setup_timeinterp(domain,this_patch,alpha);
    }
#endif

    if (fclaw2d_timeinterp_has_finegrid_neighbors(domain,this_patch))
    {
        double &alpha = *((double*) user);
        fclaw2d_clawpatch_setup_timeinterp(domain,this_patch,alpha);
    }

#if 0
    /* Seems better to just fill in each grid without checking to see if this
       grid will be needed for interplation to finer grids */
    double &alpha = *((double*) user);
    fclaw2d_clawpatch_setup_timeinterp(domain,this_patch,alpha);
#endif
}



/* ----------------------------------------------------------------------
   Main routine in this file.  This file assumes that both coarse and
   fine grids have valid interior data;  they only need to exchange (
   via interpolating and averaging) ghost cell values.
   -------------------------------------------------------------------- */

void fclaw2d_timeinterp(fclaw2d_domain_t *domain,
                        int level,double alpha)
{
    /* Store time interpolated data into m_griddata_time_sync. */
    fclaw2d_domain_iterate_level(domain, level,cb_setup_time_interp,
                                     (void *) &alpha);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
