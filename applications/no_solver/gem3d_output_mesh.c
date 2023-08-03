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

#include "gem3d_output_mesh.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>


void cb_gem3d_output_mesh (fclaw2d_domain_t * domain,
                           fclaw2d_patch_t * this_patch,
                           int this_block_idx, int this_patch_idx,
                           void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;

    int patch_num;
    int level;

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);

    FILE *f;
    f = (FILE*) s->user;

    /* Write out 8 possible neighbors for each mesh */
    /* <n>  <level>  n00 n01  n20 n21  n10  n11  n30 n31  */
    /* 63     3   62    62   -1   -1   61   61   -1  -1 */
    // fprintf(f,"%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n",patch_num,level, ...)

    fprintf(f,"%8d%6d",patch_num,level);

    int iface;
    for(iface = 0; iface < 4; iface++)
    {
        int rproc[2];
        int rblockno;
        int rpatchno[2];
        int rfaceno;
        fclaw2d_patch_relation_t neighbor_type =
                      fclaw2d_patch_face_neighbors(domain,
                                                   this_block_idx,
                                                   this_patch_idx,
                                                   iface,
                                                   rproc,
                                                   &rblockno,
                                                   rpatchno,
                                                   &rfaceno);
        int ir[2];
        if (neighbor_type == FCLAW2D_PATCH_BOUNDARY)
        {
            ir[0] = -1;
            ir[1] = -1;  
        }
        else if (neighbor_type == FCLAW2D_PATCH_SAMESIZE)
        {
            ir[0] = rpatchno[0]; 
            ir[1] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_PATCH_DOUBLESIZE)
        {
            ir[0] = rpatchno[0];
            ir[1] = rpatchno[0];
        }
        else if (neighbor_type == FCLAW2D_PATCH_HALFSIZE)
        {
            ir[0] = rpatchno[0];
            ir[1] = rpatchno[1];
        }
        else
        {
            fclaw_global_essentialf("Invalid neighbor type.\n");
            exit(0);
        }

        fprintf(f,"%8d %8d",ir[0], ir[1]);
    }
    fprintf(f,"\n");
}


/* This function isn't virtualized;  should it be? */
static
void gem3d_output_header_mesh(fclaw2d_global_t* glob, int iframe)
{
    char matname1[11];  /* One extra for the null-termination character */
 
    sprintf(matname1,"mesh.h%04d",iframe); /* number of grids */

    int ngrids = glob->domain->global_num_patches;

    FILE *f;
    f = fopen(matname1,"w");
    fprintf(f,"%d",ngrids);
    fclose(f);
}

    /*--------------------------------------------------------------------
    Public interface
    Use this function as follows : 
           fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
           vt->output_frame = &fclaw2d_clawpatch_output_ascii;
    -------------------------------------------------------------------- */

void gem3d_output_mesh(fclaw2d_global_t* glob,int iframe)
{

    fclaw2d_domain_t *domain = glob->domain;

    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        gem3d_output_header_mesh(glob,iframe);
    }

    char matname2[11];
    sprintf(matname2,"mesh.m%04d",iframe); /* mesh connectivity */    

    FILE *f;
    f = fopen(matname2,"w");
    fclaw2d_global_iterate_patches (glob, cb_gem3d_output_mesh, f);
    fclose(f);

    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */


}

