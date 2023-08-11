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

#include "amr_forestclaw.H"
#include "ClawPatch.H"
#include "circles_user.H"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void sphere_link_patch(fclaw_domain_t *domain)
{
}

void sphere_setprob(fclaw_domain_t* domain)
{
    setprob_();
}



/* -----------------------------------------------------------------
   Default routine for tagging patches for refinement and coarsening
   ----------------------------------------------------------------- */
int sphere_patch_tag4refinement(fclaw_domain_t *domain,
                                       fclaw_patch_t *this_patch,
                                       int blockno, int this_patch_idx,
                                       int initflag)
{
    /* ----------------------------------------------------------- */
    // Global parameters
    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* xp = cp->xp();
    double* yp = cp->yp();

    int tag_patch = 0;  // == 0 or 1
    tag4refinement_(mx,my,mbc,meqn,xlower,ylower,dx,dy,initflag,
                    blockno, tag_patch,xp,yp);
    return tag_patch == 1;
}

int sphere_patch_tag4coarsening(fclaw_domain_t *domain,
                                       fclaw_patch_t *this_patch,
                                       int blockno,
                                       int patchno)
{
    /* ----------------------------------------------------------- */
    // Global parameters
    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    int tag_patch = 1;  // == 0 or 1
    tag4coarsening_(mx,my,mbc,meqn,xlower,ylower,dx,dy,q,tag_patch);
    return tag_patch == 0;
}

void metric_write_header(fclaw_domain_t* domain, int iframe,int ngrids)
{
    const fclaw_options_t *gparms = get_domain_parms(domain);
    double time = get_domain_time(domain);

    printf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);

    // Write out header file containing global information for 'iframe'
    int meqn = gparms->meqn + 1;  // Write out an extra field
    int maux = 0;
    write_tfile_(iframe,time,meqn,ngrids,maux);

    // This opens file 'fort.qXXXX' for replace (where XXXX = <zero padding><iframe>, e.g. 0001,
    // 0010, 0114), and closes the file.
    new_qfile_(iframe);
}


void metric_write_output(fclaw_domain_t *domain, fclaw_patch_t *this_patch,
                         int this_block_idx, int this_patch_idx,
                         int iframe,int num,int level)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const fclaw_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    // Other input arguments
    int maxmx = mx;
    int maxmy = my;

    double *xp = cp->xp();
    double *yp = cp->yp();

    // This opens a file for append.  Now, the style is in the 'clawout' style.
    int matlab_level = level + 1;

    /* ------------------------------------------------------------- */
    // This opens a file for append.  Now, the style is in the 'clawout' style.
    output_metric_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                   iframe,num,matlab_level,this_block_idx,xp,yp);
}




#ifdef __cplusplus
#if 0
{
#endif
}
#endif
