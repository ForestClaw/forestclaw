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
#include "clawpack_fort.H"

// This is called if you want to only compute the right hand side for the
// single step routine.
double waveprop_rhs(fclaw2d_domain_t *domain,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    double t,
                    double *rhs)
{
    // This should evaluate the right hand side, but not actually do the update.
    // This will be useful in cases where we want to use something other than
    // a single step method.  For example, in a RK scheme, one might want to
    // call the right hand side to evaluate stages.
    return 0;
}

// This is called from the single_step callback.
// and is of type 'flaw_single_step_t'
double waveprop_update(fclaw2d_domain_t *domain,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    const amr_options_t* gparms = get_domain_parms(domain);
    int level = this_patch->level;

    ClawPatch *cp = get_clawpatch(this_patch);

    set_block_(&this_block_idx);

    double* qold = cp->current_data_ptr();
    double* aux = cp->aux_data_ptr();

    cp->save_current_step();  // Save for time interpolation

    // Global to all patches
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int maux = gparms->maux;
    int mwaves = gparms->mwaves;

    // Specific to the patch
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();


    // We also call a 'b4step2' in clawpatch2, below.  But it won't
    // do anything in the mapped case.
    if (gparms->manifold)
    {
        double *xp = cp->xp();
        double *yp = cp->yp();
        double *zp = cp->zp();
        double *xd = cp->xd();
        double *yd = cp->yd();
        double *zd = cp->zd();

        b4step2_mapped_(mx,my, mbc,meqn,qold,dx,dy,xp,yp,zp,xd,yd,zd,
                        t, dt, maux, aux);
    }

    int maxm = max(mx,my);

    double cflgrid;

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    clawpatch2_(maxm, meqn, maux, mbc, gparms->method,
                gparms->mthlim, gparms->mcapa, mwaves, mx, my, qold,
                aux, dx, dy, dt, cflgrid, work, mwork, xlower, ylower,level,
                t, fp, fm, gp, gm);

    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}
