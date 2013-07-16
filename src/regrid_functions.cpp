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

#include "amr_utils.H"
#include "amr_regrid.H"


/* -----------------------------------------------------------------
   User defined routines for linking in customized tagging and
   interpolation/averaging routines
   ----------------------------------------------------------------- */
    void link_regrid_functions(fclaw2d_domain_t* domain,
                               fclaw2d_patch_tag4refinement_t patch_tag4refinement,
                               fclaw2d_patch_tag4coarsening_t patch_tag4coarsening)
{
    fclaw2d_regrid_functions_t *rf = get_regrid_functions(domain);

    rf->f_patch_tag4refinement = patch_tag4refinement;
    rf->f_patch_tag4coarsening = patch_tag4coarsening;
    rf->f_patch_interpolate2fine = &patch_interpolate2fine;
    rf->f_patch_average2coarse = &patch_average2coarse;
    rf->f_patch_copy2samesize = &patch_copy2samesize;
};

fclaw_bool patch_tag4coarsening_dummy(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno,
                                      int patchno)
{
    return fclaw_true;  /* don't coarsen */
}

fclaw_bool patch_tag4refinement_dummy(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int this_block_idx, int this_patch_idx,
                                      int initflag)
{
    return fclaw_false;  /* don't coarsen */
}


void patch_copy2samesize(fclaw2d_domain_t* domain, fclaw2d_patch_t *old_patch,
                         fclaw2d_patch_t* new_patch, int blockno, int old_patchno,
                         int new_patchno)
{
    ClawPatch *cp_old = get_clawpatch(old_patch);
    ClawPatch *cp_new = get_clawpatch(new_patch);

    cp_new->copyFrom(cp_old);
}


void patch_interpolate2fine(fclaw2d_domain_t* domain, fclaw2d_patch_t *coarse_patch,
                            fclaw2d_patch_t* fine_patch,
                            int this_blockno, int coarse_patchno,
                            int fine_patchno, int igrid)

{
    const amr_options_t* gparms = get_domain_parms(domain);

    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int refratio = gparms->refratio;

    ClawPatch *cp_coarse = get_clawpatch(coarse_patch);
    double *qcoarse = cp_coarse->q();

    ClawPatch *cp_fine = get_clawpatch(fine_patch);
    double *qfine = cp_fine->q();

    // Use linear interpolation with limiters.
    interpolate_to_fine_patch_(mx,my,mbc,meqn,qcoarse,qfine,
                               p4est_refineFactor, refratio,igrid);
    if (gparms->manifold)
    {
        double *areacoarse = cp_coarse->area();
        double *areafine = cp_fine->area();

        fixcapaq2_(mx, my, mbc, meqn, qcoarse, qfine, areacoarse, areafine,
                   p4est_refineFactor, refratio, igrid);
    }
}


void patch_average2coarse(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *fine_siblings,
                          fclaw2d_patch_t *coarse_patch,
                          int blockno, int fine_patchno,
                          int coarse_patchno)

{
    const amr_options_t* gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int refratio = gparms->refratio;

    ClawPatch *cp_coarse = get_clawpatch(coarse_patch);
    double *qcoarse = cp_coarse->q();

    for(int igrid = 0; igrid < NumSiblings; igrid++)
    {
        ClawPatch *cp_fine = get_clawpatch(&fine_siblings[igrid]);
        double *qfine = cp_fine->q();
        if (gparms->manifold)
        {
            double *areacoarse = cp_coarse->area();
            double *areafine = cp_fine->area();
            average_to_coarse_mapped_(mx,my,mbc,meqn,qcoarse,qfine,
                                      areacoarse, areafine,
                                      p4est_refineFactor,
                                      refratio, igrid);
        }
        else
        {
            average_to_coarse_patch_(mx,my,mbc,meqn,qcoarse,qfine,
                                     p4est_refineFactor,refratio,igrid);
        }
    }
}

void initialize_regrid_functions(fclaw2d_regrid_functions_t* regrid_functions)
{
    fclaw2d_regrid_functions_t *rf = regrid_functions;

    rf->f_patch_tag4refinement = &patch_tag4refinement_dummy;
    rf->f_patch_tag4coarsening = &patch_tag4coarsening_dummy;
    rf->f_patch_interpolate2fine = &patch_interpolate2fine;
    rf->f_patch_average2coarse = &patch_average2coarse;
    rf->f_patch_copy2samesize = &patch_copy2samesize;

}

void copy_regrid_functions(fclaw2d_regrid_functions_t* old_regrid_functions,
                           fclaw2d_regrid_functions_t* new_regrid_functions)
{
    fclaw2d_regrid_functions_t* oldrf = old_regrid_functions;
    fclaw2d_regrid_functions_t* newrf = new_regrid_functions;

    newrf->f_patch_tag4refinement = oldrf->f_patch_tag4refinement;
    newrf->f_patch_tag4coarsening = oldrf->f_patch_tag4coarsening;
    newrf->f_patch_interpolate2fine = oldrf->f_patch_interpolate2fine;
    newrf->f_patch_average2coarse = oldrf->f_patch_average2coarse;
    newrf->f_patch_copy2samesize = oldrf->f_patch_copy2samesize;
}

fclaw2d_regrid_functions_t* get_regrid_functions(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    return ddata->regrid_functions;
}
