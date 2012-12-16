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
#include "amr_mol.H"

extern "C"
{
    void apply_lb_(const int& mx, const int& my,const int& mbc,
                   const int& meqn,const double &dx, const double& dy,
                   double q[], double diff[]);

    void src_diffusion_(const int& mx,const int& my,const int& mbc,const int& meqn,
                        const double& t, double diff[]);

    void src_reaction_(const int& mx,const int& my,const int& mbc,const int& meqn,
                       const double& t, double q[], double reac[]);

    void compute_rhs_exp_(const int& mx,const int& my,const int& mbc,
                          const int& meqn,double diff[],
                          double reac[], double rhs[]);

    void set_zeros_(const int& mx,const int& my, const int& mbc,
                    const int& meqn, double rhs[]);

}


void allencahn_rhs(fclaw2d_domain_t *domain,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   double t,
                   double *rhs)
{
    ClawPatch *cp = get_clawpatch(this_patch);
    const amr_options_t *gparms = get_domain_parms(domain);

    double dx = cp->dx();
    double dy = cp->dy();
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;
    int mx = gparms->mx;
    int my = gparms->my;

    int size = (mx+2*mbc)*(my+2*mbc)*meqn;

    double *diff = new double[size];
    double *reac = new double[size];
    double *q = cp->current_data_ptr();

    apply_lb_(mx,my,mbc,meqn,dx,dy,q,diff);
    src_diffusion_(mx,my,mbc,meqn,t,diff);
    src_reaction_(mx,my,mbc,meqn,t,q,reac);
    compute_rhs_exp_(mx,my,mbc,meqn,diff,reac,rhs);

    set_zeros_(mx,my,mbc,meqn,rhs);

    delete [] diff;
    delete [] reac;
}
