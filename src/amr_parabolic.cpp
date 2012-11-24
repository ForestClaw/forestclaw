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
#include "amr_parabolic.H"

void parabolic_step_rkc(fclaw_mol_rhs_t f_exp, int neqn,double q[], double dx,
                        double t,double dt)
{
    extern f_rhs;
    double rtol = 1e-4;
    double atol = 1e-4;
    int info[4];
    info[0] = 1;
    info[1] = 0;
    info[2] = 0;
    info[3] = 0;

    double tip1 = t_level + dt;
    int idid = 0;

    double *rkc_work = new double[8+5*neqn];

    rkc_(neqn, f_rhs, q, t_level, tip1, rtol, atol, info, rkc_work, idid);

    if (idid > 2)
    {
        // 1 and 2 are okay
        printf("RKC : ierror > 2;  ierror = %d\n",idid);
        exit(1);
    }
    printf("RKC : stages used : %d\n",get_stages_rkc_());
    delete [] rkc_work;
}


void parabolic_step_feuler(fclaw_mol_rhs_t f_exp,int neqn,double q[], double dx,
                           double t,double dt)
{
    double *feuler_work = new double[neqn];

    // Take step from ti to tip1 using forward Euler
    // Storage for RHS is supplied.  The vector q gets updated.
    feuler_(f_exp,neqn,q,dx,t_level,dt,feuler_work);
    delete [] feuler_work;
}

void single_grid_rhs(ClawPatch *cp, double t, double *rhs)
{
    double dx = cp->dx();
    double dy = cp->dy();
    int mbc = cp->mbc();
    int meqn = cp->meqn();
    int mx = cp->mx();
    int my = cp->my();

    int size = (mx+2*mbc)*(my+2*mbc)*meqn;

    double *diff = new double[size];
    double *reac = new double[size];
    double *q = cp->current_data_ptr();

    apply_lb_(mx,my,mbc,meqn,dx,dy,q,diff);
    src_diffusion_(mx,my,mbc,meqn,t,diff);
    src_reaction_(mx,my,mbc,meqn,t,q,reac);
    compute_rhs_exp_(mx,my,mbc,meqn,diff,reac,rhs);

    set_zeros_(mx,my,mbc,meqn,rhs);
}
