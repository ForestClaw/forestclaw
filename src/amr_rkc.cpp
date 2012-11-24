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

// #include "amr_forestclaw.H"
#include "amr_mol.H"


extern "C"
{
    void rkc_stats_(const double& tip1, const double& dt);
}


// This has the header required by the MOL solver, and then it should call
// the specific fclaw_mol_rhs function, which restores the data to the
// tree structure, sets the boundary conditions, iterates over patches,
// to evaluate the rhs hand side for each patch, and returns the rhs vector.
static
void f_rkc(const int& neqn, const double& t_inner, double q[], double rhs[])
{
    // This uses lots of static variables stored in amr_mol.cpp.
    fclaw_mol_rhs(t_inner, q, rhs);
}

// Take a step from t to t+dt
void parabolic_step_rkc(int neqn, double q[], double t, double dt)
{

    double rtol, atol;
    rtol = 1e-4;
    atol = 1e-4;

    int info[4];
    info[0] = 1;
    info[1] = 0;
    info[2] = 0;
    info[3] = 0;

    double tip1 = t + dt;
    int idid = 0;

    double *rkc_work = new double[8+5*neqn];

    // Fortran routine which does time stepping.
    rkc_(neqn, (void*) f_rkc, q, t, tip1, rtol, atol, info, rkc_work, idid);

    if (idid > 2)
    {
        // 1 and 2 are okay
        printf("RKC : ierror > 2;  ierror = %d\n",idid);
        exit(1);
    }
    rkc_stats_(t, tip1);

    delete [] rkc_work;
}
