#include "../radial_user.h"

//#include <fc2d_cudaclaw5.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */
//#include <cassert>

typedef double real;
static __device__ real claw_zero = 0.0;
static __device__ real rho = 1.0;
static __device__ real bulk = 4.0;

__device__ void radial_rpn2acoustics(int idir, int meqn, int mwaves, 
                              int maux, double ql[], double qr[], 
                              double auxl[], double auxr[],
                              double wave[], double s[], 
                              double amdq[], double apdq[])
{
    /* wave[mwaves][meqn] */
    /* idir in 0,1 : needed to get correct  */

    // TODO: this should be replaced with acoustics riemann solver

    //FCLAW_ASSERT(mwaves == 2);
    //FCLAW_ASSERT(meqn == 3);

    // TODO: pass in bulk and rho
    real c = sqrt(bulk/rho);
    real z = c*rho;

    // if we use template, we don't have to have this branching here
    if (0 == idir) // x-direction
    {
        s[0] = -c;
        s[1] =  c;

        real alpha1 = ( ql[0] - qr[0] + z*(qr[1] - ql[1])) / (2*z);
        real alpha2 = ( qr[0] - ql[0] + z*(qr[1] - ql[1])) / (2*z);

        // TODO: might want to replace double[] wave in argument list with
        // double* wave
        wave[0] = -alpha1*z;
        wave[1] = alpha1;
        wave[2] = claw_zero;

        wave[3] = alpha2*z;
        wave[4] = alpha2;
        wave[5] = claw_zero;
    }
    else if (1 == idir) // y-direction
    {
        s[0] = -c;
        s[1] =  c;

        real alpha1 = ( ql[0] - qr[0] + z*(qr[2] - ql[2])) / (2*z);
        real alpha2 = ( qr[0] - ql[0] + z*(qr[2] - ql[2])) / (2*z);

        wave[0] = -alpha1*z;
        wave[1] = claw_zero;
        wave[2] = alpha1;

        wave[3] = alpha2*z;
        wave[4] = claw_zero;
        wave[5] = alpha2;
    }
    else printf("Invalid value for idir in riemann solver\n");

    for (int mq = 0; mq < meqn; mq++)
    {
        amdq[mq] = s[0]*wave[mq];
        apdq[mq] = s[1]*wave[meqn+mq];
    }


}

__device__ cudaclaw5_cuda_rpn2_t radial_rpn2 = radial_rpn2acoustics;

void radial_assign_rpn2(cudaclaw5_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, radial_rpn2, sizeof(cudaclaw5_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radial_rpn2adv): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}
