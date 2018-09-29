#include "../radial_user.h"

#include <fc2d_cudaclaw5.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */



__device__ void radial_rpn2acoustics(int idir, int meqn, int mwaves, 
                              int maux, double ql[], double qr[], 
                              double auxl[], double auxr[],
                              double wave[], double s[], 
                              double amdq[], double apdq[])
{
    /* wave[mwaves][meqn] */
    /* idir in 0,1 : needed to get correct  */

    // TODO: this should be replaced with acoustics riemann solver
    wave[0] = qr[0] - ql[0];
    s[0] = auxr[idir];
    amdq[0] = SC_MIN(auxr[idir], 0) * wave[0];
    apdq[0] = SC_MAX(auxr[idir], 0) * wave[0];
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
