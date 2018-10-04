#include "../swirl_user.h"

#include <fc2d_cudaclaw.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */

__device__ double psi(double x, double y)
{
    return (pow(sin(M_PI*x), 2) * pow(sin(M_PI*y),2)) / M_PI;
} 

__device__ void swirl_b4step2_test(int mbc, int mx, int my, int meqn, double q[],
                              double xlower, double ylower, double dx, double dy, 
                              double time, double dt, int maux, 
                              double aux[], int ipatch, int jpatch, double tperiod)
{
    double vt;
    double xll, yll;

    xll = xlower + (ipatch-1)*dx;
    yll = ylower + (jpatch-1)*dy;

    vt = cos(2*M_PI*(time+dt/2.0)/tperiod);
    aux[0] = (psi(xll,yll+dy) - psi(xll,yll)) / dy;
    aux[1] = - (psi(xll+dx,yll) - psi(xll,yll)) / dx;

    aux[0] = vt * aux[0];
    aux[1] = vt * aux[1];
}

__device__ cudaclaw_cuda_b4step2_t swirl_b4step2 = swirl_b4step2_test;

void swirl_assign_b4step2(cudaclaw_cuda_b4step2_t *b4step2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(b4step2, swirl_b4step2, sizeof(cudaclaw_cuda_b4step2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (swirl_b4step2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}