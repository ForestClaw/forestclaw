#include "../swirl_user.h"

#include <fc2d_cudaclaw.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */

__managed__ double s_tperiod;

__device__ double psi(double x, double y)
{
    return (pow(sin(M_PI*x),2) * pow(sin(M_PI*y),2)) / M_PI;
} 
      
void swirl_setprob(double period_in)
{
    s_tperiod = period_in;
}


__device__ void swirl_b4step2_test(int mbc, int mx, int my, int meqn, double q[],
                                   double xlower, double ylower, double dx, double dy, 
                                   double time, double dt, int maux, 
                                   double aux[], int i, int j)
{
    double vt;
    double xll, yll;
    double p1,p2,p3;    

    xll = xlower + (i-1)*dx;
    yll = ylower + (j-1)*dy;
    vt = cos(2*M_PI*(time+dt/2.0)/s_tperiod);

    p1 = psi(xll,yll+dy);
    p2 = psi(xll,yll);
    p3 = psi(xll+dx,yll);

    aux[0] = (p1-p2) / dy;
    aux[1] = - (p3-p2) / dx;

    aux[0] *= vt;
    aux[1] *= vt;
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
