#include "../radial_user.h"

#include <fc2d_cudaclaw_check.cu>
    
__constant__ double s_rho;
__constant__ double s_bulk;
__constant__ double s_c;
__constant__ double s_z;

void radial_setprob_cuda(double rho, double bulk)
{
    double c,z;
    c = sqrt(bulk/rho);
    z = c*rho;
    CHECK(cudaMemcpyToSymbol(s_rho,  &rho, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_bulk, &rho, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_c,    &c,   sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_z,    &z,   sizeof(double)));
}


__device__ void radial_rpn2acoustics(int idir, int meqn, int mwaves, 
                                     int maux, double ql[], double qr[], 
                                     double auxl[], double auxr[],
                                     double wave[], double s[], 
                                     double amdq[], double apdq[])
{
    double alpha1, alpha2, delta[2];
    int mu,mv,mq;

    mu = 1+idir;
    mv = 2-idir;    

    s[0] = -s_c;
    s[1] =  s_c;

    delta[0] = qr[0] - ql[0];
    delta[1] = qr[mu] - ql[mu];
    alpha1 = ( -delta[0] + s_z*delta[1]) / (2.0*s_z);
    alpha2 = (  delta[0] + s_z*delta[1]) / (2.0*s_z);

    /* left going wave */
    wave[0]  = -alpha1*s_z;
    wave[mu] = alpha1;
    wave[mv] = 0;

    /* Right going wave */
    wave[3]       = alpha2*s_z;
    wave[meqn+mu] = alpha2;
    wave[meqn+mv] = 0;

    for (mq = 0; mq < meqn; mq++)
    {
        amdq[mq] = s[0]*wave[mq];
        apdq[mq] = s[1]*wave[meqn + mq];
    }
}

__device__ cudaclaw_cuda_rpn2_t radial_rpn2 = radial_rpn2acoustics;

void radial_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, radial_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radial_rpn2acoustics): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void radial_rpt2acoustics(int idir, int meqn, int mwaves, int maux,
                                     double ql[], double qr[], 
                                     double aux1[], double aux2[], double aux3[],
                                     int imp, int pm, double asdq[],
                                     double bmasdq[], double bpasdq[])
{

    double alpha1, alpha2, delta[2];
    int mu,mv;

    mu = 1+idir;
    mv = 2-idir;    

    delta[0] = asdq[0];
    delta[1] = asdq[mv];
    alpha1 = ( -delta[0] + s_z*delta[1]) / (2.0*s_z);
    alpha2 = (  delta[0] + s_z*delta[1]) / (2.0*s_z);

    /* Down going wave */
    bmasdq[0]  = s_c * alpha1 * s_z;
    bmasdq[mu] = 0;
    bmasdq[mv] = -s_c * alpha1;

    /* Up going wave */
    bpasdq[0]       = s_c * alpha2 * s_z;
    bpasdq[mu] = 0;
    bpasdq[mv] = s_c * alpha2;
}

__device__ cudaclaw_cuda_rpt2_t radial_rpt2 = radial_rpt2acoustics;

void radial_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, radial_rpt2, sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radial_rpt2acoustics): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}

