#include "../swirl_user.h"

#include <fc2d_cudaclaw.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */



__device__ void swirl_rpn2adv(int idir, int meqn, int mwaves, 
                              int maux, double ql[], double qr[], 
                              double auxl[], double auxr[],
                              double wave[], double s[], 
                              double amdq[], double apdq[])
{

    /* Solve q_t + D q_x = 0, where D = diag([u,u,...,u]), D in R^{meqn x meqn} */
    int mq;

    for(mq = 0; mq < meqn; mq++)
    {
        wave[mq] = qr[mq] - ql[mq];        
    }

    s[0] = auxr[idir];    /* Assume all waves move at the same speed */

    for(mq = 0; mq < meqn; mq++)
    {
        amdq[mq] = SC_MIN(s[0], 0) * wave[mq];
        apdq[mq] = SC_MAX(s[0], 0) * wave[mq];            
    }
}

__device__ cudaclaw_cuda_rpn2_t swirl_rpn2 = swirl_rpn2adv;

void swirl_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, swirl_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (swirl_rpn2adv): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}