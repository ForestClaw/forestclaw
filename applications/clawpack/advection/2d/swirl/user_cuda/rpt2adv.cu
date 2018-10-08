#include "../swirl_user.h"

#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */


__device__ void swirl_rpt2adv(int idir, int meqn, int mwaves, int maux,
                              double ql[], double qr[], 
                              double aux1[], double aux2[], double aux3[],
                              int imp, double asdq[],
                              double bmasdq[], double bpasdq[])
{
    int mq, kv;

    kv = 1-idir;
    for(mq = 0; mq < meqn; mq++)
    {
        bmasdq[mq] = SC_MIN(aux2[imp*maux + kv], 0) * asdq[mq];
        bpasdq[mq] = SC_MAX(aux3[imp*maux + kv], 0) * asdq[mq];            
    }
}

__device__ cudaclaw_cuda_rpt2_t swirl_rpt2 = swirl_rpt2adv;

void swirl_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, swirl_rpt2, sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (swirl_rpt2adv): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}