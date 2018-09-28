#include "cudaclaw5_flux2.h"
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */

#include <math.h>

__device__ void rpn2adv_cuda2(int idir, int meqn, int mwaves, int maux,
     double ql[], double qr[], double auxl[], double auxr[],
     double wave[], double s[], double amdq[], double apdq[])
{
    /* wave[mwaves][meqn] */
    /* idir in 0,1 : needed to get correct  */

    wave[0] = qr[0] - ql[0];
    s[0] = auxr[idir];
    amdq[0] = SC_MIN(auxr[idir], 0) * wave[0];
    apdq[0] = SC_MAX(auxr[idir], 0) * wave[0];
}


#define MEQN   10
#define MAUX   20 
#define MWAVES 10

extern "C"
{

__host__ int cudaclaw5_check_dims(int meqn, int maux, int mwaves)
{
    int check;
    check = (meqn <= MEQN) && (maux <= MAUX) && (mwaves <= MWAVES);
    return check;
}
}


__global__ void cudaclaw5_flux2(int idir, int mx, int my, int meqn, int mbc,
                                int maux, double* qold, double* aux, double dx,
                                double dy, double dt, double* cflgrid,
                                double* fm, double* fp, double* gm, double* gp,
                                double* waves, double *speeds,
                                cudaclaw5_cuda_rpn2_t rpn2, void* rpt2,
                                int mwaves) 
{
    /* Static memory seems much faster than dynamic memory */
    double ql[MEQN];
    double qr[MEQN];
    double auxl[MAUX];
    double auxr[MAUX];
    double s[MWAVES];
    double wave[MEQN*MWAVES];
    double amdq[MEQN];
    double apdq[MEQN];

    int mq, mw,m;
    int x_stride_q, y_stride_q, I_q;
    int x_stride_aux, y_stride_aux, I_aux;
    int x_stride_waves, y_stride_waves, I_waves;
    int x_stride_s, y_stride_s, I_speeds;

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;

    if (ix < mx + 2*mbc-1 && iy < my + 2*(mbc-1)) 
    {
        x_stride_q = meqn;
        y_stride_q = (2 * mbc + mx) * x_stride_q;
        I_q = (ix + mbc-1) * x_stride_q + (iy + mbc-1) * y_stride_q;

        x_stride_aux = maux;
        y_stride_aux = (2 * mbc + mx) * x_stride_aux;
        I_aux = (ix + mbc-1) * x_stride_aux + (iy + mbc-1) * y_stride_aux;

        x_stride_waves = mwaves*meqn;
        y_stride_waves = (2 * mbc + mx) * x_stride_waves;
        I_waves = (ix + mbc-1) * x_stride_waves + (iy + mbc-1) * y_stride_waves;

        x_stride_s = mwaves;
        y_stride_s = (2 * mbc + mx) * x_stride_s;
        I_speeds = (ix + mbc-1) * x_stride_s + (iy + mbc-1) * y_stride_s;

        ql[0] = qold[I_q - x_stride_q];
        qr[0] = qold[I_q];
        auxl[0] = aux[I_aux - x_stride_aux];
        auxl[1] = aux[I_aux - x_stride_aux + 1];
        auxr[0] = aux[I_aux];
        auxr[1] = aux[I_aux + 1];

        //rpn2adv_cuda2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);
        rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        for (mq = 0; mq < meqn; mq++) 
        {
            int i = I_q + mq;
            fp[i] = -apdq[mq]; 
            fm[i] = amdq[mq];
        }

        for (m = 0; m < meqn*mwaves; m++)
        {
            int i = I_waves + m;
            waves[i] = wave[m];
        }

        for (mw = 0; mw < mwaves; mw++)
        {
            int i = I_speeds + mw;
            speeds[i] = s[mw];
        }

        /* Limit waves here */

        /* Compute second order corrections */
        double dtdx = dt/dx;
        for(mq = 0; mq < meqn; mq++)
        {
#if 0            
            cqxx(m,i) = cqxx(m,i) + abs_sign
     &             * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
#endif      
            double cqxx = 0;
            for(mw = 0; mw < mwaves; mw++)
            {
                int m = mw*meqn + mq;
                cqxx += fabs(s[mw])*(1. - fabs(s[mw])*dtdx)*wave[m];
            }

        }
    }

    if (ix < mx + 2*(mbc-1) && iy < my + 2*mbc-1) 
    {
        int x_stride = meqn;
        int y_stride = (2 * mbc + mx) * x_stride;
        int I = (ix + mbc-1) * x_stride + (iy + mbc-1) * y_stride;

        int x_stride_aux = maux;
        int y_stride_aux = (2 * mbc + mx) * x_stride_aux;
        int I_aux = (ix + mbc-1) * x_stride_aux + (iy + mbc-1) * y_stride_aux;
        int mq;

        ql[0] = qold[I - y_stride];
        qr[0] = qold[I];
        auxl[0] = aux[I_aux - y_stride_aux];
        auxl[1] = aux[I_aux - y_stride_aux + 1];
        auxr[0] = aux[I_aux];
        auxr[1] = aux[I_aux + 1];

        //rpn2adv_cuda2(1, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);
        rpn2(1, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        for (mq = 0; mq < meqn; mq++) 
        {
            int i = I + mq;
            gp[i] = -apdq[mq]; 
            gm[i] = amdq[mq];
        }
    }
}
