#include "cudaclaw5_flux2.h"
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
__global__ void cudaclaw5_flux2(int idir, int mx, int my, int meqn, int mbc,
                                int maux, double* qold, double* aux, double dx,
                                double dy, double dt, double* cflgrid,
                                double* fm, double* fp, double* gm, double* gp,
                                cudaclaw5_cuda_rpn2_t rpn2, void* rpt2,
                                int mwaves) 
{
    double* ql = new double[meqn];
    double* qr = new double[meqn];
    double* auxl = new double[maux];
    double* auxr = new double[maux];
    double* s = new double[mwaves];
    double* wave = new double[meqn * mwaves];
    double* amdq = new double[meqn];
    double* apdq = new double[meqn];

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;

    if (ix < mx + 2*mbc-1 && iy < my + 2*(mbc-1)) {
        int x_stride = meqn;
        int y_stride = (2 * mbc + mx) * x_stride;
        int I = (ix + mbc - 1) * x_stride + (iy + mbc - 1) * y_stride;
        int x_stride_aux = maux;
        int y_stride_aux = (2 * mbc + mx + 1) * x_stride_aux;
        int I_aux = (ix + mbc - 1) * x_stride_aux + (iy + mbc - 1) * y_stride_aux;
        int mq;

        ql[0] = qold[I - x_stride];
        qr[0] = qold[I];
        auxl[0] = aux[I_aux - x_stride_aux];
        auxl[1] = aux[I_aux - x_stride_aux + 1];
        auxr[0] = aux[I_aux];
        auxr[1] = aux[I_aux + 1];

        rpn2adv_cuda2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        for (mq = 0; mq < meqn; mq++) {
            int i = I + mq;
            qold[i] = qold[i] - dt / dx * (amdq[mq] + apdq[mq]);
            //qold[i] = 1;
        }
    }
    delete[] ql;
    delete[] qr;
    delete[] auxl;
    delete[] auxr;
    delete[] s;
    delete[] wave;
    delete[] amdq;
    delete[] apdq;
}
