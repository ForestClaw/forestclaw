#include "cudaclaw5_flux2.h"
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */

#include <math.h>


/* Use this version (in swirl example) to test performance hit with function
   pointers */
__device__ void rpn2adv(int idir, int meqn, int mwaves, 
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

#if 0
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
#endif


#define MEQN   32
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
                                int maux, int mwaves, double* qold, double* aux, 
                                double dx, double dy, double dt, double* cflgrid,
                                double* fm, double* fp, double* gm, double* gp,
                                double* waves, double *speeds,
                                cudaclaw5_cuda_rpn2_t rpn2, void* rpt2)
{
    int mq, mw, m;
    int x_stride_q, y_stride_q, I_q;
    int x_stride_aux, y_stride_aux, I_aux;
    int x_stride_waves, y_stride_waves, I_waves;
    int x_stride_s, y_stride_s, I_speeds;

    /* Static memory seems much faster than dynamic memory */
    double ql[MEQN];
    double qr[MEQN];
    double auxl[MAUX];
    double auxr[MAUX];
    double s[MWAVES];
    double wave[MEQN*MWAVES];
    double amdq[MEQN];
    double apdq[MEQN];

    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;

    /* Array of structures */
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


    if (idir == 0 && (ix < mx + 2*mbc-1 && iy < my + 2*(mbc-1)))
    {
        for(mq = 0; mq < meqn; mq++)
        {
            ql[mq] = qold[I_q - x_stride_q + mq];
            qr[mq] = qold[I_q + mq];            
        }
        for(m = 0; m < maux; m++)
        {
            auxl[m] = aux[I_aux - x_stride_aux + m];
            auxr[m] = aux[I_aux+m];
        }

        //rpn2adv(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);
        rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        /* Set value at left interface of cell I */
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
            speeds[i] = s[0];
        }        
    }
    else if (idir == 1 && (ix < mx + 2*(mbc-1) && iy < my + 2*mbc-1))
    {
        for(mq = 0; mq < meqn; mq++)
        {
            ql[mq] = qold[I_q - y_stride_q + mq];
            qr[mq] = qold[I_q + mq];            
        }
        for(m = 0; m < maux; m++)
        {
            auxl[m] = aux[I_aux - y_stride_aux + m];
            auxr[m] = aux[I_aux + m];
        }

        //rpn2adv(1, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);
        rpn2(1, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        /* Set value at bottom interface of cell I */
        for (mq = 0; mq < meqn; mq++) 
        {
            int j = I_q + mq;
            gp[j] = -apdq[mq]; 
            gm[j] = amdq[mq];
        }
        /* Assumes array of structures */
        for (m = 0; m < meqn*mwaves; m++)
        {
            int i = I_waves + m;
            waves[i] = wave[m];
        }

        for (mw = 0; mw < mwaves; mw++)
        {
            int i = I_speeds + mw;
            speeds[i] = s[0];
        }
    }
}

__global__ void cudaclaw5_compute_cfl(int idir, int mx, int my, int meqn, int mwaves, 
                                      int mbc, double dx, double dy, double dt, 
                                      double *speeds, double* cflgrid)
{
#if 0    
      # from fortran_source/cudaclaw5_flux2.f */

c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do 50 mw=1,mwaves
         do 50 i=1,mx+1
c          # if s>0 use dtdx1d(i) to compute CFL,
c          # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(mw,i),
     &                          -dtdx1d(i-1)*s(mw,i))
   50       continue
#endif   
    /* Compute largest waves speeds, scaled by dt/dx,  on grid */


}


__device__ void cudaclaw5_second_order(int idir, int mx, int my, int meqn, int mbc,
                                       int maux, double* qold, double* aux, double dx,
                                       double dy, double dt, double* cflgrid,
                                       double* fm, double* fp, double* gm, double* gp,
                                       double* waves, double *speeds,
                                       cudaclaw5_cuda_rpn2_t rpn2, void* rpt2,
                                       int mwaves) 
{    
    int mq, mw, m;

    /* TODO : Limit waves here */


    /* Compute second order corrections */
    double dtdx = dt/dx;
    for(mq = 0; mq < meqn; mq++)
    {
        double cqxx = 0;
        for(mw = 0; mw < mwaves; mw++)
        {
            m = mw*meqn + mq;
            cqxx += fabs(speeds[mw])*(1.0 - fabs(speeds[mw])*dtdx)*waves[m];
        }
    }
}

