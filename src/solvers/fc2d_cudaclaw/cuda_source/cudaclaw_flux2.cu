#include "cudaclaw_flux2.h"
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */
#include "cudaclaw_allocate.h"

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

#define MEQN   32
#define MAUX   20 
#define MWAVES 10

extern "C"
{

__host__ int cudaclaw_check_dims(int meqn, int maux, int mwaves)
{
    int check;
    check = (meqn <= MEQN) && (maux <= MAUX) && (mwaves <= MWAVES);
    return check;
}
}


__global__
void cudaclaw_flux2_and_update_batch (int mx, int my, int meqn, int mbc, 
                                int maux, int mwaves, double dt,
                                cudaclaw_fluxes_t* array_fluxes_struct_dev,
                                cudaclaw_cuda_rpn2_t rpn2)
{
    // TODO: check this device function does not depend on blockIdx.z inside
    cudaclaw_flux2_and_update(mx,my,meqn,mbc,maux,mwaves,
            dt/array_fluxes_struct_dev[blockIdx.z].dx,
            dt/array_fluxes_struct_dev[blockIdx.z].dy,
                                    &(array_fluxes_struct_dev[blockIdx.z]),
                                    rpn2);
}

__device__  
void cudaclaw_flux2_and_update (int mx, int my, int meqn, int mbc, 
                           int maux, int mwaves, double dtdx, double dtdy,
                           cudaclaw_fluxes_t* fluxes,
                           cudaclaw_cuda_rpn2_t rpn2)
{
    //TODO: get this from Scott
}


__global__ void cudaclaw_flux2(int mx, int my, int meqn, int mbc,
                                int maux, int mwaves, 
                                double dtdx, double dtdy,
                                double* qold, double* aux, 
                                double* fm, double* fp, double* gm, double* gp,
                                double* waves, double *speeds,
                                cudaclaw_cuda_rpn2_t rpn2)
{
    int mq, mw, m;
    int xs, ys, zs;
    int I, I_q, I_aux, I_waves, I_speeds;

    /* Static memory seems much faster than dynamic memory */
    extern __shared__ double shared_mem[];
    double* ql   = shared_mem+threadIdx.x*(5*meqn+3*maux+mwaves+meqn*mwaves);//meqn
    double* qr   = ql+meqn;         //meqn
    double* qd   = qr+meqn;         //meqn
    double* auxl = qd+meqn;         //maux
    double* auxr = auxl+maux;       //maux
    double* auxd = auxr+maux;       //maux
    double* s    = auxd+maux;       //mwaves
    double* wave = s+mwaves;        //meqn*mwaves
    double* amdq = wave+meqn*mwaves;//meqn
    double* apdq = amdq+meqn;       //meqn

    int ifaces_x = mx+2*mbc-1;
    int ifaces_y = my+2*mbc-1;
    int num_ifaces = ifaces_x*ifaces_y;

    /* Compute strides */
    xs = 1;
    ys = (2*mbc + mx)*xs;
    zs = (2*mbc + my)*xs*ys;

    for(int thread_index = threadIdx.x; thread_index<num_ifaces; thread_index+=blockDim.x){

        int ix = thread_index%ifaces_x;
        int iy = thread_index/ifaces_y;

        /* (i,j) index */
        I = (iy + mbc-1)*ys + (ix + mbc-1)*xs;

        if (ix < mx + 2*mbc-1 && iy < my + 2*mbc-1)
        {
            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                ql[mq] = qold[I_q - xs];
                qr[mq] = qold[I_q];  
                qd[mq] = qold[I_q - ys];          
            }
            for(m = 0; m < maux; m++)
            {
                I_aux = I + m*zs;
                auxl[m] = aux[I_aux - xs];
                auxr[m] = aux[I_aux];
                auxd[m] = aux[I_aux - ys];
            }

            //rpn2adv(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);
            rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

            /* Set value at left interface of cell I */
            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                fp[I_q] = -apdq[mq]; 
                fm[I_q] = amdq[mq];
            }
#if 0        
            for (m = 0; m < meqn*mwaves; m++)
            {
                I_waves = I + m*zs;
                waves[I_waves] = wave[m];
            }
#endif        
            for (mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + mw*zs;
                speeds[I_speeds] = s[mw];
            } 


            rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, amdq, apdq);

            /* Set value at bottom interface of cell I */
            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                gp[I_q] = -apdq[mq]; 
                gm[I_q] = amdq[mq];
            }
            for (mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + (mwaves+mw)*zs;
                speeds[I_speeds] = s[mw];
            } 

#if 0        
            for (m = 0; m < meqn*mwaves; m++)
            {
                I_waves = I + m*zs;
                waves[I_waves] = wave[m];
            }
            for (mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + mw*zs;
                speeds[I_speeds] = s[mw];
            } 
#endif

        }
    }
    __syncthreads();
    for(int thread_index = threadIdx.x; thread_index<mx*my; thread_index+=blockDim.x){

        int ix = thread_index%mx;
        int iy = thread_index/my;

        I = (ix+mbc)*xs + (iy+mbc)*ys;

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx * (fm[I_q + xs] - fp[I_q]) 
                                  - dtdy * (gm[I_q + ys] - gp[I_q]);
        }        

    }
}

__global__ void cudaclaw_compute_cfl(int idir, int mx, int my, int meqn, int mwaves, 
                                      int mbc, double dx, double dy, double dt, 
                                      double *speeds, double* cflgrid)
{
#if 0    
      # from fortran_source/cudaclaw_flux2.f */

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


__device__ void cudaclaw_second_order(int idir, int mx, int my, int meqn, int mbc,
                                       int maux, double* qold, double* aux, double dx,
                                       double dy, double dt, double* cflgrid,
                                       double* fm, double* fp, double* gm, double* gp,
                                       double* waves, double *speeds,
                                       cudaclaw_cuda_rpn2_t rpn2, void* rpt2,
                                       int mwaves) 
{    
    int mq, mw, m;

    /* TODO : Limit waves here */


    /* TODO : Compute second order corrections */
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

