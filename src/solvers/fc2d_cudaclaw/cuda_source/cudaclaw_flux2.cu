#include "cudaclaw_flux2.h"
#include "cudaclaw_allocate.h"
#include <math.h>
#include <cub/cub.cuh>   // or equivalently <cub/block/block_reduce.cuh>

static
__device__ double minmod(double r)
{
    return max(0.0,min(1.0,r));
}

static
__device__ double limiter(double r)
{
    return minmod(r);
}

__device__
void cudaclaw_flux2_and_update(int mx, int my, int meqn, int mbc,
                                int maux, int mwaves, int mwork,
                                double xlower, double ylower, 
                                double dx, double dy,
                                double* qold, double* aux, 
                                double* fm, double* fp, 
                                double* gm, double* gp,
                                double* waves, double *speeds,
								double * maxcflblocks_dev,
                                cudaclaw_cuda_rpn2_t rpn2,
                                cudaclaw_cuda_b4step2_t b4step2,
                                double t,double dt)
{
    /* Does this 128 have to match the 128 grid size used to launch this kernel? */
    typedef cub::BlockReduce<double,128> BlockReduce;

    __shared__ typename BlockReduce::TempStorage temp_storage;

    int mq, mw, m;
    int xs, ys, zs;
    int I, I_q, I_aux, I_waves, I_speeds;
    int thread_index;
    int ix,iy,ifaces_x, ifaces_y, num_ifaces;

    int i,j; /* Used for (i,j) indexing in patches numbers */
    double dtdx, dtdy;
    double maxcfl, cfl;
    double wnorm2,dotr,dotl, wlimitr,r;
    double cqxx;
    double cqyy;

    extern __shared__ double shared_mem[];
    double* start = shared_mem + mwork*threadIdx.x;
    double* ql   = start;             //meqn
    double* qr   = ql+meqn;           //meqn
    double* qd   = qr+meqn;           //meqn
    double* auxl = qd+meqn;           //maux
    double* auxr = auxl+maux;         //maux
    double* auxd = auxr+maux;         //maux
    double* s    = auxd+maux;         //mwaves
    double* wave = s+mwaves;          //meqn*mwaves
    double* amdq = wave+meqn*mwaves;  //meqn
    double* apdq = amdq+meqn;         //meqn
    double* bmdq = apdq+meqn;         //meqn
    double* bpdq = bmdq+meqn;         //meqn

    ifaces_x = mx+2*mbc-1;
    ifaces_y = my+2*mbc-1;
    num_ifaces = ifaces_x*ifaces_y;

    dtdx = dt/dx;
    dtdy = dt/dy;

    /* Compute strides */
    xs = 1;
    ys = (2*mbc + mx)*xs;
    zs = (2*mbc + my)*xs*ys;

    maxcfl = 0;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

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
            
            if (b4step2 != NULL)
            {
                i = ix-(mbc-2);  /* i,j for index in the grid */
                j = iy-(mbc-2);
                b4step2(mbc,mx,my,meqn,qr,xlower,ylower,dx,dy, 
                    t,dt,maux,auxr,i,j);//tperiod = 4.0                
            }

            rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

            /* Set value at left interface of cell I */
            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                fp[I_q] = -apdq[mq]; 
                fm[I_q] = amdq[mq];
            }

            for (m = 0; m < meqn*mwaves; m++)
            {
                I_waves = I + m*zs;
                waves[I_waves] = wave[m];
            }

            for (mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + mw*zs;
                speeds[I_speeds] = s[mw];
                cfl = abs(s[mw]*dtdx);
                if (cfl > maxcfl)
                {
                    maxcfl = cfl;
                }
            } 

            rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, amdq, apdq);

            /* Set value at bottom interface of cell I */
            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                gp[I_q] = -apdq[mq]; 
                gm[I_q] = amdq[mq];
            }
            for (m = 0; m < meqn*mwaves; m++)
            {
                I_waves = I + (meqn*mwaves+m)*zs;
                waves[I_waves] = wave[m];
            }

            for (mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + (mwaves + mw)*zs;
                speeds[I_speeds] = s[mw];
                cfl = fabs(s[mw])*dtdy;
                if (cfl > maxcfl)
                {
                    maxcfl = cfl;
                }
            } 
        }
    }

    maxcflblocks_dev[blockIdx.z] = BlockReduce(temp_storage).Reduce(maxcfl,cub::Max());

    __syncthreads();

    /* ---------------------------------- Limit waves --------------------------------------*/  
    
    ifaces_x = mx + 1;
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx + 1 && iy < my + 1)
        {
            for(mq=0; mq < meqn; mq++)
            {
                I_q = I + mq;
                amdq[mq] = fm[I_q];
                apdq[mq] = fp[I_q];
                bmdq[mq] = gm[I_q];
                bpdq[mq] = gp[I_q];
            }
            for(mw = 0; mw < mwaves; mw++)
            {
                /* x-faces */
                wnorm2 = dotl = dotr = 0;
                for(mq = 0; mq < meqn; mq++)
                {
                    I_waves = I + (mw*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                    wnorm2 += pow(wave[mq],2);
                    dotl += wave[mq]*waves[I_waves-1];
                    dotr += wave[mq]*waves[I_waves+1];
                }
                I_speeds = I + mw*zs;
                s[mw] = speeds[I_speeds];
                wlimitr = 1;
                if (wnorm2 != 0)
                {
                    r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                    wlimitr = limiter(r);  /* allow for selection */
                }

                for(mq = 0; mq < meqn; mq++)
                {
                    I_q = I + mq*zs;
                    cqxx = fabs(s[mw])*(1.0 - fabs(s[mw])*dtdx)*wlimitr*wave[mq];
                    fm[I_q] += 0.5*cqxx;   /* amdq + cqxx */
                    fp[I_q] += 0.5*cqxx;   /* apdq - cqxx */
                }


                /* y-faces */
                wnorm2 = dotl = dotr = 0;
                for(mq = 0; mq < meqn; mq++)
                {
                    I_waves = I + ((mwaves+mw)*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                    wnorm2 += pow(wave[mq],2);
                    dotl += wave[mq]*waves[I_waves-ys];
                    dotr += wave[mq]*waves[I_waves+ys];
                }
                I_speeds = I + (mwaves + mw)*zs;
                s[mw] = speeds[I_speeds];
                wlimitr = 1;
                if (wnorm2 != 0)
                {
                    r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                    wlimitr = limiter(r);  /* allow for selection */
                }

                for(mq = 0; mq < meqn; mq++)
                {
                    I_q = I + mq*zs;
                    cqyy = fabs(s[mw])*(1.0 - fabs(s[mw])*dtdy)*wlimitr*wave[mq];
                    gm[I_q] += 0.5*cqyy;   /* amdq + cqxx */
                    gp[I_q] += 0.5*cqyy;   /* apdq - cqxx */
                }
            }
        }
    }

    __syncthreads();

    for(thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    {
        ix = thread_index % mx;
        iy = thread_index/my;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx * (fm[I_q + xs] - fp[I_q]) 
                                  - dtdy * (gm[I_q + ys] - gp[I_q]);
        }        
    }
}



__global__
void cudaclaw_flux2_and_update_batch (int mx, int my, int meqn, int mbc, 
                                int maux, int mwaves, int mwork,
                                double dt, double t,
                                cudaclaw_fluxes_t* array_fluxes_struct_dev,
                                double * maxcflblocks_dev,
                                cudaclaw_cuda_rpn2_t rpn2,
                                cudaclaw_cuda_b4step2_t b4step2)
{
    // TODO: check this device function does not depend on blockIdx.z inside
    cudaclaw_flux2_and_update(mx,my,meqn,mbc,maux,mwaves,mwork,
            array_fluxes_struct_dev[blockIdx.z].xlower,
            array_fluxes_struct_dev[blockIdx.z].ylower,
            array_fluxes_struct_dev[blockIdx.z].dx,
            array_fluxes_struct_dev[blockIdx.z].dy,
            array_fluxes_struct_dev[blockIdx.z].qold_dev,
            array_fluxes_struct_dev[blockIdx.z].aux_dev,
            array_fluxes_struct_dev[blockIdx.z].fm_dev,
            array_fluxes_struct_dev[blockIdx.z].fp_dev,
            array_fluxes_struct_dev[blockIdx.z].gm_dev,
            array_fluxes_struct_dev[blockIdx.z].gp_dev,
            array_fluxes_struct_dev[blockIdx.z].waves_dev,
            array_fluxes_struct_dev[blockIdx.z].speeds_dev,
            maxcflblocks_dev, rpn2, b4step2,
            t,dt);
}



#if 0
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
#endif

