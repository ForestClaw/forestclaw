/*
  Copyright (c) 2018 Carsten Burstedde, Donna Calhoun, Melody Shih, Scott Aiton, 
  Xinsheng Qin.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



#include "../fc2d_cudaclaw_cuda.h"

#include <math.h>
#include <cub/cub.cuh>   // or equivalently <cub/block/block_reduce.cuh>

#include "cudaclaw_allocate.h"

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

static
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
                                cudaclaw_cuda_rpt2_t rpt2,
                                cudaclaw_cuda_b4step2_t b4step2,
                                double t,double dt)
{
    /* Does this 128 have to match the 128 grid size used to launch this kernel? */
    typedef cub::BlockReduce<double,128> BlockReduce;

    __shared__ typename BlockReduce::TempStorage temp_storage;

    int mq, mw, m, k;
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
    double gupdate;
    int imp;

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
    double* aux1 = bpdq+meqn;         //2*maux
    double* aux2 = aux1+2*maux;       //2*maux
    double* aux3 = aux2+2*maux;       //2*maux
    double* bmasdq = aux3+2*maux;     //meqn
    double* bpasdq = bmasdq+meqn;     //meqn

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
                ql[mq] = qold[I_q - 1];
                qr[mq] = qold[I_q];  
                qd[mq] = qold[I_q - ys];          
            }

            for(m = 0; m < maux; m++)
            {
                I_aux = I + m*zs;
                auxl[m] = aux[I_aux - 1];
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
                I_q = I + mq*zs;
                amdq[mq] = fm[I_q];
                apdq[mq] = fp[I_q];
                bmdq[mq] = gm[I_q];
                bpdq[mq] = gp[I_q];
            }

            /* Limit waves */
            for(mw = 0; mw < mwaves; mw++)
            {
                /* X-faces */
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
                    fm[I_q] += 0.5*cqxx;   
                    fp[I_q] += 0.5*cqxx;   
                    amdq[mq] += cqxx;   /* Used for transverse waves */
                    apdq[mq] -= cqxx;      
                }

                /* Y-faces */
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
                    gm[I_q] += 0.5*cqyy;   
                    gp[I_q] += 0.5*cqyy;   
                    bmdq[mq] += cqyy;     /* Used for transverse waves */
                    bpdq[mq] -= cqyy;      
                }                
            }  /* End of mwaves loop */


            __syncthreads();

            /* ----------------------- Transverse : X-faces --------------------------- */

            for(imp = 0; imp < 2; imp++)
            {
                for(m = 0; m < maux; m++)
                {
                    I_aux = I + m*zs;
                    k = imp*maux + m;
                    aux1[k] = aux[I_aux - ys + (imp - 1)];
                    aux2[k] = aux[I_aux      + (imp - 1)];
                    aux3[k] = aux[I_aux + ys + (imp - 1)];
                }
            }

            /* idir = 0; imp = 0 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bmasdq[mq];
                gm[I_q - 1] -= gupdate;        /* Subtract 1 when imp=0 */
                gp[I_q - 1] -= gupdate;

                gupdate = 0.5*dtdx*bpasdq[mq];
                gm[I_q - 1 + ys] -= gupdate;
                gp[I_q - 1 + ys] -= gupdate;
            }

            /* idir = 0; imp = 1 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bmasdq[mq];
                gm[I_q] -= gupdate;        
                gp[I_q] -= gupdate;

                gupdate = 0.5*dtdx*bpasdq[mq];
                gm[I_q + ys] -= gupdate;
                gp[I_q + ys] -= gupdate;
            }
            

            /* ----------------------- Transverse : Y-faces --------------------------- */

            for(imp = 0; imp < 2; imp++)
            {
                for(m = 0; m < maux; m++)
                {
                    I_aux = I + m*zs;
                    k = imp*maux + m;
                    aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                    aux2[k] = aux[I_aux     + ys*(imp - 1)];
                    aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
                }
            }

            /* idir = 1; imp = 0;  Re-use amdq, apdq */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdy*bmasdq[mq];
                fm[I_q - ys] -= gupdate;        /* Subtract 1 when imp=0 */
                fp[I_q - ys] -= gupdate;

                gupdate = 0.5*dtdy*bpasdq[mq];
                fm[I_q - ys + 1] -= gupdate;
                fp[I_q - ys + 1] -= gupdate;
            }

            /* idir = 1; imp = 1;  Re-use amdp, apdq */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdy*bmasdq[mq];
                fm[I_q] -= gupdate;        
                fp[I_q] -= gupdate;

                gupdate = 0.5*dtdy*bpasdq[mq];
                fm[I_q + 1] -= gupdate;
                fp[I_q + 1] -= gupdate;
            }   

        } /* Thread conditional */
    } /* Thread loop */

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


/* ---------------------------------------------------------------------------------------
   PUBLIC function  
   ------------------------------------------------------------------------------------ */
__global__
void cudaclaw_flux2_and_update_batch (int mx, int my, int meqn, int mbc, 
                                int maux, int mwaves, int mwork,
                                double dt, double t,
                                cudaclaw_fluxes_t* array_fluxes_struct_dev,
                                double * maxcflblocks_dev,
                                cudaclaw_cuda_rpn2_t rpn2,
                                cudaclaw_cuda_rpt2_t rpt2,
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
            maxcflblocks_dev, rpn2, rpt2, b4step2,
            t,dt);
}



