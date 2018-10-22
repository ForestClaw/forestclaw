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

#include "../fc2d_cudaclaw_check.cu"

#include <cub/block/block_reduce.cuh>  

#include "cudaclaw_allocate.h"  /* Needed to for definition of 'fluxes' */

__constant__ int order[2];
__constant__ int mthlim[FC2D_CUDACLAW_MWAVES];

extern "C"
{
int cudaclaw_check_parameters(int mwaves)
{
    return mwaves <= FC2D_CUDACLAW_MWAVES;
}

void cudaclaw_set_method_parameters(int *order_in, int *mthlim_in, int mwaves)
{
    CHECK(cudaMemcpyToSymbol(order,order_in,2*sizeof(int)));
    CHECK(cudaMemcpyToSymbol(mthlim,mthlim_in,mwaves*sizeof(int)));
}

}

/* Include this here so we don't include device code in fc2d_cudaclaw_cuda.h */
__device__ double cudaclaw_limiter(int lim_choice, double r);


static
__device__
void cudaclaw_flux2_and_update(int mx, int my, int meqn, int mbc,
                               int maux, int mwaves, int mwork,
                               double xlower, double ylower, 
                               double dx, double dy,
                               double* qold, double* aux, 
                               volatile double* fm, volatile double* fp, 
                               volatile double* gm, volatile double* gp,
                               double* amdq_trans, double* apdq_trans, 
                               double* bmdq_trans, double* bpdq_trans,
                               double* waves, double *speeds,
                               double * maxcflblocks,
                               cudaclaw_cuda_rpn2_t rpn2,
                               cudaclaw_cuda_rpt2_t rpt2,
                               cudaclaw_cuda_b4step2_t b4step2,
                               double t,double dt)
{
    typedef cub::BlockReduce<double,FC2D_CUDACLAW_BLOCK_SIZE> BlockReduce;

    __shared__ typename BlockReduce::TempStorage temp_storage;


    extern __shared__ double shared_mem[];

    double* start  = shared_mem + mwork*threadIdx.x;

    double* ql     = start;                 /* meqn        */
    double* qr     = ql     + meqn;         /* meqn        */
    double* qd     = qr     + meqn;         /* meqn        */
    double* auxl   = qd     + meqn;         /* maux        */
    double* auxr   = auxl   + maux;         /* maux        */
    double* auxd   = auxr   + maux;         /* maux        */
    double* s      = auxd   + maux;         /* mwaves      */
    double* wave   = s      + mwaves;       /* meqn*mwaves */
    double* amdq   = wave   + meqn*mwaves;  /* meqn        */
    double* apdq   = amdq   + meqn;         /* meqn        */
    double* bmdq   = apdq   + meqn;         /* meqn        */
    double* bpdq   = bmdq   + meqn;         /* meqn        */
    double* aux1   = bpdq   + meqn;         /* 2*maux      */
    double* aux2   = aux1   + 2*maux;       /* 2*maux      */
    double* aux3   = aux2   + 2*maux;       /* 2*maux      */
    double* bmasdq = aux3   + 2*maux;       /* meqn        */
    double* bpasdq = bmasdq + meqn;         /* meqn        */

    double dtdx, dtdy;

    int mq, mw, m, k;
    int xs, ys, zs;
    int I, I_q, I_aux, I_waves, I_speeds;
    int ix,iy,ifaces_x, ifaces_y, num_ifaces;
    int thread_index;

    int i,j; /* Used for (i,j) indexing in patches  */
    double maxcfl;
    double wnorm2,dotr,dotl, wlimitr,r;
    double cqxx;
    double cqyy;
    double gupdate;
    int imp;


    /* --------------------------------- Start code ----------------------------------- */


    dtdx = dt/dx;
    dtdy = dt/dy;

    /* Compute strides */
    xs = 1;
    ys = (2*mbc + mx)*xs;
    zs = (2*mbc + my)*xs*ys;

    maxcfl = 0;

    ifaces_x = mx + 2*mbc-1;
    ifaces_y = my + 2*mbc-1;
    num_ifaces = ifaces_x*ifaces_y;

    if (b4step2 != NULL)
    {
        for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        {
            ix = thread_index % ifaces_x;
            iy = thread_index/ifaces_x;

            I = (iy + 1)*ys + (ix + 1)*xs;  /* Start at one cell from left/bottom */

            //if (ix < mx + 2*mbc-1 && iy < my + 2*mbc-1)
            {
                for(mq = 0; mq < meqn; mq++)
                {
                    I_q = I + mq*zs;
                    qr[mq] = qold[I_q];  
                }

                for(m = 0; m < maux; m++)
                {
                    /* In case aux is already set */
                    I_aux = I + m*zs;
                    auxr[m] = aux[I_aux];
                }          

                /* Compute (i,j) for patch index (i,j) in (1-mbc:mx+mbc,1-mbc:my+mbc)

                      i + (mbc-2) == ix   Check : i  = 1 --> ix = mbc-1
                      j + (mbc-2) == iy   Check : ix = 0 --> i  = 2-mbc
                */
                i = ix-(mbc-2);  
                j = iy-(mbc-2);
                b4step2(mbc,mx,my,meqn,qr,xlower,ylower,dx,dy, 
                        t,dt,maux,auxr,i,j);

                for(m = 0; m < maux; m++)
                {
                    /* In case aux is set by b4step2 */
                    I_aux = I + m*zs;
                    aux[I_aux] = auxr[m];
                }
            } 
        } 
        __syncthreads(); /* Needed to be sure all aux variables are available below */ 
    } 

    /* -------------------------- Compute fluctuations -------------------------------- */

    ifaces_x = mx + 2*mbc-1;
    ifaces_y = my + 2*mbc-1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I = (iy + 1)*ys + (ix + 1);  /* Start one cell from left/bottom edge */

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;
            qr[mq] = qold[I_q];        /* Right */
            ql[mq] = qold[I_q - 1];    /* Left  */
            qd[mq] = qold[I_q - ys];   /* Down  */  
        }

        for(m = 0; m < maux; m++)
        {
            I_aux = I + m*zs;
            auxl[m] = aux[I_aux - 1];
            auxr[m] = aux[I_aux];
            auxd[m] = aux[I_aux - ys];
        }                        

        /* ------------------------ Normal solve in X direction ------------------- */
            
        rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

        for (mq = 0; mq < meqn; mq++) 
        {
            I_q = I + mq*zs;
            fm[I_q] = amdq[mq];
            fp[I_q] = -apdq[mq]; 
            if (order[1] > 0)
            {
                amdq_trans[I_q] = amdq[mq];                                        
                apdq_trans[I_q] = apdq[mq];  
            }
        }

        for(mw = 0; mw < mwaves; mw++)
        {
            maxcfl = max(maxcfl,abs(s[mw]*dtdx));

            if (order[0] == 2)
            {                    
                I_speeds = I + mw*zs;
                speeds[I_speeds] = s[mw];
                for(mq = 0; mq < meqn; mq++)
                {
                    k = mw*meqn + mq;
                    I_waves = I + k*zs;
                    waves[I_waves] = wave[k];
                }
            }
        }
        
        /* ------------------------ Normal solve in Y direction ------------------- */
        rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, bmdq, bpdq);

        /* Set value at bottom interface of cell I */
        for (mq = 0; mq < meqn; mq++) 
        {
            I_q = I + mq*zs;
            gm[I_q] = bmdq[mq];
            gp[I_q] = -bpdq[mq]; 
            if (order[1] > 0)
            {
                bpdq_trans[I_q] = bpdq[mq];
                bmdq_trans[I_q] = bmdq[mq];                                                   
            }
        }

        for(mw = 0; mw < mwaves; mw++)
        {
            maxcfl = max(maxcfl,fabs(s[mw])*dtdy);

            if (order[0] == 2)
            {                    
                I_speeds = I + (mwaves + mw)*zs;
                speeds[I_speeds] = s[mw];
                for(mq = 0; mq < meqn; mq++)
                {
                    I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                    waves[I_waves] = wave[mw*meqn + mq];
                }
            }
        }
    }
    

    maxcflblocks[blockIdx.z] = BlockReduce(temp_storage).Reduce(maxcfl,cub::Max());

    //__syncthreads();  /* Does block reduce take care of this sync? */


    /* ---------------------- Second order corrections and limiters --------------------*/  
    
    if (order[0] == 2)
    {
        ifaces_x = mx + 1;  
        ifaces_y = my + 2;
        num_ifaces = ifaces_x*ifaces_y;

        for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            ix = thread_index % ifaces_x;
            iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            I = (iy + mbc-1)*ys + (ix + mbc);

            /* ------------------------------- X-directions --------------------------- */
            for(mw = 0; mw < mwaves; mw++)
            {
                I_speeds = I + mw*zs;
                s[mw] = speeds[I_speeds];

                for(mq = 0; mq < meqn; mq++)
                {
                    I_waves = I + (mw*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    wnorm2 = dotl = dotr = 0;
                    for(mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);

                        I_waves = I + (mw*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-1];
                        dotr += wave[mq]*waves[I_waves+1];
                    }
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;
                                               
                    r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                    wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(mq = 0; mq < meqn; mq++)
                {
                    I_q = I + mq*zs;
                    cqxx = abs(s[mw])*(1.0 - abs(s[mw])*dtdx)*wave[mq];

                    fm[I_q] += 0.5*cqxx;   
                    fp[I_q] += 0.5*cqxx;  
                    if (order[1] > 0)
                    {                         
                        amdq_trans[I_q] += cqxx;   
                        apdq_trans[I_q] -= cqxx;                                 
                    }  
                }
            }
        }

        //__syncthreads();

        ifaces_x = mx + 2;  
        ifaces_y = my + 1;
        num_ifaces = ifaces_x*ifaces_y;

        for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            ix = thread_index % ifaces_x;
            iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            I = (iy + mbc)*ys + ix + mbc - 1;

            /* ------------------------------- X-directions --------------------------- */
            for(mw = 0; mw < mwaves; mw++)
            {
                /* ------------------------------- Y-directions --------------------------- */
                I_speeds = I + (mwaves + mw)*zs;
                s[mw] = speeds[I_speeds];

                for(mq = 0; mq < meqn; mq++)
                {
                    I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    wnorm2 = dotl = dotr = 0;
                    for(mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);
                        I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-ys];
                        dotr += wave[mq]*waves[I_waves+ys];
                    }  
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;  

                    r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;

                    wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(mq = 0; mq < meqn; mq++)
                {
                    I_q = I + mq*zs;
                    cqyy = abs(s[mw])*(1.0 - abs(s[mw])*dtdy)*wave[mq];

                    gm[I_q] += 0.5*cqyy;   
                    gp[I_q] += 0.5*cqyy;  

                    if (order[1] > 0)
                    {                            
                        bmdq_trans[I_q] += cqyy;     
                        bpdq_trans[I_q] -= cqyy;      
                    } 
                }   
            }  
        }  
        __syncthreads();
    }  /* Done with second order corrections */


    if (order[1] == 0)
    {
        /* No transverse propagation; Update the solution and exit */
        for(thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
        {
            ix = thread_index % mx;
            iy = thread_index/mx;

            I = (ix + mbc)*xs + (iy + mbc)*ys;

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
                                      - dtdy * (gm[I_q + ys] - gp[I_q]);
            }        
        }
        return;
    }


    /* ------------------------ Transverse Propagation : X-faces ---------------------- */


    ifaces_x = mx + 1;  /* Visit x - edges of all non-ghost cells */
    ifaces_y = my + 2;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

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


        /*     transverse-x
    
                |     |     | 
                |     |     | 
            ----|-----|-----|-----
                |     X     | 
                |     X  q  |
                |  v--X     |
            ----|--O--|-----|-----
                |     |     |
                |     |     |
    
        */              

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdx*bmasdq[mq];
            gm[I_q - 1] -= gupdate;       
            gp[I_q - 1] -= gupdate;   
        }            
    }

    __syncthreads();


    ifaces_x = mx + 1;  /* Visit x - edges of all non-ghost cells */
    ifaces_y = my + 2;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

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


        /*     transverse-x
    
                |     |     | 
                |     |     | 
            ----|--0--|-----|-----
                |  ^--X     | 
                |     X  q  |
                |     X     |
            ----|--O--|-----|-----
                |     |     |
                |     |     |
    
        */              

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q - 1 + ys] -= gupdate;
            gp[I_q - 1 + ys] -= gupdate;
        }        
    }


    __syncthreads();

    ifaces_x = mx + 1;  /* Visit x - edges of all non-ghost cells */
    ifaces_y = my + 2;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

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

        /*     transverse-x
    
                |     |     | 
                |     |     | 
            ----|-----|-----|-----
                |     X     | 
                |     X  q  |
                |     X--v  |
            ----|-----|--0--|-----
                |     |     |
                |     |     |
    
        */              
        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdx*bmasdq[mq];
            gm[I_q] -= gupdate;       
            gp[I_q] -= gupdate;
        }
    }

    __syncthreads();

    ifaces_x = mx + 1;  /* Visit x - edges of all non-ghost cells */
    ifaces_y = my + 2;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

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


        /*     transverse-x
    
                |     |     | 
                |     |     | 
            ----|-----|--0--|-----
                |     X--^  | 
                |     X  q  |
                |     X     |
            ----|-----|-----|-----
                |     |     |
                |     |     |
    
        */              

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q + ys] -= gupdate;
            gp[I_q + ys] -= gupdate;
        }
        
    } 

    /* May not the synchthreads(), below, since gm/gp updated above, but only fm/gp
       updated below */
    __syncthreads();  

    
    /* ----------------------------- Transverse : Y-faces ----------------------------- */


    /*  transverse-y

             |     |     
        -----|-----|-----
             |     |     
             O-- --O      
             | ^ ^ |     
        -----|-XXX-|-----
             | v v |     
             0-- --0     
             |     |     
        -----|-----|-----
             |     |     
    */              

    ifaces_x = mx + 2;  /* Visit edges of all non-ghost cells */
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc)*ys + (ix + mbc-1);


        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

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

        /*  transverse-y
    
                 |     |     
            -----|-----|-----
                 |     |     
                 |  q  |      
                 |     |     
            -----|-XXX-|-----
                 | v   |     
                 0--   0     
                 |     |     
            -----|-----|-----
                 |     |     
        */              


        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q - ys] -= gupdate;        
            fp[I_q - ys] -= gupdate;
        }

    }

    ifaces_x = mx + 2;  /* Visit edges of all non-ghost cells */
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc)*ys + (ix + mbc-1);


        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

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


        /*  transverse-y
    
                 |     |     
            -----|-----|-----
                 |     |     
                 |  q  |      
                 |     |     
            -----|-XXX-|-----
                 |   v |     
                 0   --0     
                 |     |     
            -----|-----|-----
                 |     |     
        */              
        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q - ys + 1] -= gupdate;
            fp[I_q - ys + 1] -= gupdate;                
        }
        
    }

    __syncthreads();

    ifaces_x = mx + 2;  /* Visit edges of all non-ghost cells */
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc)*ys + (ix + mbc-1);


        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

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


        /*  transverse-y
    
                 |     |     
            -----|-----|-----
                 |  q  |     
                 O--   |           
                 | ^   |     
            -----|-XXX-|-----
                 |     |     
                 |     | 
                 |     |     
            -----|-----|-----
                 |     |     
        */              

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q] -= gupdate;        
            fp[I_q] -= gupdate;
        }
    }

    __syncthreads();

    ifaces_x = mx + 2;  /* Visit edges of all non-ghost cells */
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_x;

        I =  (iy + mbc)*ys + (ix + mbc-1);


        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

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

        /*  transverse-y
    
                 |     |     
            -----|-----|-----
                 |  q  |     
                 |   --0           
                 |   ^ |     
            -----|-XXX-|-----
                 |     |     
                 |     | 
                 |     |     
            -----|-----|-----
                 |     |     
        */              


        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq);

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;  
            gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q + 1] -= gupdate;
            fp[I_q + 1] -= gupdate;
        }   
        
    } 

    __syncthreads();

    /* ------------------------------- Final update ----------------------------------- */

    for(thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    {
        ix = thread_index % mx;
        iy = thread_index/my;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        for(mq = 0; mq < meqn; mq++)
        {
            I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
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
                                      cudaclaw_fluxes_t* array_fluxes_struct,
                                      double * maxcflblocks,
                                      cudaclaw_cuda_rpn2_t rpn2,
                                      cudaclaw_cuda_rpt2_t rpt2,
                                      cudaclaw_cuda_b4step2_t b4step2)
    {
        cudaclaw_flux2_and_update(mx,my,meqn,mbc,maux,mwaves,mwork,
                                  array_fluxes_struct[blockIdx.z].xlower,
                                  array_fluxes_struct[blockIdx.z].ylower,
                                  array_fluxes_struct[blockIdx.z].dx,
                                  array_fluxes_struct[blockIdx.z].dy,
                                  array_fluxes_struct[blockIdx.z].qold_dev,
                                  array_fluxes_struct[blockIdx.z].aux_dev,
                                  array_fluxes_struct[blockIdx.z].fm_dev,
                                  array_fluxes_struct[blockIdx.z].fp_dev,
                                  array_fluxes_struct[blockIdx.z].gm_dev,
                                  array_fluxes_struct[blockIdx.z].gp_dev,
                                  array_fluxes_struct[blockIdx.z].amdq_dev,
                                  array_fluxes_struct[blockIdx.z].apdq_dev,
                                  array_fluxes_struct[blockIdx.z].bmdq_dev,
                                  array_fluxes_struct[blockIdx.z].bpdq_dev,
                                  array_fluxes_struct[blockIdx.z].waves_dev,
                                  array_fluxes_struct[blockIdx.z].speeds_dev, 
                                  maxcflblocks, rpn2, rpt2, b4step2, t,dt);
}



