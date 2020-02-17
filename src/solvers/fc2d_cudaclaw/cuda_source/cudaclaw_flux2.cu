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
__constant__ int use_fwaves;

extern "C"
{
int cudaclaw_check_parameters(int mwaves)
{
    return mwaves <= FC2D_CUDACLAW_MWAVES;
}

void cudaclaw_set_method_parameters(int *order_in, int *mthlim_in, int mwaves, 
                                    int use_fwaves_in)
{
    CHECK(cudaMemcpyToSymbol(order,order_in,2*sizeof(int)));
    CHECK(cudaMemcpyToSymbol(use_fwaves,&use_fwaves_in,sizeof(int)));
    CHECK(cudaMemcpyToSymbol(mthlim,mthlim_in,mwaves*sizeof(int)));
}

}

/* Include this here so we don't include device code in fc2d_cudaclaw_cuda.h */
__device__ double cudaclaw_limiter(int lim_choice, double r);


static
__device__
void cudaclaw_flux2_and_update(const int mx,   const int my, 
                               const int meqn, const int mbc,
                               const int maux, const int mwaves, 
                               const int mwork,
                               const double xlower, const double ylower, 
                               const double dx,     const double dy,
                               double *const qold,       double *const aux, 
                               double *const fm,         double *const fp, 
                               double *const gm,         double *const gp,
                               double *const amdq_trans, double *const apdq_trans, 
                               double *const bmdq_trans, double *const bpdq_trans,
                               double *const waves,      double *const speeds,
                               double *const maxcflblocks,
                               cudaclaw_cuda_rpn2_t rpn2,
                               cudaclaw_cuda_rpt2_t rpt2,
                               cudaclaw_cuda_b4step2_t b4step2,
                               double t,double dt)
{
    typedef cub::BlockReduce<double,FC2D_CUDACLAW_BLOCK_SIZE> BlockReduce;
    
    __shared__ typename BlockReduce::TempStorage temp_storage;

    extern __shared__ double shared_mem[];

    double* start  = shared_mem + mwork*threadIdx.x;

#if 0
    double *const ql     = start;                 /* meqn        */
    double *const qr     = ql     + meqn;         /* meqn        */
    double *const qd     = qr     + meqn;         /* meqn        */
    double *const auxl   = qd     + meqn;         /* maux        */
    double *const auxr   = auxl   + maux;         /* maux        */
    double *const auxd   = auxr   + maux;         /* maux        */
    double *const s      = auxd   + maux;         /* mwaves      */
    double *const wave   = s      + mwaves;       /* meqn*mwaves */
    double *const amdq   = wave   + meqn*mwaves;  /* meqn        */
    double *const apdq   = amdq   + meqn;         /* meqn        */
    double *const bmdq   = apdq   + meqn;         /* meqn        */
    double *const bpdq   = bmdq   + meqn;         /* meqn        */
    double *const aux1   = bpdq   + meqn;         /* 2*maux      */
    double *const aux2   = aux1   + 2*maux;       /* 2*maux      */
    double *const aux3   = aux2   + 2*maux;       /* 2*maux      */
    double *const bmasdq = aux3   + 2*maux;       /* meqn        */
    double *const bpasdq = bmasdq + meqn;         /* meqn        */
#endif    


    /* --------------------------------- Start code ----------------------------------- */

    __shared__ double dtdx, dtdy;
    __shared__ int xs,ys,zs;
    if (threadIdx.x == 0)
    {
        dtdx = dt/dx;
        dtdy = dt/dy;        

        /* Compute strides */
        xs = 1;
        ys = (2*mbc + mx)*xs;
        zs = (2*mbc + my)*xs*ys;
    }
    
    __syncthreads();

    double maxcfl = 0;

    int ifaces_x = mx + 2*mbc-1;
    int ifaces_y = my + 2*mbc-1;
    int num_ifaces = ifaces_x*ifaces_y;

#if 0
    if (b4step2 != NULL)
    {
        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        {
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            int I = (iy + 1)*ys + (ix + 1);  /* Start at one cell from left/bottom */

            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qr[mq] = qold[I_q];  
            }

            for(int m = 0; m < maux; m++)
            {
                /* In case aux is already set */
                int I_aux = I + m*zs;
                auxr[m] = aux[I_aux];
            }          

            /* Compute (i,j) for patch index (i,j) in (1-mbc:mx+mbc,1-mbc:my+mbc)
    
                i + (mbc-2) == ix   Check : i  = 1 --> ix = mbc-1
                j + (mbc-2) == iy   Check : ix = 0 --> i  = 2-mbc
            */
            int i = ix-(mbc-2);  
            int j = iy-(mbc-2);
            b4step2(mbc,mx,my,meqn,qr,xlower,ylower,dx,dy, 
                    t,dt,maux,auxr,i,j);

            for(int m = 0; m < maux; m++)
            {
                /* In case aux is set by b4step2 */
                int I_aux = I + m*zs;
                aux[I_aux] = auxr[m];
            }
        }      
        __syncthreads(); /* Needed to be sure all aux variables are available below */ 
    } 
#endif    

    /* -------------------------- Compute fluctuations -------------------------------- */

#if 0
    ifaces_x = mx + 2*mbc-1;
    ifaces_y = my + 2*mbc-1;
    num_ifaces = ifaces_x*ifaces_y;
#endif    

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I = (iy + 1)*ys + (ix + 1);  /* Start one cell from left/bottom edge */

        {
            /* ------------------------ Normal solve in X direction ------------------- */
            double *const ql     = start;                 /* meqn        */
            double *const qr     = ql     + meqn;         /* meqn        */
            double *const auxl   = qr     + meqn;         /* maux        */
            double *const auxr   = auxl   + maux;         /* maux        */
            double *const s      = auxr   + maux;         /* mwaves      */
            double *const wave   = s      + mwaves;       /* meqn*mwaves */
            double *const amdq   = wave   + meqn*mwaves;  /* meqn        */
            double *const apdq   = amdq   + meqn;         /* meqn        */
            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qr[mq] = qold[I_q];        /* Right */
                ql[mq] = qold[I_q - 1];    /* Left  */
            }

            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                auxr[m] = aux[I_aux];
                auxl[m] = aux[I_aux - 1];
            }               

            
            rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

            for (int mq = 0; mq < meqn; mq++) 
            {
                int I_q = I + mq*zs;
                fm[I_q] = amdq[mq];
                fp[I_q] = -apdq[mq]; 
                if (order[1] > 0)
                {
                    amdq_trans[I_q] = amdq[mq];                                        
                    apdq_trans[I_q] = apdq[mq];  
                }
            }
        

            for(int mw = 0; mw < mwaves; mw++)
            {
                maxcfl = max(maxcfl,abs(s[mw]*dtdx));

                if (order[0] == 2)
                {                    
                    int I_speeds = I + mw*zs;
                    speeds[I_speeds] = s[mw];
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        int k = mw*meqn + mq;
                        int I_waves = I + k*zs;
                        waves[I_waves] = wave[k];
                    }
                }
            }
        }

        {
            /* ------------------------ Normal solve in Y direction ------------------- */
            double *const qr     = start;                 /* meqn        */
            double *const qd     = qr     + meqn;         /* meqn        */
            double *const auxr   = qd     + meqn;         /* maux        */
            double *const auxd   = auxr   + maux;         /* maux        */
            double *const s      = auxd   + maux;         /* mwaves      */
            double *const wave   = s      + mwaves;       /* meqn*mwaves */
            double *const bmdq   = wave   + meqn*mwaves;  /* meqn        */
            double *const bpdq   = bmdq   + meqn;         /* meqn        */

            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qr[mq] = qold[I_q];        /* Right */
                qd[mq] = qold[I_q - ys];   /* Down  */  
            }

            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                auxr[m] = aux[I_aux];
                auxd[m] = aux[I_aux - ys];
            }               

            rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, bmdq, bpdq);

            /* Set value at bottom interface of cell I */
            for (int mq = 0; mq < meqn; mq++) 
            {
                int I_q = I + mq*zs;
                gm[I_q] = bmdq[mq];
                gp[I_q] = -bpdq[mq]; 
                if (order[1] > 0)
                {
                    bmdq_trans[I_q] = bmdq[mq];                                                   
                    bpdq_trans[I_q] = bpdq[mq];
                }
            }

            for(int mw = 0; mw < mwaves; mw++)
            {
                maxcfl = max(maxcfl,fabs(s[mw])*dtdy);

                if (order[0] == 2)
                {                    
                    int I_speeds = I + (mwaves + mw)*zs;
                    speeds[I_speeds] = s[mw];
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                        waves[I_waves] = wave[mw*meqn + mq];
                    }
                }
            }
        }
    }


    maxcflblocks[blockIdx.z] = BlockReduce(temp_storage).Reduce(maxcfl,cub::Max());

    //__syncthreads();  /* Does block reduce take care of this sync? */


    /* ---------------------- Second order corrections and limiters --------------------*/  
    
    if (order[0] == 2)
    {
        int ifaces_x = mx + 1;  
        int ifaces_y = my + 2;
        int num_ifaces = ifaces_x*ifaces_y;

        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            int I = (iy + mbc-1)*ys + (ix + mbc);

            /* ------------------------------- X-directions --------------------------- */

            double *s = start;
            double *wave = s + mwaves;
            for(int mw = 0; mw < mwaves; mw++)
            {                
                int I_speeds = I + mw*zs;
                s[mw] = speeds[I_speeds];

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_waves = I + (mw*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    double wnorm2 = 0, dotl=0, dotr = 0;
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);

                        int I_waves = I + (mw*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-1];
                        dotr += wave[mq]*waves[I_waves+1];
                    }
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;
                                               
                    double r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                    double wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (int mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(int mq = 0; mq < meqn; mq++)
                {
                    double cqxx = (1.0 - fabs(s[mw])*dtdx)*wave[mq];
                    cqxx *= (use_fwaves) ? copysign(1.,s[mw]) : fabs(s[mw]);

                    int I_q = I + mq*zs;
                    fm[I_q] += 0.5*cqxx;   
                    fp[I_q] += 0.5*cqxx;  
                    if (order[1] > 0)
                    {                         
                        amdq_trans[I_q] += cqxx;   
                        apdq_trans[I_q] -= cqxx;      /* Subtract cqxx later */                         
                    }  
                }
            }
        }

        ifaces_x = mx + 2;  
        ifaces_y = my + 1;
        num_ifaces = ifaces_x*ifaces_y;

        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            int I = (iy + mbc)*ys + ix + mbc - 1;

            double *const s = start;
            double *const wave = s + mwaves;
            for(int mw = 0; mw < mwaves; mw++)
            {
                /* ------------------------------- Y-directions --------------------------- */
                int I_speeds = I + (mwaves + mw)*zs;
                s[mw] = speeds[I_speeds];

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    double wnorm2 = 0, dotl = 0, dotr = 0;
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);
                        int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-ys];
                        dotr += wave[mq]*waves[I_waves+ys];
                    }  
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;  

                    double r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;

                    double wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (int mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_q = I + mq*zs;
                    double cqyy = (1.0 - fabs(s[mw])*dtdx)*wave[mq];
                    cqyy *= (use_fwaves) ? copysign(1.,s[mw]) : fabs(s[mw]);

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
    }  

    /* ------------------------- First order final update ----------------------------- */

    if (order[1] == 0)
    {
        /* No transverse propagation; Update the solution and exit */
        for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
        {
            int ix = thread_index % mx;
            int iy = thread_index/mx;

            int I = (ix + mbc) + (iy + mbc)*ys;

            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
                                      - dtdy * (gm[I_q + ys] - gp[I_q]);
            }        
        }
        return;
    }


    /* ------------------------ Transverse Propagation : X-faces ---------------------- */

    ifaces_x = mx + 1;            /* Match edges visited by Clawpack */
    ifaces_y = my + 2;
    num_ifaces = ifaces_x*ifaces_y;


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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */
        

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const amdq   = ql + meqn;      /* meqn   */
        double *const aux1   = amdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            //apdq[mq] = apdq_trans[I_q];
        }            
        
        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }            

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq);
        

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bmasdq[mq];
            gm[I_q - 1] -= gupdate;       
            gp[I_q - 1] -= gupdate;   
        }            
    }
    __syncthreads();


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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const amdq   = ql + meqn;      /* meqn   */
        double *const aux1   = amdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
            //apdq[mq] = apdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q - 1 + ys] -= gupdate;
            gp[I_q - 1 + ys] -= gupdate;
        }        
    }
    __syncthreads();


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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const apdq   = ql + meqn;      /* meqn   */
        double *const aux1   = apdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            //amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bmasdq[mq];
            gm[I_q] -= gupdate;       
            gp[I_q] -= gupdate;
        }
    }
    __syncthreads();


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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const apdq   = ql + meqn;      /* meqn   */
        double *const aux1   = apdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            //amdq[mq] = amdq_trans[I_q];
            apdq[mq] = apdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q + ys] -= gupdate;
            gp[I_q + ys] -= gupdate;
        }
        
    } 
    __syncthreads();

    
    /* ----------------------------- Transverse : Y-faces ----------------------------- */

    ifaces_x = mx + 2;  /* Visit edges of all non-ghost cells */
    ifaces_y = my + 1;
    num_ifaces = ifaces_x*ifaces_y;

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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bmdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bmdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            //bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q - ys] -= gupdate;        
            fp[I_q - ys] -= gupdate;
        }

    }
    __syncthreads();


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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bmdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bmdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
            //bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            int gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q - ys + 1] -= gupdate;
            fp[I_q - ys + 1] -= gupdate;                
        }
    }
    __syncthreads();



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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bpdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bpdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            //bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q] -= gupdate;        
            fp[I_q] -= gupdate;
        }
    }
    __syncthreads();

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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bpdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bpdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            //bmdq[mq] = bmdq_trans[I_q];
            bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q + 1] -= gupdate;
            fp[I_q + 1] -= gupdate;
        }   
    } 

    __syncthreads();

    /* ------------------------------- Final update ----------------------------------- */

    for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    {
        int ix = thread_index % mx;
        int iy = thread_index/mx;

        int I = (ix + mbc) + (iy + mbc)*ys;

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
                                  - dtdy * (gm[I_q +ys] - gp[I_q]);
        }        
    }
}


/* ---------------------------------------------------------------------------------------
   PUBLIC function  
   ------------------------------------------------------------------------------------ */
__global__
void cudaclaw_flux2_and_update_batch (const int mx,    const int my, 
                                      const int meqn,  const int mbc, 
                                      const int maux,  const int mwaves, 
                                      const int mwork,
                                      const double dt, const double t,
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



