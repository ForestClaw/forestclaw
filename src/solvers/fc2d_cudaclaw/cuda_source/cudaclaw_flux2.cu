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

#include <sm_60_atomic_functions.h>

static
__device__ double minmod(double r)
{
    return max(0.0,min(1.0,r));
}

static
__device__ double limiter(int lim_choice, double r)
{
    switch(lim_choice)
    {
        case 1:
            return minmod(r);

        default:
            return 1;  /* No limiting */
    }
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
                                double* amdq_g, double* apdq_g, 
                                double* bmdq_g, double* bpdq_g,
                                double* waves, double *speeds,
								double * maxcflblocks_dev,
                                cudaclaw_cuda_rpn2_t rpn2,
                                cudaclaw_cuda_rpt2_t rpt2,
                                cudaclaw_cuda_b4step2_t b4step2,
                                int* order, int* mthlim,
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


    if (b4step2 != NULL)
    {
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
                    qr[mq] = qold[I_q];  
                }

                for(m = 0; m < maux; m++)
                {
                    /* In case aux is already set */
                    I_aux = I + m*zs;
                    auxr[m] = aux[I_aux];
                }                        
                
                i = ix-(mbc-2);  /* i,j for index in the grid */
                j = iy-(mbc-2);
                b4step2(mbc,mx,my,meqn,qr,xlower,ylower,dx,dy, 
                        t,dt,maux,auxr,i,j);

                for(m = 0; m < maux; m++)
                {
                    /* In case aux is set by b4step2 */
                    I_aux = I + m*zs;
                    aux[I_aux] = auxr[m];
                }
            } /* Thread conditional */
        } /* Thread loop */

        __syncthreads();

    } /* b4step2 conditional */


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
                ql[mq] = qold[I_q - 1];          /* Left  */
                qr[mq] = qold[I_q];              /* Right */
                qd[mq] = qold[I_q - ys];         /* Down  */  
            }

            for(m = 0; m < maux; m++)
            {
                /* How is this getting set? */
                I_aux = I + m*zs;
                auxl[m] = aux[I_aux - 1];
                auxr[m] = aux[I_aux];
                auxd[m] = aux[I_aux - ys];
            }                        

            rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq);

            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                fp[I_q] = -apdq[mq]; 
                fm[I_q] = amdq[mq];
                if (order[0] == 2)
                {
                    apdq_g[I_q] = apdq[mq];  
                    amdq_g[I_q] = amdq[mq];
                }
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

            rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, bmdq, bpdq);

            /* Set value at bottom interface of cell I */
            for (mq = 0; mq < meqn; mq++) 
            {
                I_q = I + mq*zs;
                gm[I_q] = bmdq[mq];
                gp[I_q] = -bpdq[mq]; 
                if (order[0] == 2)
                {
                    bpdq_g[I_q] = bpdq[mq];
                    bmdq_g[I_q] = bmdq[mq];                    
                }
            }

            if (order[0] == 2)
            {
                /* Waves are only needed for second order corrections */
                for (m = 0; m < meqn*mwaves; m++)
                {
                    I_waves = I + (meqn*mwaves+m)*zs;
                    waves[I_waves] = wave[m];
                }                
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

    if (order[0] == 2)
    {
        for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            ix = thread_index % ifaces_x;
            iy = thread_index/ifaces_y;

            I = (ix + mbc)*xs + (iy + mbc)*ys;

            if (ix < mx + 1 && iy < my + 1)   /* Is this needed? */
            {
                /* Limit waves */
                for(mw = 0; mw < mwaves; mw++)
                {
                    /* X-faces */

                    I_speeds = I + mw*zs;
                    s[mw] = speeds[I_speeds];

                    wlimitr = 1;
                    if (mthlim[mw] > 0)
                    {
                        wnorm2 = dotl = dotr = 0;
                        for(mq = 0; mq < meqn; mq++)
                        {
                            I_waves = I + (mw*meqn + mq)*zs;
                            wave[mq] = waves[I_waves];
                            wnorm2 += pow(wave[mq],2);
                            dotl += wave[mq]*waves[I_waves-1];
                            dotr += wave[mq]*waves[I_waves+1];
                        }
                        if (wnorm2 != 0)
                        {
                            r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                            wlimitr = limiter(mthlim[mw],r);  
                        }
                    }

 
                    for(mq = 0; mq < meqn; mq++)
                    {
                        I_q = I + mq*zs;
                        cqxx = fabs(s[mw])*(1.0 - fabs(s[mw])*dtdx)*wlimitr*wave[mq];
                        fm[I_q] += 0.5*cqxx;   
                        fp[I_q] += 0.5*cqxx;                               
                        if (order[1] == 2)
                        {
                            /* Propagate second order corrections 
                               in transverse dir. */
                            amdq_g[I_q] += cqxx;   
                            apdq_g[I_q] -= cqxx;      
                        }
                    }

                    /* Y-faces */

                    I_speeds = I + (mwaves + mw)*zs;
                    s[mw] = speeds[I_speeds];

                    wlimitr = 1;
                    if (mthlim[mw] > 0)
                    {
                        wnorm2 = dotl = dotr = 0;
                        for(mq = 0; mq < meqn; mq++)
                        {
                            I_waves = I + ((mwaves+mw)*meqn + mq)*zs;
                            wave[mq] = waves[I_waves];
                            wnorm2 += pow(wave[mq],2);
                            dotl += wave[mq]*waves[I_waves-ys];
                            dotr += wave[mq]*waves[I_waves+ys];
                        }
                        if (wnorm2 != 0)
                        {
                            r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                            wlimitr = limiter(mthlim[mw],r);  
                        }
                    }

                    for(mq = 0; mq < meqn; mq++)
                    {
                        I_q = I + mq*zs;
                        cqyy = fabs(s[mw])*(1.0 - fabs(s[mw])*dtdy)*wlimitr*wave[mq];
                        gm[I_q] += 0.5*cqyy;   
                        gp[I_q] += 0.5*cqyy;   
                        if (order[1] == 2)
                        {
                            /* Propagate second order corrections 
                               in transverse dir. */
                            bmdq_g[I_q] += cqyy;     
                            bpdq_g[I_q] -= cqyy;      
                        }
                    }                
                }  /* End of mwaves loop */
            } /* End of thread conditional */
        } /* End of thread loop */
    } /* End of check on order[0] == 2 */


    __syncthreads();


    if (order[1] == 0)
    {
        /* Update the solution and exit */
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
        /* Return from here and don't do any transverse propagation */
        return;
    }


    /* ------------------------ Transverse Propagation : X-faces ---------------------- */
    
    /*     transverse-x

            |     |     | 
            |     |     | 
        ----|-----|-----|-----
            |     X     | 
            |  v--X     |
        ----|--O--|-----|-----
            |     |     |
            |     |     |

    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        /* Lower left face */
        if (0 < ix && ix < mx + 1 && iy < my)   /* Is this needed? */
        {
            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                amdq[mq] = amdq_g[I_q];
                apdq[mq] = apdq_g[I_q];
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

            /* idir = 0; imp = 0 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,0,amdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bmasdq[mq];
                gm[I_q - 1] -= gupdate;       
                gp[I_q - 1] -= gupdate;
            }
        } /* thread conditional */
    } /* Thread loop */

    __syncthreads();

    /*   transverse-x  
            |     |     | 
            |     |     | 
        ----|--O--|-----|----
            |  ^__X     | 
            |     X     |
        ----|-----|-----|----
            |     |     |
            |     |     |

    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (0 < ix && ix < mx + 1 && iy < my)   /* Is this needed? */
        {
            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                amdq[mq] = amdq_g[I_q];
                apdq[mq] = apdq_g[I_q];
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

            /* idir = 0; imp = 0 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,1,amdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bpasdq[mq];
                gm[I_q - 1 + ys] -= gupdate;
                gp[I_q - 1 + ys] -= gupdate;
            }
        } /* thread conditional */
    } /* Thread loop */

    __syncthreads();

    /*  transverse-x
            |     |     | 
            |     |     | 
        ----|-----|-----|----
            |     X     | 
            |     X--v  |
        ----|-----|--O--|----
            |     |     |
            |     |     |

    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && iy < my)   
        {

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                amdq[mq] = amdq_g[I_q];
                apdq[mq] = apdq_g[I_q];
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

            /* idir = 0; imp = 0 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,0,apdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bmasdq[mq];
                gm[I_q] -= gupdate;       
                gp[I_q] -= gupdate;
            }
        } /* thread conditional */
    } /* Thread loop */

    __syncthreads();

    /*  transverse-x 
            |     |     | 
            |     |     | 
        ----|-----|--O--|----
            |     X__^  | 
            |     X     |
        ----|-----|-----|----
            |     |     |
            |     |     |

    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && iy < my)   
        {

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                amdq[mq] = amdq_g[I_q];
                apdq[mq] = apdq_g[I_q];
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

            /* idir = 0; imp = 0 */
            rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,1,apdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdx*bpasdq[mq];
                gm[I_q + ys] -= gupdate;
                gp[I_q + ys] -= gupdate;
            }
        } /* thread conditional */
    } /* Thread loop */

    __syncthreads();

    /* ----------------------------- Transverse : Y-faces ----------------------------- */


    /*  transverse-y

             |     |     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
             |     |     
        -----|-XXX-|-----
             |  v  |     
             0--   |     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && 0 < iy && iy < my+1)   /* Is this needed? */
        {

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                bmdq[mq] = bmdq_g[I_q];
                bpdq[mq] = bpdq_g[I_q];
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

            /* idir = 1; imp = 0 */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,0,bmdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdy*bmasdq[mq];
                fm[I_q - ys] -= gupdate;        
                fp[I_q - ys] -= gupdate;
            }
        } /* Thread conditional */
    } /* Thread loop */

    __syncthreads();

    /*  transverse-y

             |     |     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
             |     |     
        -----|-XXX-|-----
             |  v  |     
             |   --O     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && 0 < iy && iy < my+1)   /* Is this needed? */
        {
            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                bmdq[mq] = bmdq_g[I_q];
                bpdq[mq] = bpdq_g[I_q];
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

            /* idir = 1; imp = 0 */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,1,bmdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdy*bpasdq[mq];
                fm[I_q - ys + 1] -= gupdate;
                fp[I_q - ys + 1] -= gupdate;
            }
        } /* Thread conditional */
    } /* Thread loop */

    __syncthreads();


    /*  transverse-y

             |     |     
             |     |     
        -----|-----|-----
             |     |     
             O---  |     
             |  ^  |     
        -----|-XXX-|-----
             |     |     
             |     |     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
    */ 

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && iy < my)   /* Is this needed? */
        {
            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                bmdq[mq] = bmdq_g[I_q];
                bpdq[mq] = bpdq_g[I_q];
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

            /* idir = 1; imp = 1 */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,0,bpdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
                gupdate = 0.5*dtdy*bmasdq[mq];
                fm[I_q] -= gupdate;        
                fp[I_q] -= gupdate;
            }   
        } /* Thread conditional */
    } /* Thread loop */

    __syncthreads();

    /*  transverse-y

             |     |     
             |     |     
        -----|-----|-----
             |     |     
             |  ---O     
             |  ^  |     
        -----|-XXX-|-----
             |     |     
             |     |     
             |     |     
        -----|-----|-----
             |     |     
             |     |     
    */              

    for(thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        ix = thread_index % ifaces_x;
        iy = thread_index/ifaces_y;

        I = (ix + mbc)*xs + (iy + mbc)*ys;

        if (ix < mx && iy < my)   /* Is this needed? */
        {

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;
                bmdq[mq] = bmdq_g[I_q];
                bpdq[mq] = bpdq_g[I_q];
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

            /* idir = 1; imp = 1 */
            rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,1,bpdq,bmasdq,bpasdq);

            for(mq = 0; mq < meqn; mq++)
            {
                I_q = I + mq*zs;  
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
                                int* order, int* mthlim,
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
            array_fluxes_struct_dev[blockIdx.z].amdq_dev,
            array_fluxes_struct_dev[blockIdx.z].apdq_dev,
            array_fluxes_struct_dev[blockIdx.z].bmdq_dev,
            array_fluxes_struct_dev[blockIdx.z].bpdq_dev,
            array_fluxes_struct_dev[blockIdx.z].waves_dev,
            array_fluxes_struct_dev[blockIdx.z].speeds_dev,
            maxcflblocks_dev, rpn2, rpt2, b4step2,
            order, mthlim, t,dt);
}



