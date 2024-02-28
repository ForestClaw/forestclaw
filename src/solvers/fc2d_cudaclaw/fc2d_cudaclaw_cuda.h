/*
  Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#ifndef FC2D_CUDACLAW_CUDA_H
#define FC2D_CUDACLAW_CUDA_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* linear block length (set by configuration)  */
#if !defined(FC2D_CUDACLAW_BLOCK_SIZE)
#define FC2D_CUDACLAW_BLOCK_SIZE 128
#endif    


#define FC2D_CUDACLAW_MWAVES  10           /* Used to set shared memory (checked) */ 


struct fclaw2d_global;
struct fclaw2d_patch;
struct cudaclaw_fluxes;


/* ------------------------------ Typdefs for CUDA device functions --------------------*/

typedef void (*cudaclaw_cuda_b4step2_t)(int mbc, int mx, int my, int meqn, double q[],
                                        double xlower, double ylower, 
                                        double dx, double dy, 
                                        double time, double dt, int maux, 
                                        double aux[], int i, int j);

typedef void (*cudaclaw_cuda_src2_t)(int meqn, int maux, 
                                     double xlower, double ylower,
                                     double dx, double dy, 
                                     double q[], double aux[],
                                     double time, double dt, 
                                     int i, int j);

typedef void (*cudaclaw_cuda_rpn2_t)(int idir, int meqn, int mwaves, int maux,
                                     double ql[], double qr[], 
                                     double auxl[], double auxr[],
                                     double wave[], double s[], 
                                     double amdq[], double apdq[]);

typedef void (*cudaclaw_cuda_rpt2_t)(int idir, int meqn, int mwaves, int maux,
                                     double ql[], double qr[], 
                                     double aux1[], double aux2[], double aux3[],
                                     int imp, double dsdq[],
                                     double bmasdq[], double bpasdq[]);

typedef void (*cudaclaw_cuda_speeds_t)(int idir, int meqn, int mwaves, int maux,
                                       double ql[], double  qr[],
                                       double auxl[], double auxr[],
                                       double s[]);


/* ------------------------------------- Function headers ------------------------------*/

void cudaclaw_allocate_fluxes(struct fclaw2d_global *glob,
                              struct fclaw2d_patch *patch);
    
void cudaclaw_deallocate_fluxes(struct fclaw2d_global *glob,
                                struct fclaw2d_patch *patch);


double cudaclaw_step2_batch(struct fclaw2d_global* glob,
                            struct cudaclaw_fluxes* fluxes_array,
                            int patch_buffer_len, double t, double dt);

void cudaclaw_store_buffer(struct fclaw2d_global* glob,
                           struct fclaw2d_patch *this_patch,
                           int patchno, int blockno,
                           int total, int iter,
                           struct cudaclaw_fluxes* flux_array,
                           struct fclaw2d_patch** patch_array,
                           int* patchno_array,
                           int* blockno_array);
                          

double *cudaclaw_get_cpu_membuffer();

double *cudaclaw_get_gpu_membuffer();

struct cudaclaw_fluxes* cudaclaw_get_flux_buffer();

/* --------------------------- Function headers (used outside) -------------------------*/

void fc2d_cudaclaw_allocate_buffers(struct fclaw2d_global *glob);  /* Done once */

void fc2d_cudaclaw_deallocate_buffers(struct fclaw2d_global *glob);

void fc2d_cudaclaw_initialize_GPUs(struct fclaw2d_global *glob);


void  cudaclaw_get_method_parameters(int** order, int** mthlim);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
