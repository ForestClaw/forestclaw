#include "../fc2d_cudaclaw.h"
#include "cudaclaw_allocate.h"
#include "cudaclaw_update_q.h"
#include "cudaclaw_flux2.h"

#include <fc2d_cudaclaw_options.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include "../fc2d_cudaclaw_check.cu"  /* CHECK defined here */

#include <cublas_v2.h>

    
double cudaclaw_step2_batch(fclaw2d_global_t *glob,
        cudaclaw_fluxes_t* array_fluxes_struct, 
        int batch_size, double dt)
{
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds;
    int i;

    double maxcfl = 0.0;
    double dtdx, dtdy, s;

    FCLAW_ASSERT(batch_size !=0);

    /* To get patch-independent parameters */
    fc2d_cudaclaw_options_t *clawopt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    clawopt = fc2d_cudaclaw_get_options(glob);
    int mwaves = clawopt->mwaves;

    fc2d_cudaclaw_vtable_t*  cuclaw_vt = fc2d_cudaclaw_vt();
    FCLAW_ASSERT(cuclaw_vt->cuda_rpn2 != NULL);


    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int maux = clawpatch_opt->maux;
    int meqn = clawpatch_opt->meqn;

    cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[0]);
    size_t bytes = batch_size*(fluxes->num_bytes + fluxes->num_bytes_aux);

    for(i = 0; i < batch_size; i++)
    {
        cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
        cudaMemcpy(fluxes->qold_dev, fluxes->qold, fluxes->num_bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(fluxes->aux_dev, fluxes->aux, fluxes->num_bytes_aux, cudaMemcpyHostToDevice);
    }

    cudaclaw_fluxes_t* array_fluxes_struct_dev = NULL;
    cudaMalloc(&array_fluxes_struct_dev, batch_size*sizeof(cudaclaw_fluxes_t));

    CHECK(cudaMemcpy(array_fluxes_struct_dev, array_fluxes_struct, batch_size*sizeof(cudaclaw_fluxes_t), cudaMemcpyHostToDevice));

    dim3 block(128,1,1);
    dim3 grid(1,1,batch_size);

    size_t bytes_per_thread = sizeof(double)*(5*meqn+3*maux+mwaves+meqn*mwaves);

    cudaclaw_flux2_and_update_batch<<<grid,block,128*bytes_per_thread>>>(mx,my,meqn,
                                                                     mbc,maux,mwaves,dt,
                                                                     array_fluxes_struct_dev,
                                                                     cuclaw_vt->cuda_rpn2);

    cudaDeviceSynchronize();
    CHECK(cudaPeekAtLastError());

    /* -------------------------------- Compute CFL ------------------------------------*/ 
    int n = 2*mwaves*(2*mbc+mx)*(2*mbc+my);
    int maxidx;

    // TODO: batch cublas call
    cublasStatus_t stat;
    cublasHandle_t handle;
    cublasCreate(&handle);
    for (int i = 0; i < batch_size; ++i)
    {
        cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
        stat = cublasIdamax(handle,n, fluxes->speeds_dev,1,&maxidx);
        if (stat != CUBLAS_STATUS_SUCCESS) {
                printf ("cublasIdamax failed");
                cublasDestroy(handle);
                return EXIT_FAILURE;
        }
        dtdx = dt/fluxes->dx;
        dtdy = dt/fluxes->dy;

        double maxabsspeed_patch = 0.0;
        double maxcfl_patch = 0.0;
        cudaMemcpy(&maxabsspeed_patch,fluxes->speeds_dev+maxidx-1,sizeof(double),
                   cudaMemcpyDeviceToHost);
        s = fabs(maxabsspeed_patch);
        maxcfl_patch = maxidx < n/2 ? s*dtdx : s*dtdy;        

        maxcfl = max(maxcfl_patch,maxcfl);

        /* -------------------------- Copy q back to host ----------------------------------*/ 
        cudaEventRecord(start);

        cudaMemcpy(fluxes->qold, fluxes->qold_dev, fluxes->num_bytes, cudaMemcpyDeviceToHost);

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY].cumulative += milliseconds*1e-3;    
        /* ------------------------------ Clean up -----------------------------------------*/ 

    }
    cublasDestroy(handle);

    cudaFree(array_fluxes_struct_dev);
    return maxcfl;
}

