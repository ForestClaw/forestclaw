#include "../fc2d_cudaclaw.h"
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>

void fc2d_cudaclaw_store_buffer(fclaw2d_global_t* glob,
                                fclaw2d_patch_t *this_patch,
                                int this_patch_idx,
                                int count, int iter, 
                                cudaclaw_fluxes_t* flux_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);
    double *qold, *aux;
    int meqn, maux;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds;

    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
    fclaw2d_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    cudaEventRecord(start);

    {
        PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer:cudaMemcpy",5);
        cudaMemcpy(fluxes->qold_dev, qold, fluxes->num_bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(fluxes->aux_dev, aux, fluxes->num_bytes_aux, cudaMemcpyHostToDevice);
    }

    flux_array[iter % FC2D_CUDACLAW_BUFFER_LEN] = *fluxes;

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop); 
    glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY].cumulative += milliseconds*1e-3;

}
