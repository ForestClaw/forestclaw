#include "cudaclaw5_allocate.h"

#include <fc2d_cudaclaw5_options.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_timer.h>



void cudaclaw5_allocate_fluxes(struct fclaw_global *glob,
                               struct fclaw2d_patch *patch)
{
    const fclaw_clawpatch_options_t *claw_opt = fclaw_clawpatch_get_options(glob);
    int mx = claw_opt->mx;
    int my = claw_opt->my;
    int mbc = claw_opt->mbc;
    int meqn = claw_opt->meqn;
    int maux = claw_opt->maux;

    fc2d_cudaclaw5_options_t* cuda_opt = fc2d_cudaclaw5_get_options(glob);
    int mwaves = cuda_opt->mwaves;

    cudaclaw5_fluxes_t *fluxes = FCLAW_ALLOC(cudaclaw5_fluxes,1);

    size_t size = (2*mbc+mx)*(2*mbc+my)*sizeof(double);
    fluxes->num_bytes        = meqn*size;
    fluxes->num_bytes_aux    = maux*size;
    fluxes->num_bytes_waves  = mwaves*meqn*size;
    fluxes->num_bytes_speeds = mwaves*size;

    /* Assumption here is that cudaMalloc is a synchronous call */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]); 
          
    cudaMalloc((void**)&fluxes->qold_dev,   fluxes->num_bytes);
    cudaMalloc((void**)&fluxes->fm_dev,     fluxes->num_bytes);
    cudaMalloc((void**)&fluxes->fp_dev,     fluxes->num_bytes);
    cudaMalloc((void**)&fluxes->gm_dev,     fluxes->num_bytes);
    cudaMalloc((void**)&fluxes->gp_dev,     fluxes->num_bytes);
    cudaMalloc((void**)&fluxes->aux_dev,    fluxes->num_bytes_aux);
    cudaMalloc((void**)&fluxes->waves_dev,  fluxes->num_bytes_waves);
    cudaMalloc((void**)&fluxes->speeds_dev, fluxes->num_bytes_speeds);

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]);    

    fclaw2d_patch_set_user_data(glob,patch,fluxes);
}

void cudaclaw5_deallocate_fluxes(fclaw_global_t *glob,
                                 fclaw_patch_t *patch)
{
    cudaclaw5_fluxes_t *fluxes = (cudaclaw5_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,patch);

    FCLAW_ASSERT(fluxes != NULL);

    /* Assumption here is that cudaFree is a synchronous call */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]);       
    cudaFree(fluxes->qold_dev);
    cudaFree(fluxes->fm_dev);
    cudaFree(fluxes->fp_dev);
    cudaFree(fluxes->gm_dev);
    cudaFree(fluxes->gp_dev);
    cudaFree(fluxes->aux_dev);
    cudaFree(fluxes->waves_dev);
    cudaFree(fluxes->speeds_dev);
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]);    

    FCLAW_FREE((void*) fluxes);
}

