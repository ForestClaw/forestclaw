#include "../fc2d_cudaclaw_cuda.h"
#include "cudaclaw_allocate.h" /* Needed for definition of fluxes */

#include <fc2d_cudaclaw_options.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <fc2d_cuda_profiler.h>

#include "../fc2d_cudaclaw_check.cu"


void cudaclaw_allocate_fluxes(fclaw2d_global_t *glob,
                               fclaw2d_patch_t *patch)
{
    PROFILE_CUDA_GROUP("Allocate patch data in memory device",4);       
    int mx,my,mbc;
    double xlower,ylower,dx,dy;

    const fclaw2d_clawpatch_options_t *claw_opt = fclaw2d_clawpatch_get_options(glob);
    int meqn = claw_opt->meqn;
    int maux = claw_opt->maux;

    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_cudaclaw_options_t* cuda_opt = fc2d_cudaclaw_get_options(glob);
    int mwaves = cuda_opt->mwaves;

    cudaclaw_fluxes_t *fluxes = FCLAW_ALLOC(cudaclaw_fluxes,1);

    size_t size = (2*mbc+mx)*(2*mbc+my);
 
    fluxes->num        = meqn*size;
    fluxes->num_aux    = maux*size;
    fluxes->num_waves  = 2*mwaves*meqn*size;
    fluxes->num_speeds = 2*mwaves*size;
 

    fluxes->num_bytes        = meqn*size*sizeof(double);
    fluxes->num_bytes_aux    = maux*size*sizeof(double);
    fluxes->num_bytes_waves  = 2*mwaves*meqn*size*sizeof(double);
    fluxes->num_bytes_speeds = 2*mwaves*size*sizeof(double);
    
    fluxes->dx = dx;
    fluxes->dy = dy;

    fluxes->xlower = xlower;
    fluxes->ylower = ylower;
          
#if 0    
    CHECK(cudaMalloc((void**)&fluxes->qold_dev,   fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->aux_dev,    fluxes->num_bytes_aux));
#endif    

    CHECK(cudaMalloc((void**)&fluxes->fm_dev,     fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->fp_dev,     fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->gm_dev,     fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->gp_dev,     fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->waves_dev,  fluxes->num_bytes_waves));
    CHECK(cudaMalloc((void**)&fluxes->speeds_dev, fluxes->num_bytes_speeds));
    CHECK(cudaMemset((void*)fluxes->speeds_dev, 0, fluxes->num_bytes_speeds));

    fclaw2d_patch_set_user_data(glob,patch,fluxes);
}

void cudaclaw_deallocate_fluxes(fclaw2d_global_t *glob,
                                 fclaw2d_patch_t *patch)
{
    PROFILE_CUDA_GROUP("De-allocate patch data",4);       
    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,patch);

    FCLAW_ASSERT(fluxes != NULL);

    /* Assumption here is that cudaFree is a synchronous call */
#if 0    
    CHECK(cudaFree(fluxes->qold_dev));
    CHECK(cudaFree(fluxes->aux_dev));
#endif    
    CHECK(cudaFree(fluxes->fm_dev));
    CHECK(cudaFree(fluxes->fp_dev));
    CHECK(cudaFree(fluxes->gm_dev));
    CHECK(cudaFree(fluxes->gp_dev));
    CHECK(cudaFree(fluxes->waves_dev));
    CHECK(cudaFree(fluxes->speeds_dev));

    FCLAW_FREE((void*) fluxes);
}

