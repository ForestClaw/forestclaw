#include "../fc2d_cudaclaw_cuda.h"
#include "cudaclaw_allocate.h" /* Needed for definition of fluxes */

#include <fc2d_cudaclaw_options.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <fc2d_cuda_profiler.h>

#include "../fc2d_cudaclaw_check.cu"


/* Static buffers for transfering to the GPU */
static double* s_membuffer;
static double* s_membuffer_dev;
static cudaclaw_fluxes_t* s_array_fluxes_struct_dev;


void cudaclaw_allocate_fluxes(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *patch)
{
    PROFILE_CUDA_GROUP("Allocate patch data in memory device",4);       
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    double value;

    cudaclaw_fluxes_t *fluxes = FCLAW_ALLOC(cudaclaw_fluxes,1);

    const fclaw2d_clawpatch_options_t *claw_opt = fclaw2d_clawpatch_get_options(glob);
    int meqn = claw_opt->meqn;
    int maux = claw_opt->maux;

    /* Set values needed in batch node */
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fluxes->dx = dx;
    fluxes->dy = dy;

    fluxes->xlower = xlower;
    fluxes->ylower = ylower;
          
    /* Set global arrays on device */
    fc2d_cudaclaw_options_t* cuda_opt = fc2d_cudaclaw_get_options(glob);
    int mwaves = cuda_opt->mwaves;


    size_t size = (2*mbc+mx)*(2*mbc+my);
 
    fluxes->num        = meqn*size;
    fluxes->num_aux    = maux*size;
    fluxes->num_waves  = 2*mwaves*meqn*size;
    fluxes->num_speeds = 2*mwaves*size;
 

    fluxes->num_bytes        = meqn*size*sizeof(double);
    fluxes->num_bytes_aux    = maux*size*sizeof(double);
    fluxes->num_bytes_waves  = 2*mwaves*meqn*size*sizeof(double);
    fluxes->num_bytes_speeds = 2*mwaves*size*sizeof(double);
    
    CHECK(cudaMalloc((void**)&fluxes->fm_dev,   fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->fp_dev,   fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->gm_dev,   fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->gp_dev,   fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->amdq_dev, fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->apdq_dev, fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->bmdq_dev, fluxes->num_bytes));
    CHECK(cudaMalloc((void**)&fluxes->bpdq_dev, fluxes->num_bytes));

    CHECK(cudaMalloc((void**)&fluxes->waves_dev,  fluxes->num_bytes_waves));
    CHECK(cudaMalloc((void**)&fluxes->speeds_dev, fluxes->num_bytes_speeds));

    /* Set all values to 0 so max works, even if boundary edge values are not assigned 
       speeds in kernel */
    value = 0;
    CHECK(cudaMemset((void*)fluxes->speeds_dev, value, fluxes->num_bytes_speeds));

    fclaw2d_patch_set_user_data(glob,patch,fluxes);
}

void cudaclaw_deallocate_fluxes(fclaw2d_global_t *glob,
                                     fclaw2d_patch_t *patch)
{
    PROFILE_CUDA_GROUP("De-allocate patch device memory",4);       
    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,patch);

    FCLAW_ASSERT(fluxes != NULL);

    /* Assumption here is that cudaFree is a synchronous call */
    CHECK(cudaFree(fluxes->fm_dev));
    CHECK(cudaFree(fluxes->fp_dev));
    CHECK(cudaFree(fluxes->gm_dev));
    CHECK(cudaFree(fluxes->gp_dev));

    CHECK(cudaFree(fluxes->amdq_dev));
    CHECK(cudaFree(fluxes->apdq_dev));
    CHECK(cudaFree(fluxes->bmdq_dev));
    CHECK(cudaFree(fluxes->bpdq_dev));
    
    CHECK(cudaFree(fluxes->waves_dev));
    CHECK(cudaFree(fluxes->speeds_dev));

    FCLAW_FREE((void*) fluxes);
}



void fc2d_cudaclaw_allocate_buffers(fclaw2d_global_t *glob)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int maux = clawpatch_opt->maux;
    int meqn = clawpatch_opt->meqn;  

    int batch_size = FC2D_CUDACLAW_BUFFER_LEN;
    size_t size = (2*mbc+mx)*(2*mbc+my);
    size_t bytes = batch_size*size*(meqn + maux)*sizeof(double);

    CHECK(cudaMallocHost(&s_membuffer,bytes));    
    CHECK(cudaMalloc(&s_membuffer_dev, bytes)); 
    CHECK(cudaMalloc(&s_array_fluxes_struct_dev, 
                     batch_size*sizeof(cudaclaw_fluxes_t)));
}

void fc2d_cudaclaw_deallocate_buffers(fclaw2d_global_t *glob)
{
    cudaFreeHost(s_membuffer);
    cudaFree(s_membuffer_dev);
    cudaFree(s_array_fluxes_struct_dev);
}

double *cudaclaw_get_cpu_membuffer()
{
    FCLAW_ASSERT(s_membuffer != NULL);
    return s_membuffer;
}

double *cudaclaw_get_gpu_membuffer()
{
    FCLAW_ASSERT(s_membuffer_dev != NULL);
    return s_membuffer_dev;
}

cudaclaw_fluxes_t* cudaclaw_get_flux_buffer()
{
    FCLAW_ASSERT(s_array_fluxes_struct_dev != NULL);
    return s_array_fluxes_struct_dev;
}






