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
#include <fc2d_cuda_profiler.h>
#include <cub/cub.cuh>

    
double cudaclaw_step2_batch(fclaw2d_global_t *glob,
        cudaclaw_fluxes_t* array_fluxes_struct, 
        int batch_size, double dt)
{
    PROFILE_CUDA_GROUP("cudaclaw_step2_batch",5);
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

    /* ---------------------------------- Merge Memory ---------------------------------*/ 
    cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[0]);
    size_t size = batch_size*(fluxes->num + fluxes->num_aux);
    size_t bytes = size*sizeof(double);

    double *membuffer = FCLAW_ALLOC(double,size);

    double* membuffer_dev;
    CHECK(cudaMalloc((void**)&membuffer_dev, bytes));

    size_t memoffset = 0;
    for(i = 0; i < batch_size; i++)   
    {
        cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
        int I_q = memoffset;
        int I_aux = I_q + fluxes->num;

        memcpy(&membuffer[I_q],fluxes->qold,fluxes->num_bytes);
        memcpy(&membuffer[I_aux],fluxes->aux,fluxes->num_bytes_aux);

        /* Assign gpu pointers */
        fluxes->qold_dev = membuffer_dev + memoffset;
        fluxes->aux_dev = fluxes->qold_dev + fluxes->num;

        memoffset += (fluxes->num + fluxes->num_aux);
    }        
    FCLAW_ASSERT(memoffset == size);

    CHECK(cudaMemcpy(membuffer_dev, membuffer, bytes, cudaMemcpyHostToDevice));

    FCLAW_FREE(membuffer);

    /* -------------------------------- Work with array --------------------------------*/ 

    cudaclaw_fluxes_t* array_fluxes_struct_dev = NULL;
    CHECK(cudaMalloc(&array_fluxes_struct_dev, batch_size*sizeof(cudaclaw_fluxes_t)));

    // printf("sizeof(struct) = %d\n",sizeof(cudaclaw_fluxes_t));

    CHECK(cudaMemcpy(array_fluxes_struct_dev, array_fluxes_struct, 
                     batch_size*sizeof(cudaclaw_fluxes_t), 
                     cudaMemcpyHostToDevice));

    dim3 block(128,1,1);
    dim3 grid(1,1,batch_size);

    size_t bytes_per_thread = sizeof(double)*(5*meqn+3*maux+mwaves+meqn*mwaves);
    
    double* maxcflblocks_dev;
    cudaMalloc(&maxcflblocks_dev,batch_size*sizeof(double)); 
    cudaclaw_flux2_and_update_batch<<<grid,block,128*bytes_per_thread >>>(mx,my,meqn,
                                                                     mbc,maux,mwaves,dt,
                                                                     array_fluxes_struct_dev,
								                                     maxcflblocks_dev,
                                                                     cuclaw_vt->cuda_rpn2);
    cudaDeviceSynchronize();
    CHECK(cudaPeekAtLastError());
	
    /* -------------------------------- Finish CFL ------------------------------------*/ 
    void    *temp_storage_dev = NULL;
    size_t  temp_storage_bytes = 0;
    double  *cflgrid_dev;

    cudaMalloc(&cflgrid_dev, sizeof(double));  
    CubDebugExit(cub::DeviceReduce::Max(temp_storage_dev,temp_storage_bytes,
                                        maxcflblocks_dev,cflgrid_dev,batch_size));
    cudaMalloc(&temp_storage_dev, temp_storage_bytes);
    CubDebugExit(cub::DeviceReduce::Max(temp_storage_dev,temp_storage_bytes,
                                        maxcflblocks_dev,cflgrid_dev,batch_size));
    cudaMemcpy(&maxcfl, cflgrid_dev, sizeof(double),cudaMemcpyDeviceToHost);
    cudaFree(temp_storage_dev);
    cudaFree(cflgrid_dev);
    /* ------------------------------ Done with CFL ------------------------------------*/ 

    /* -------------------------- Copy q back to host ----------------------------------*/ 
    for (int i = 0; i < batch_size; ++i)    
    {      
        cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
        cudaMemcpy(fluxes->qold, fluxes->qold_dev, fluxes->num_bytes, 
                   cudaMemcpyDeviceToHost);
    }

    /* ------------------------------ Clean up -----------------------------------------*/ 
    cudaFree(array_fluxes_struct_dev);
    cudaFree(membuffer_dev);
    return maxcfl;
}

