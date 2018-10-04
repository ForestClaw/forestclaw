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
        int batch_size, double t, double dt)
{
    PROFILE_CUDA_GROUP("cudaclaw_step2_batch",5);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds;
    int i;

    double maxcfl = 0.0;
    //double dtdx, dtdy, s;

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
    size_t size = batch_size*(fluxes->num + fluxes->num_aux);
    size_t bytes = size*sizeof(double);
    double *membuffer;
    double* membuffer_dev;

    /* ---------------------------------- Merge Memory ---------------------------------*/ 
    {
        PROFILE_CUDA_GROUP("cudaclaw_copy_loop",7);    
        CHECK(cudaMallocHost((void**)&membuffer,bytes));

        CHECK(cudaMalloc((void**)&membuffer_dev, bytes));

        for(i = 0; i < batch_size; i++)   
        {
            cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);    

            int I_q = i*fluxes->num;
            int I_aux = batch_size*fluxes->num + i*fluxes->num_aux;

            memcpy(&membuffer[I_q]  ,fluxes->qold ,fluxes->num_bytes);
            memcpy(&membuffer[I_aux],fluxes->aux  ,fluxes->num_bytes_aux);

            /* Assign gpu pointers */
            fluxes->qold_dev = &membuffer_dev[I_q];
            fluxes->aux_dev  = &membuffer_dev[I_aux];
        }        

        CHECK(cudaMemcpy(membuffer_dev, membuffer, bytes, cudaMemcpyHostToDevice));
    }        

    /* -------------------------------- Work with array --------------------------------*/ 

    cudaclaw_fluxes_t* array_fluxes_struct_dev = NULL;
    CHECK(cudaMalloc(&array_fluxes_struct_dev, batch_size*sizeof(cudaclaw_fluxes_t)));

    CHECK(cudaMemcpy(array_fluxes_struct_dev, array_fluxes_struct, 
                     batch_size*sizeof(cudaclaw_fluxes_t), 
                     cudaMemcpyHostToDevice));

    dim3 block(128,1,1);
    dim3 grid(1,1,batch_size);

    size_t bytes_per_thread = sizeof(double)*(5*meqn+3*maux+mwaves+meqn*mwaves);
    
    double* maxcflblocks_dev;
    CHECK(cudaMalloc(&maxcflblocks_dev,batch_size*sizeof(double))); 
    cudaclaw_flux2_and_update_batch<<<grid,block,128*bytes_per_thread >>>(mx,my,meqn,
                                                                     mbc,maux,mwaves,dt,t,
                                                                     array_fluxes_struct_dev,
								                                     maxcflblocks_dev,
                                                                     cuclaw_vt->cuda_rpn2,
                                                                     cuclaw_vt->cuda_b4step2);
    cudaDeviceSynchronize();
    CHECK(cudaPeekAtLastError());
	
    /* -------------------------------- Finish CFL ------------------------------------*/ 
    {
        PROFILE_CUDA_GROUP("Finish CFL",1);
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
    }
    /* ------------------------------ Done with CFL ------------------------------------*/ 

    /* -------------------------- Copy q back to host ----------------------------------*/ 
    CHECK(cudaMemcpy(membuffer, membuffer_dev, batch_size*fluxes->num_bytes, 
                     cudaMemcpyDeviceToHost));

    {
        PROFILE_CUDA_GROUP("Copy back to patches loop",2);
        for (int i = 0; i < batch_size; ++i)    
        {      

            cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
            int I_q = i*fluxes->num;

            memcpy(fluxes->qold,&membuffer[I_q],fluxes->num_bytes);
        }        
    }

    /* ------------------------------ Clean up -----------------------------------------*/ 
    cudaFree(array_fluxes_struct_dev);
    cudaFree(membuffer_dev);
    cudaFreeHost(membuffer);

    return maxcfl;
}

