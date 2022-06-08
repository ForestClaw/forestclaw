#include "../fc2d_cudaclaw.h"

#include "../fc2d_cudaclaw_cuda.h"

#include "cudaclaw_allocate.h"  /* Needed for def of cudaclaw_fluxes_t */


#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fc2d_cudaclaw_options.h>

#include "../fc2d_cudaclaw_check.h"  /* CHECK defined here */

#include <fc2d_cuda_profiler.h>
#include <cub/cub.cuh>

#define thread_count 224


/* Put header here so it doesn't have to go in *.h file */
__global__
void cudaclaw_flux2_and_update_batch (const int mx,    const int my, 
                                      const int meqn,  const int mbc, 
                                      const int maux,  const int mwaves, 
                                      const int mwork,
                                      const double dt, const double t,
                                      struct cudaclaw_fluxes* array_fluxes_struct_dev,
                                      double * maxcflblocks,
                                      cudaclaw_cuda_rpn2_t rpn2,
                                      cudaclaw_cuda_rpt2_t rpt2,
                                      cudaclaw_cuda_b4step2_t b4step2);

__global__
void cudaclaw_compute_speeds_batch (const int mx,    const int my, 
                                    const int meqn,  const int mbc, 
                                    const int maux,  const int mwaves, 
                                    const int mwork,
                                    const double dt, const double t,
                                    cudaclaw_fluxes_t* array_fluxes_struct,
                                    double * maxcflblocks,
                                    cudaclaw_cuda_speeds_t compute_speeds,
                                    cudaclaw_cuda_b4step2_t b4step2);


double cudaclaw_step2_batch(fclaw2d_global_t *glob,
        cudaclaw_fluxes_t* array_fluxes_struct, 
        int batch_size, double t, double dt)
{
    PROFILE_CUDA_GROUP("cudaclaw_step2_batch",1);

    size_t size, bytes, bytes_per_thread;
    float bytes_kb;
    int I_q, I_aux, mwork;
    int i;

    int mx,my,mbc,maux,meqn,mwaves;
    double maxcfl;

    double* maxcflblocks_dev;    

    double *membuffer_cpu, *membuffer_dev;
    cudaclaw_fluxes_t *array_fluxes_struct_dev;

    cudaclaw_fluxes_t* fluxes;

    /* To get patch-independent parameters */
    fc2d_cudaclaw_options_t *clawopt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    /* ---------------------------------- start code ---------------------------------- */
    FCLAW_ASSERT(batch_size > 0);

    clawopt = fc2d_cudaclaw_get_options(glob);
    mwaves = clawopt->mwaves;

    fc2d_cudaclaw_vtable_t*  cuclaw_vt = fc2d_cudaclaw_vt(glob);
    FCLAW_ASSERT(cuclaw_vt->cuda_rpn2 != NULL);
    if (clawopt->order[1] > 0)
    {
        FCLAW_ASSERT(cuclaw_vt->cuda_rpt2 != NULL);        
    }

    clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;
    mbc = clawpatch_opt->mbc;
    maux = clawpatch_opt->maux;
    meqn = clawpatch_opt->meqn;  

    fluxes = &(array_fluxes_struct[0]);
    size = batch_size*(fluxes->num + fluxes->num_aux);
    bytes = size*sizeof(double);

    /* ---------------------------------- Merge Memory ---------------------------------*/ 
    membuffer_cpu = cudaclaw_get_cpu_membuffer();
    membuffer_dev = cudaclaw_get_gpu_membuffer();
    fclaw2d_timer_start_threadsafe (&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2H]);       
    {
        PROFILE_CUDA_GROUP("Copy data on patches to CPU memory buffer",5);    
        for(i = 0; i < batch_size; i++)   
        {
            fluxes = &(array_fluxes_struct[i]);    

            I_q = i*fluxes->num;
            memcpy(&membuffer_cpu[I_q]  ,fluxes->qold ,fluxes->num_bytes);
            fluxes->qold_dev = &membuffer_dev[I_q];

            if (fluxes->num_aux > 0)
            {
                I_aux = batch_size*fluxes->num + i*fluxes->num_aux;
                memcpy(&membuffer_cpu[I_aux],fluxes->aux  ,fluxes->num_bytes_aux);                
                fluxes->aux_dev  = &membuffer_dev[I_aux];
            }
        }  
    }     
    fclaw2d_timer_stop_threadsafe(&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2H]);       
  
    

    fclaw2d_timer_start_threadsafe(&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2D]);       
    {
        PROFILE_CUDA_GROUP("Copy CPU buffer to device memory",3);              
        CHECK(cudaMemcpy(membuffer_dev, membuffer_cpu, bytes, cudaMemcpyHostToDevice));            
    }            
    fclaw2d_timer_stop_threadsafe(&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2D]);       


    /* -------------------------------- Work with array --------------------------------*/ 

    {
        PROFILE_CUDA_GROUP("Copy fluxes to device memory",3);    

        array_fluxes_struct_dev = cudaclaw_get_flux_buffer();

        CHECK(cudaMemcpy(array_fluxes_struct_dev, array_fluxes_struct, 
                         batch_size*sizeof(cudaclaw_fluxes_t), 
                         cudaMemcpyHostToDevice));
    }        

    {
        PROFILE_CUDA_GROUP("Malloc for CFL computation",2);    

        /* Data needed to reduce CFL number */
        CHECK(cudaMalloc(&maxcflblocks_dev,batch_size*sizeof(double)));         
    }


#if 0
    {
        PROFILE_CUDA_GROUP("Configure and call to compute speeds",6);  

        /* Compute speeds */
        int block_size = FC2D_CUDACLAW_BLOCK_SIZE;

        dim3 block(block_size,1,1);
        dim3 grid(1,1,batch_size);

        mwork = 2*(meqn + maux) + mwaves;
        bytes_per_thread = sizeof(double)*mwork;
        bytes = bytes_per_thread*block_size;
        bytes_kb = bytes/1024.0;

        cudaclaw_compute_speeds_batch <<<grid,block,bytes>>>(mx,my,meqn, mbc, maux, mwaves,
                                                             mwork, dt, t, 
                                                             array_fluxes_struct_dev,
                                                             maxcflblocks_dev,
                                                             cuclaw_vt->cuda_speeds,
                                                             cuclaw_vt->cuda_b4step2);

        cudaDeviceSynchronize();

        cudaError_t code = cudaPeekAtLastError();
        if (code != cudaSuccess) 
        {
            fclaw_global_essentialf("ERROR (cudaclaw_step2.cu (compute_speeds)) : %s\n\n", 
                                    cudaGetErrorString(code));
            exit(code);
        }
    }        
#endif    

    {
        PROFILE_CUDA_GROUP("Configure and call main kernel",6);  

        /* Determine shared memory size */
        int block_size = FC2D_CUDACLAW_BLOCK_SIZE;
        //int block_size = thread_count;
        dim3 block(block_size,1,1);
        dim3 grid(1,1,batch_size);

        int mwork1 = 4*meqn + 2*maux + mwaves + meqn*mwaves;
        int mwork2 = 5*meqn + 6*maux;
        mwork = (mwork1 > mwork2) ? mwork1 : mwork2;
        bytes_per_thread = sizeof(double)*mwork;
        bytes = bytes_per_thread*block_size;

        bytes_kb = bytes/1024.0;
        //fclaw_global_essentialf("[fclaw] Shared memory  : %0.2f kb\n\n",bytes_kb);

        cudaclaw_flux2_and_update_batch<<<grid,block,bytes>>>(mx,my,meqn,mbc,maux,mwaves,
                                                              mwork, dt,t,
                                                              array_fluxes_struct_dev,
                                                              maxcflblocks_dev,
                                                              cuclaw_vt->cuda_rpn2,
                                                              cuclaw_vt->cuda_rpt2,
                                                              cuclaw_vt->cuda_b4step2);
        cudaDeviceSynchronize();

        
        cudaError_t code = cudaPeekAtLastError();

        if (code != cudaSuccess) 
        {
            fclaw_global_essentialf("ERROR (cudaclaw_step2.cu) : %s\n", 
                                    cudaGetErrorString(code));
            exit(code);
        }
    }

    /* -------------------------------- Finish CFL ------------------------------------*/ 
    {
        PROFILE_CUDA_GROUP("Finish CFL",2);
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


	
    /* -------------------------- Copy q back to host ----------------------------------*/ 
    fclaw2d_timer_start_threadsafe (&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_D2H]);       

    {
        PROFILE_CUDA_GROUP("Copy device memory buffer back to CPU",3);

        CHECK(cudaMemcpy(membuffer_cpu, membuffer_dev, batch_size*fluxes->num_bytes, 
                         cudaMemcpyDeviceToHost));
    }
    fclaw2d_timer_stop_threadsafe (&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_D2H]);       

    fclaw2d_timer_start_threadsafe (&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2H]);       
    {
        PROFILE_CUDA_GROUP("Copy CPU buffer back to patches",5);
        for (i = 0; i < batch_size; ++i)    
        {      
            fluxes = &(array_fluxes_struct[i]);
            I_q = i*fluxes->num;

            memcpy(fluxes->qold,&membuffer_cpu[I_q],fluxes->num_bytes);
        }        
    }
    fclaw2d_timer_stop_threadsafe (&glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY_H2H]);       

    return maxcfl;
}

