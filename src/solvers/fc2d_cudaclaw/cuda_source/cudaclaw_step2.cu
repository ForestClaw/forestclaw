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

static double* s_membuffer;
static double* s_membuffer_dev;

cudaclaw_fluxes_t* s_array_fluxes_struct_dev;


void cudaclaw_allocate_buffers(fclaw2d_global_t *glob)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    int maux = clawpatch_opt->maux;
    int meqn = clawpatch_opt->meqn;  

    size_t size = (2*mbc+mx)*(2*mbc+my);
    size_t bytes = size*(meqn + maux)*sizeof(double);

    CHECK(cudaMallocHost((void**)&s_membuffer,bytes));    
    CHECK(cudaMalloc((void**)&s_membuffer_dev, bytes)); 

    int batch_size = FC2D_CUDACLAW_BUFFER_LEN;
    CHECK(cudaMalloc(&s_array_fluxes_struct_dev, batch_size*sizeof(cudaclaw_fluxes_t)));
}

void cudaclaw_deallocate_buffers(fclaw2d_global_t *glob)
{
    cudaFreeHost(s_membuffer);
    cudaFree(s_membuffer_dev);
    cudaFree(s_array_fluxes_struct_dev);
}


double cudaclaw_step2_batch(fclaw2d_global_t *glob,
        cudaclaw_fluxes_t* array_fluxes_struct, 
        int batch_size, double t, double dt)
{
    PROFILE_CUDA_GROUP("cudaclaw_step2_batch",5);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    //float milliseconds;
    int i;

    double maxcfl = 0.0;

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
    //double *membuffer;
    //double* membuffer_dev;

    /* ---------------------------------- Merge Memory ---------------------------------*/ 
    {
        {
            PROFILE_CUDA_GROUP("Malloc buffer on the host and device",2);    
            FCLAW_ASSERT(s_membuffer != NULL);
            FCLAW_ASSERT(s_membuffer_dev != NULL);
            //CHECK(cudaMallocHost((void**)&membuffer,bytes));

            //CHECK(cudaMalloc((void**)&membuffer_dev, bytes));            
        }

        PROFILE_CUDA_GROUP("cudaclaw_copy_loop",3);    
        for(i = 0; i < batch_size; i++)   
        {
            cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);    

            int I_q = i*fluxes->num;
            int I_aux = batch_size*fluxes->num + i*fluxes->num_aux;

            printf("here ...0\n");
            memcpy(&s_membuffer[I_q]  ,fluxes->qold ,fluxes->num_bytes);
            printf("here ...1\n");
            memcpy(&s_membuffer[I_aux],fluxes->aux  ,fluxes->num_bytes_aux);

            /* Assign gpu pointers */
            fluxes->qold_dev = &s_membuffer_dev[I_q];
            fluxes->aux_dev  = &s_membuffer_dev[I_aux];
        }        

        {
            PROFILE_CUDA_GROUP("Copy buffer to device",7);              
            printf("here ...2\n");
            CHECK(cudaMemcpy(s_membuffer_dev, s_membuffer, bytes, cudaMemcpyHostToDevice));            
        }
    }        

    printf("here ...3\n");

    /* -------------------------------- Work with array --------------------------------*/ 

    //cudaclaw_fluxes_t* array_fluxes_struct_dev = NULL;
    //CHECK(cudaMalloc(&array_fluxes_struct_dev, batch_size*sizeof(cudaclaw_fluxes_t)));

    FCLAW_ASSERT(s_array_fluxes_struct_dev != NULL);

    printf("here ...4\n");
    CHECK(cudaMemcpy(s_array_fluxes_struct_dev, array_fluxes_struct, 
                     batch_size*sizeof(cudaclaw_fluxes_t), 
                     cudaMemcpyHostToDevice));

    printf("here ...5\n");

    dim3 block(128,1,1);
    dim3 grid(1,1,batch_size);

    size_t bytes_per_thread = sizeof(double)*(5*meqn+3*maux+mwaves+meqn*mwaves);
    
    double* maxcflblocks_dev;
    CHECK(cudaMalloc(&maxcflblocks_dev,batch_size*sizeof(double))); 
    cudaclaw_flux2_and_update_batch<<<grid,block,128*bytes_per_thread >>>(mx,my,meqn,
                                                                     mbc,maux,mwaves,dt,t,
                                                                     s_array_fluxes_struct_dev,
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
    CHECK(cudaMemcpy(s_membuffer, s_membuffer_dev, batch_size*fluxes->num_bytes, 
                     cudaMemcpyDeviceToHost));

    {
        PROFILE_CUDA_GROUP("Copy back to patches loop",2);
        for (int i = 0; i < batch_size; ++i)    
        {      

            cudaclaw_fluxes_t* fluxes = &(array_fluxes_struct[i]);
            int I_q = i*fluxes->num;

            printf("here ...5\n");
            memcpy(fluxes->qold,&s_membuffer[I_q],fluxes->num_bytes);
            printf("here ...6\n");

        }        
    }

    /* ------------------------------ Clean up -----------------------------------------*/ 
    //cudaFree(array_fluxes_struct_dev);
    //cudaFree(membuffer_dev);
    //cudaFreeHost(membuffer);

    return maxcfl;
}

