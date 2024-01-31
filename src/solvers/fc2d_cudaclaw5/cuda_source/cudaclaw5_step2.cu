#include "../fc2d_cudaclaw5.h"
#include "cudaclaw5_allocate.h"

#include <fclaw_patch.h>
#include <fclaw_global.h>
#include <fclaw2d_clawpatch.h>

#include "../fc2d_cudaclaw5_options.h"

#include "cudaclaw5_update_q.h"
#include "cudaclaw5_flux2.h"

// Extra Source
#include <cub/cub.cuh>

#if 0
__device__ void rpn2adv_cuda(int idir, int meqn, int mwaves, int maux,
     double ql[], double qr[], double auxl[], double auxr[],
     double wave[], double s[], double amdq[], double apdq[])
{
    /* wave[mwaves][meqn] */
    /* idir in 0,1 : needed to get correct  */

    wave[0] = qr[0] - ql[0];
    s[0] = auxr[idir];
    amdq[0] = SC_MIN(auxr[idir], 0) * wave[0];
    apdq[0] = SC_MAX(auxr[idir], 0) * wave[0];
}
#endif

//__device__ cudaclaw5_cuda_rpn2_t rpn2_dev =  rpn2adv_cuda;

// CustomMax functor
struct AbsMax
{
    template <typename T>
    __device__ __forceinline__
    T operator()(const T &a, const T &b) const {
        return (abs(b) > abs(a)) ? abs(b) : abs(a);
    }
};

static
double cudaclaw5_compute_cfl(int idir, int mx, int my, int mbc, int mwaves, 
			     double dx, double dy, double dt, double* speeds_dev)
{
	int size = (2*mbc+mx)*(2*mbc+my)*mwaves;
	void    *temp_storage_dev = NULL;
	size_t  temp_storage_bytes = 0;
	double  *cflgrid_dev;
	double cflgrid;

	AbsMax maxop;
	double init=0.0;
	cudaMalloc(&cflgrid_dev, sizeof(double));	
	CubDebugExit(cub::DeviceReduce::Reduce(temp_storage_dev,temp_storage_bytes,
                                         speeds_dev,cflgrid_dev,size,maxop,init));
	cudaMalloc(&temp_storage_dev, temp_storage_bytes);  
	CubDebugExit(cub::DeviceReduce::Reduce(temp_storage_dev,temp_storage_bytes,
                                       speeds_dev,cflgrid_dev,size,maxop,init));
	cudaMemcpy(&cflgrid, cflgrid_dev, sizeof(double),cudaMemcpyDeviceToHost);

	cudaFree(temp_storage_dev);
	cudaFree(cflgrid_dev);
	if(idir == 0)
	{
		cflgrid = cflgrid*dt/dx;
	}
	else
	{
		cflgrid =  cflgrid*dt/dy;
	}

	return cflgrid;
}

double cudaclaw5_step2(fclaw_global_t *glob,
                       fclaw_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds;
    cudaError_t code;
    
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();

    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;
    double cflgrid;

    fc2d_cudaclaw5_options_t* cuda_opt = fc2d_cudaclaw5_get_options(glob);

    cudaclaw5_fluxes_t *fluxes = (cudaclaw5_fluxes_t*) 
               fclaw_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    FCLAW_ASSERT(cuclaw5_vt->fort_rpn2 != NULL);
    FCLAW_ASSERT(cuclaw5_vt->fort_rpt2 != NULL);

    FCLAW_ASSERT(cuclaw5_vt->cuda_rpn2 != NULL);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    fclaw2d_clawpatch_save_current_step(glob, this_patch);
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

#if 0
    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];
#endif    


    int ierror = 0;
    // cudaclaw5_fort_flux2_t flux2 = CUDACLAW5_FLUX2;

    int* block_corner_count = fclaw_patch_block_corner_count(glob,this_patch);

    size_t size = fclaw2d_clawpatch_size(glob);

    /* -------------------------- Construct fluctuations -------------------------------*/ 
    cudaEventRecord(start);

    cudaMemcpy(fluxes->qold_dev, qold, fluxes->num_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(fluxes->aux_dev, aux, fluxes->num_bytes_aux, cudaMemcpyHostToDevice);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop); 
    glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY].cumulative += milliseconds*1e-3;

    {
        dim3 block(32,32);  
        dim3 grid((mx+2*mbc-1+block.x-1)/block.x,(my+2*(mbc-1)+block.y-1)/block.y);

        int mwaves = cuda_opt->mwaves;
        cflgrid = 0.0;

        cudaEventRecord(start);

        /* ---------------------------------------------------------------------------- */
        /* X direction */
        /* ---------------------------------------------------------------------------- */
        cudaclaw5_flux2<<<grid, block>>>(0,mx,my,meqn,mbc,maux,mwaves, 
                                         fluxes->qold_dev,
                                         fluxes->aux_dev, dx,dy,dt,&cflgrid,
                                         fluxes->fm_dev,fluxes->fp_dev,
                                         fluxes->gm_dev,fluxes->gp_dev,
                                         fluxes->waves_dev, fluxes->speeds_dev,
                                         cuclaw5_vt->cuda_rpn2, NULL);

        code = cudaPeekAtLastError();
        if(code != cudaSuccess)
        {
            printf("ERROR (1): %s\n",cudaGetErrorString(code));
            exit(code);
        }

        cflgrid = cudaclaw5_compute_cfl(0,mx,my,mbc,mwaves,dx,dy,dt,fluxes->speeds_dev);
        //cudaclaw5_compute_cfl<<<grid, block>>>(0,mx,my,meqn,mwaves, mbc,
        //                                       dx,dy,dt,fluxes->speeds_dev, &cflgrid);
        cudaDeviceSynchronize();

        /* ---------------------------------------------------------------------------- */
        /* Y direction */
        /* ---------------------------------------------------------------------------- */
        cudaclaw5_flux2<<<grid, block>>>(1,mx,my,meqn,mbc,maux,mwaves,
                                         fluxes->qold_dev,fluxes->aux_dev, 
                                         dx,dy,dt,&cflgrid,
                                         fluxes->fm_dev,fluxes->fp_dev,
                                         fluxes->gm_dev,fluxes->gp_dev,
                                         fluxes->waves_dev, fluxes->speeds_dev,
                                         cuclaw5_vt->cuda_rpn2,NULL);
        code = cudaPeekAtLastError();
        if(code != cudaSuccess)
        {
            printf("ERROR (2): %s\n",cudaGetErrorString(code));
            exit(code);
        }


        cflgrid = cudaclaw5_compute_cfl(0,mx,my,mbc,mwaves,dx,dy,dt,fluxes->speeds_dev);
        cudaDeviceSynchronize();

#if 0
        cudaclaw5_compute_cfl<<<grid, block>>>(1,mx,my,meqn,mwaves, mbc,
                                               dx,dy,dt,fluxes->speeds_dev, &cflgrid);
#endif                                               

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        glob->timers[FCLAW2D_TIMER_CUDA_KERNEL1].cumulative += milliseconds*1e-3;        
    }


    /* -------------------------- Update solution --------------------------------------*/ 
    double dtdx, dtdy;
    dtdx = dt/dx;
    dtdy = dt/dy;

    cudaEventRecord(start);

    dim3 block(32,32);  
    dim3 grid((mx+block.x-1)/block.x,(my+block.y-1)/block.y);

    cudaclaw5_update_q_cuda2<<<grid, block>>>(mbc, mx,my,meqn,dtdx, dtdy, 
                                              fluxes->qold_dev, 
                                              fluxes->fm_dev, fluxes->fp_dev,
                                              fluxes->gm_dev, fluxes->gp_dev);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    glob->timers[FCLAW2D_TIMER_CUDA_KERNEL2].cumulative += milliseconds*1e-3;

    code = cudaPeekAtLastError();
    if(code != cudaSuccess)
    {
        printf("ERROR (3): %s\n",cudaGetErrorString(code));
        exit(code);
    }


    /* -------------------------- Copy q back to host ----------------------------------*/ 
    cudaEventRecord(start);
    cudaMemcpy(qold, fluxes->qold_dev, fluxes->num_bytes, cudaMemcpyDeviceToHost);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    glob->timers[FCLAW2D_TIMER_CUDA_MEMCOPY].cumulative += milliseconds*1e-3;
    
    /* ------------------------------ Clean up -----------------------------------------*/ 
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    FCLAW_ASSERT(ierror == 0);

    return cflgrid;
}

#if 0
/* Use for possible work arrays */
c     # Local variables
      integer i0faddm, i0faddp, i0gaddm, i0gaddp
      integer i0q1d, i0dtdx1, i0dtdy1
      integer i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl
      integer i0auxr

      integer i,j,m

c     Needed by Riemann solvers.  This should be fixed later by a 'context'
c     for a Riemann solver.
      double precision dtcom, dxcom,dycom,tcom
      integer icom, jcom
      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcomdouble dtdx, double dtdy,
                            double* qold,
                            double* fm, double* fp,
                            double* gm, double* gp);

c     # This should be set to actual time, in case the user wants it
c     # it for some reason in the Riemann solver.

c     # Set up work arrays (these are not used yet)

      i0faddm = 1
      i0faddp = i0faddm +   (maxm+2*mbc)*meqn
      i0gaddm = i0faddp +   (maxm+2*mbc)*meqn
      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d   +   (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1 +   (maxm+2*mbc)
      i0aux1  = i0dtdy1 +   (maxm+2*mbc)
      i0aux2  = i0aux1  +   (maxm+2*mbc)*maux
      i0aux3  = i0aux2  +   (maxm+2*mbc)*maux
c
c
      i0next  = i0aux3  + (maxm+2*mbc)*maux    !# next free space
      mused   = i0next - 1                    !# space already used
      mwork1  = mwork - mused           !# remaining space (passed to step2)

      if (mused.gt.mwork) then
         ierror = 1
         return
      endifid need for
c     # global array
c      call cudaclaw5_step2(maxm,maxmx,maxmy,meqn,maux, mbc,
c     &      mx,my, qold,aux,dx,dy,dt,
c     &      cfl,fm,fp,gm,gp,
c     &      work(i0faddm),work(i0faddp),
c     &      work(i0gaddm),work(i0gaddp),
c     &      work(i0q1d),work(i0dtdx1),work(i0dtdy1),
c     &      work(i0aux1),work(i0aux2),work(i0aux3),
c     &      work(i0next),mwork1,rpn2,rpt2,flux2,
c     &      mwaves,mcapa,method,mthlim,block_corner_count,ierror)
#endif

