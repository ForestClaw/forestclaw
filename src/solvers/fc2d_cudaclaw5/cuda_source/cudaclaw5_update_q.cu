#include "cudaclaw5_update_q.h"

__global__ void cudaclaw5_update_q_cuda(int mbc, 
                                        double dtdx, double dtdy,
                                        double* qold, 
                                        double* fm, double* fp, 
                                        double* gm, double* gp)
{
    int meqn = threadIdx.z;
    int x = threadIdx.x;
    int x_stride = blockDim.z;
    int y = threadIdx.y;
    int y_stride = (blockDim.x+2*mbc)*x_stride;
    int i = meqn + (x+mbc)*x_stride + (y+mbc)*y_stride;
    qold[i] = qold[i] - dtdx * (fm[i+x_stride] - fp[i]) 
                      - dtdy * (gm[i+y_stride] - gp[i]);
}


void cudaclaw5_update_q(int meqn, int mx, int my, int mbc, 
                        double dtdx, double dtdy, double qold[], 
                        double fm[], double fp[], 
                        double gm[], double gp[], int mcapa) 
{
    int size = meqn * (mx + 2 * mbc) * (my + 2 * mbc);

    double* qold_dev;
    double* fm_dev;
    double* fp_dev;
    double* gm_dev;
    double* gp_dev;    

    cudaMalloc((void**)&qold_dev, size * sizeof(double));
    cudaMalloc((void**)&fm_dev, size * sizeof(double));
    cudaMalloc((void**)&fp_dev, size * sizeof(double));
    cudaMalloc((void**)&gm_dev, size * sizeof(double));
    cudaMalloc((void**)&gp_dev, size * sizeof(double));

    cudaMemcpy(qold_dev, qold, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(fm_dev, fm, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(fp_dev, fp, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gm_dev, gm, size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gp_dev, gp, size * sizeof(double), cudaMemcpyHostToDevice);

    // this is not optimal
    dim3 dimGrid(mx, my, meqn);
    dim3 dimBlock(1, 1);
    cudaclaw5_update_q_cuda<<<dimBlock, dimGrid>>>(mbc, dtdx, dtdy, 
                                                   qold_dev, fm_dev, fp_dev, 
                                                   gm_dev, gp_dev);
    //equivalent c loop
    /*
    int x_stride = meqn;
    int y_stride = (mx + 2 * mbc)*x_stride;
    for(int m=0;m<meqn;m++){
        for(int x=0;x<mx;x++){
            for(int y=0;y<my;y++){
                int i = m+(x+mbc)*x_stride+(y+mbc)*y_stride;
                qold[i] =qold[i] -dtdx * (fm[i+x_stride] - fp[i]) -
                                  dtdy * (gm[i+y_stride] - gp[i]);

            }
        }
    }
    */

	cudaError_t code = cudaPeekAtLastError();
    if(code!=cudaSuccess){
        printf("ERROR: %s\n",cudaGetErrorString(code));
    }

    cudaMemcpy(qold, qold_dev, size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(qold_dev);
    cudaFree(fm_dev);
    cudaFree(fp_dev);
    cudaFree(gm_dev);
    cudaFree(gp_dev);
}
