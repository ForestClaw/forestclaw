#include "cudaclaw5_update_q.h"
__global__ void update_q_cuda(int x_stride, int mbc, double dtdx, double dtdy,
                              double* qold, double* fm, double* fp, double* gm,
                              double* gp) {
    int i = (blockIdx.x + mbc) * x_stride + blockIdx.y + mbc;
    qold[i] += -dtdx * (fm[i+x_stride] - fp[i]) -
                          dtdy * (gm[i+1] - gp[i]);
}
void update_q_(int& meqn, int& mx, int& my, int& mbc, double& dtdx,
               double& dtdy, double* qold, double* fm, double* fp, double* gm,
               double* gp, int& mcapa) {
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
    dim3 dimGrid(mx, my);
    dim3 dimBlock(1, 1);

    int x_stride = mx + 2 * mbc;
    update_q_cuda<<<dimGrid, dimBlock>>>(x_stride, mbc, dtdx, dtdy, qold_dev,
                                         fm_dev, fp_dev, gm_dev, gp_dev);

    cudaMemcpy(qold, qold_dev, size * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(qold_dev);
    cudaFree(fm_dev);
    cudaFree(fp_dev);
    cudaFree(gm_dev);
    cudaFree(gp_dev);
}
