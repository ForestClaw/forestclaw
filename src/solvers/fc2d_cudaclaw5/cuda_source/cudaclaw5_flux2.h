#ifndef CUDACLAW5_FLUX2_H
#define CUDACLAW5_FLUX2_H


/* Only include headers needed to get this file to compile;  all other
   headers should go in c files */

#ifdef __cplusplus
extern "C"
{
#endif

typedef void (*cudaclaw5_cuda_rpn2_t)(int idir, int meqn, int mwaves, int maux,
                                      double ql[], double qr[], 
                                      double auxl[], double auxr[],
                                      double wave[], double s[], 
                                      double amdq[], double apdq[]);



__global__ void cudaclaw5_flux2(int idir, int mx, int my, int meqn, int mbc,
                                int maux, int mwaves, 
                                double* qold, double* aux, double dx,
                                double dy, double dt, double* cflgrid,
                                double* fm, double* fp, double* gm, double* gp,
                                double* waves, double *speeds,
                                cudaclaw5_cuda_rpn2_t rpn2, void* rpt2);

__global__ void cudaclaw5_compute_cfl(int idir, int mx, int my, int meqn, int mwaves, 
                                     int mbc, double dx, double dy, double dt, 
                                     double *speeds, double* cflgrid);




#ifdef __cplusplus
}
#endif
#endif

