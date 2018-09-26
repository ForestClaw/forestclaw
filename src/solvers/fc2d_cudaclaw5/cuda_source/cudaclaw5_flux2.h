#include <fc2d_cudaclaw5.h>
__global__ void cudaclaw5_flux2(int idir, int mx, int my,int meqn, int mbc, int maux,
				double* qold, double* aux, double dt, double dx, double dy,
				double* cflgrid, double* fm, double* fp,
				double* gm, double* gp, cudaclaw5_cuda_rpn2_t rpn2, void* rpt2,
				int mwaves);
