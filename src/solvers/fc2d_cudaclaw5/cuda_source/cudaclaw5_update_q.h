#ifndef CUDACLAW5_UPDATE_Q_H
#define CUDACLAW5_UPDATE_Q_H


/* Only include headers needed to get this file to compile;  all other
   headers should go in c files */

#include "../fc2d_cudaclaw5.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

__global__  void cudaclaw5_update_q_cuda(int mbc,
                             double dtdx, double dtdy,
                             double* qold,
                             double* fm, double* fp,
                             double* gm, double* gp);

__global__ void cudaclaw5_update_q_cuda2(int mbc, int mx, int my, int meqn,
                                        double dtdx, double dtdy,
                                        double* qold, 
                                        double* fm, double* fp, 
                                        double* gm, double* gp);


void cudaclaw5_update_q(int meqn, int mx, int my, int mbc,
                         double dtdx, double dtdy,
                         double qold[],
                         double fm[], double fp[],
                         double gm[], double gp[], int mcapa);
#ifdef __cplusplus
#if 0
{
#endif
}
#endif
#endif
