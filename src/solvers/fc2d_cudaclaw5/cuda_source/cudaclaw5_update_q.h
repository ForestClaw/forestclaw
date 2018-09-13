#ifndef CUDACLAW5_UPDATE_Q_H
#define CUDACLAW5_UPDATE_Q_H


#include "../fc2d_cudaclaw5.h"

#include "../fc2d_cudaclaw5_fort.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch.hpp>

#include <fclaw2d_clawpatch_output_ascii.h>
#include <fclaw2d_clawpatch_output_vtk.h>


#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

__global__
void cudaclaw5_update_q_cuda(int x_stride, int mbc,
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
