#ifndef CUDACLAW5_STEP2_H
#define CUDACLAW5_STEP2_H
#include "../fc2d_cudaclaw5.h"
#include "../fc2d_cudaclaw5_fort.h"
#include "../fc2d_cudaclaw5_options.h"

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
void update_q_(int& meqn,int& mx, int&my,int& mbc,double& dtdx,double& dtdy,double*
        qold,double*fm,double*fp,double*gm,double*gp,int&mcapa);
#ifdef __cplusplus
#if 0
{
#endif
}
#endif
#endif
