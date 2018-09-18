#ifndef CUDACLAW5_ALLOCATE_H
#define CUDACLAW5_ALLOCATE_H


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

typedef struct cudaclaw5_fluxes
{
    size_t num_bytes;   /* All members have the same size */
    double *fp;
    double *fm;
    double *gp;
    double *gm;

    double *qold_dev;
    double *fp_dev;
    double *fm_dev;
    double *gp_dev;
    double *gm_dev;
} cudaclaw5_fluxes_t;

void cudaclaw5_allocate_fluxes(struct fclaw2d_global *glob,
                               struct fclaw2d_patch *patch);

void cudaclaw5_deallocate_fluxes(struct fclaw2d_global *glob,
                                 struct fclaw2d_patch *patch);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
#endif
