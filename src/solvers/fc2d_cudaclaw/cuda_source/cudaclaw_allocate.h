#ifndef CUDACLAW_ALLOCATE_H
#define CUDACLAW_ALLOCATE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdlib.h>   /* defines size_t */


struct fclaw2d_patch;
struct fclaw2d_global;

/* This can breaks cuda if memory is not aligne dproperly; use -malign-double flag
   in gcc */
   
typedef struct  cudaclaw_fluxes
{
    size_t num;
    size_t num_aux;
    size_t num_waves;
    size_t num_speeds;

    size_t num_bytes;   
    size_t num_bytes_aux;  
    size_t num_bytes_waves;  
    size_t num_bytes_speeds;  

    double* qold;
    double* aux;
    double *qold_dev;
    double *aux_dev;
    
    double *fp_dev;
    double *fm_dev;
    double *gp_dev;
    double *gm_dev;

    double *apdq_dev;
    double *amdq_dev;
    double *bpdq_dev;
    double *bmdq_dev;

    double *waves_dev;
    double *speeds_dev;

    double dx;
    double dy;

    double xlower;
    double ylower;
} cudaclaw_fluxes_t;

typedef struct cudaclaw_patch_data
{
    struct fclaw2d_patch **patch_array;
    cudaclaw_fluxes_t *flux_array;
    int *blockno_array;
    int *patchno_array;
}cudaclaw_patch_data_t;


#ifdef __cplusplus
}
#endif
#endif
