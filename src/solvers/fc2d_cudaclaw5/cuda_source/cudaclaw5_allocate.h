#ifndef CUDACLAW5_ALLOCATE_H
#define CUDACLAW5_ALLOCATE_H


/* Only include headers needed to get this file to compile;  all other
   headers should go in c files */

// #include "../fc2d_cudaclaw5.h"

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw_patch;
struct fclaw_global;

typedef struct cudaclaw5_fluxes
{
    size_t num_bytes;   /* All members have the same size */
    size_t num_bytes_aux;  
    size_t num_bytes_waves;  
    size_t num_bytes_speeds;  

    double *qold_dev;
    double *aux_dev;
    
    double *fp_dev;
    double *fm_dev;
    double *gp_dev;
    double *gm_dev;

    double *waves_dev;
    double *speeds_dev;
    
} cudaclaw5_fluxes_t;

void cudaclaw5_allocate_fluxes(struct fclaw_global *glob,
                               struct fclaw_patch *patch);

void cudaclaw5_deallocate_fluxes(struct fclaw_global *glob,
                                 struct fclaw_patch *patch);


#ifdef __cplusplus
}
#endif
#endif
