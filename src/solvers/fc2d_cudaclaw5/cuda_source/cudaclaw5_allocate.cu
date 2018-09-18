#include "cudaclaw5_allocate.h"

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw_timer.h>



void cudaclaw5_allocate_fluxes(struct fclaw2d_global *glob,
                               struct fclaw2d_patch *patch)
{
    size_t size = fclaw2d_clawpatch_size(glob);

    cudaclaw5_fluxes_t *fluxes = FCLAW_ALLOC(cudaclaw5_fluxes,1);
    fluxes->num_bytes = size*sizeof(double);

    /* CPU memory allocation */
    fluxes->fp = new double[size];
    fluxes->fm = new double[size];
    fluxes->gp = new double[size];
    fluxes->gm = new double[size];


    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]);
    cudaMalloc((void**)&fluxes->qold_dev, size * sizeof(double));
    cudaMalloc((void**)&fluxes->fm_dev, size * sizeof(double));
    cudaMalloc((void**)&fluxes->fp_dev, size * sizeof(double));
    cudaMalloc((void**)&fluxes->gm_dev, size * sizeof(double));
    cudaMalloc((void**)&fluxes->gp_dev, size * sizeof(double));
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_CUDA_ALLOCATE]);

    fclaw2d_patch_set_user_data(glob,patch,fluxes);
}
