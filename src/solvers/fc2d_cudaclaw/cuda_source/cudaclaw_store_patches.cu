#include "../fc2d_cudaclaw.h"
#include <fclaw_global.h>
#include <fclaw_patch.h>

#include <fclaw_clawpatch.h>
#include <fc2d_cudaclaw_options.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>

void cudaclaw_store_buffer(fclaw_global_t* glob,
                           fclaw_patch_t *this_patch,
                           int this_patch_idx,
                           int count, int iter, 
                           cudaclaw_fluxes_t* flux_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);
    double *qold, *aux;
    int meqn, maux;

    const fc2d_cudaclaw_options_t *cuclaw_opt = fc2d_cudaclaw_get_options(glob);


    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    fclaw_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    fclaw_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    fluxes->qold = qold;
    fluxes->aux = aux;

    flux_array[iter % cuclaw_opt->buffer_len] = *fluxes;
}
