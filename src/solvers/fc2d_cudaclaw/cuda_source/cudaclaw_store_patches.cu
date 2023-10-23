#include "../fc2d_cudaclaw.h"
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fc2d_cudaclaw_options.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>

void cudaclaw_store_buffer(fclaw2d_global_t* glob,
                           fclaw2d_patch_t *this_patch,
                           int this_patch_idx,
                           int this_block_idx,
                           int count, int iter, 
                           cudaclaw_fluxes_t* flux_array,
                           fclaw2d_patch_t** patch_array,
                           int* patchno_array,
                           int* blockno_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);
    double *qold, *aux;
    int meqn, maux;

    const fc2d_cudaclaw_options_t *cuclaw_opt = fc2d_cudaclaw_get_options(glob);


    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    fluxes->qold = qold;
    fluxes->aux = aux;

    int n = iter % cuclaw_opt->buffer_len;
    flux_array[n] = *fluxes;
    if (cuclaw_opt->src_term > 0)
    {
        patch_array[n] = this_patch;
        patchno_array[n] = this_patch_idx;
        blockno_array[n] = this_block_idx;
    }
        
}
