/* # check to see if value exceeds threshold */

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 3
#endif

#if 0
#if REFINE_DIM == 2 && PATCH_DIM == 2

#include "../fclaw2d_clawpatch.h"

#include "../fclaw2d_clawpatch_options.h"
#include "../fclaw2d_clawpatch_fort.h"

#elif REFINE_DIM == 2 && PATCH_DIM == 3

#include "../fclaw3dx_clawpatch.h"

#include "../fclaw3dx_clawpatch_options.h"
#include "../fclaw2d_clawpatch_fort.h"
#include "../fclaw3dx_clawpatch_fort.h"

#include <_fclaw2d_to_fclaw3dx.h>

#endif
#endif

#include <fclaw2d_global.h>

#include "../fclaw3dx_clawpatch.h"

#include "../fclaw3dx_clawpatch_options.h"
// #include "../fclaw2d_clawpatch_fort.h"
#include "../fclaw3dx_clawpatch_fort.h"




/* ------------------------------------------------------------------------------------ */

int FCLAW3DX_CLAWPATCH_TAG_CRITERIA(const int* blockno,
                                        const double *qval, 
                                        const double* qmin, 
                                        const double *qmax,
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *dz, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *zc, 
                                        const double *tag_threshold,
                                        const int* init_flag,
                                        const int* is_ghost)
{
    struct fclaw2d_global* glob = fclaw2d_global_get_global();
    fclaw3dx_clawpatch_vtable_t* clawpatch_vt = fclaw3dx_clawpatch_vt(glob);
    fclaw3dx_clawpatch_fort_exceeds_threshold_t user_exceeds_threshold = 
                                clawpatch_vt->fort_user_exceeds_threshold;

    fclaw3dx_clawpatch_options_t *clawpatch_opt = fclaw3dx_clawpatch_get_options(glob);    
    int meqn_val = clawpatch_opt->meqn, *meqn = &meqn_val;
    int ivar_val = clawpatch_opt->threshold_variable, *ivar_variable=&ivar_val;

    int exceeds_th = 1;
    int refinement_criteria = fclaw3dx_clawpatch_get_options(glob)->refinement_criteria;
    switch(refinement_criteria)
    {
        case FCLAW_REFINE_CRITERIA_VALUE:
            exceeds_th = FCLAW3DX_CLAWPATCH_VALUE_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, dz, xc, yc, zc, 
                    ivar_variable,
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_DIFFERENCE:
            exceeds_th = FCLAW3DX_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, dz, xc, yc, zc, 
                    ivar_variable,
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_MINMAX:
            exceeds_th = FCLAW3DX_CLAWPATCH_MINMAX_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, dz, xc, yc, zc, 
                    ivar_variable,
                    tag_threshold,init_flag,is_ghost);

            break;
        case FCLAW_REFINE_CRITERIA_GRADIENT:
            exceeds_th = FCLAW3DX_CLAWPATCH_GRADIENT_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, dz, xc, yc, zc, 
                    ivar_variable,
                    tag_threshold,init_flag,is_ghost);
            break;

        case FCLAW_REFINE_CRITERIA_USER:

            /* This must be define by the user, or nothing happens */
            if (user_exceeds_threshold == NULL)
            {
                fclaw_global_essentialf("\nfclaw2d_clawpatch_exceeds_threshold : " \
                                        "User specified criteria is set, but no user " \
                                        "defined\n function was found.\n\n");
            }
            FCLAW_ASSERT(user_exceeds_threshold != NULL);
            exceeds_th = user_exceeds_threshold(blockno, meqn,qval,qmin,qmax,quad, 
                                                dx, dy, dz, xc, yc, zc, 
                                                ivar_variable,
                                                tag_threshold,init_flag,is_ghost);
            break;
        default:
            fclaw_global_essentialf("fclaw3dx_clawpatch_exceeds_threshold.c): " \
                                    "No valid refinement criteria specified\n");
            exit(0);
            break;
    }
    return exceeds_th;
}
