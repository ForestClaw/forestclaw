/* # check to see if value exceeds threshold */

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif


#include<fclaw2d_global.h>

#if REFINE_DIM == 2 && PATCH_DIM == 2

#include "../fclaw_clawpatch.h"

#include "../fclaw_clawpatch_options.h"
#include "../fclaw2d_clawpatch_fort.h"

#elif REFINE_DIM == 2 && PATCH_DIM == 3

#include "../fclaw3dx_clawpatch.h"

#include "../fclaw_clawpatch_options.h"
#include "../fclaw2d_clawpatch_fort.h"
#include "../fclaw3dx_clawpatch_fort.h"

#include <_fclaw2d_to_fclaw3dx.h>

#endif


/* ------------------------------------------------------------------------------------ */

/* This signature is not the same as the actual threshold functions called below 
   This is because we can call the glob and get options. This is only called from 
   the tagging/coarsening routines in patches/clawpatch.   If users copy these 
   routines, they can then alter the call to this function in any way they wish.
*/
int FCLAW2D_CLAWPATCH_TAG_CRITERIA(const int* blockno,
                                        const double qval[], 
                                        const double qmin[], 
                                        const double qmax[],
                                        const double quad[], 
                                        const double *dx, 
                                        const double *dy, 
                                        const double *xc, 
                                        const double *yc, 
                                        const double *tag_threshold,
                                        const int* init_flag,
                                        const int* is_ghost)
{
    struct fclaw_global* glob = fclaw_global_get_static_global();
    fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

    clawpatch_fort_exceeds_threshold_t user_exceeds_threshold = 
                                clawpatch_vt->d2->fort_user_exceeds_threshold;

    fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int meqn_val = clawpatch_opt->meqn, *meqn = &meqn_val;
    int ivar_val = clawpatch_opt->threshold_variable, *ivar_threshold=&ivar_val;

    int exceeds_th = 1;
    int refinement_criteria = fclaw_clawpatch_get_options(glob)->refinement_criteria;
    switch(refinement_criteria)
    {
        case FCLAW_REFINE_CRITERIA_VALUE:
            exceeds_th = FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    ivar_threshold,
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_DIFFERENCE:
            exceeds_th = FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    ivar_threshold,
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_MINMAX:
            exceeds_th = FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    ivar_threshold,
                    tag_threshold,init_flag,is_ghost);

            break;
        case FCLAW_REFINE_CRITERIA_GRADIENT:
            exceeds_th = FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH(blockno,meqn,
                    qval,qmin,qmax,quad, dx, dy, xc, yc,     
                    ivar_threshold,
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
                                                dx, dy, xc, yc,
                                                ivar_threshold,
                                                tag_threshold,init_flag,is_ghost);
            break;
        default:
            fclaw_global_essentialf("fclaw2d_clawpatch_exceeds_threshold.c): " \
                                    "No valid refinement criteria specified\n");
            exit(0);
            break;
    }
    return exceeds_th;
}
