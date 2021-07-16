/* # check to see if value exceeds threshold */

#include <fclaw2d_clawpatch.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>

int FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD(int* blockno,
                                        double qval[], 
                                        double* qmin, double *qmax,
                                        double quad[], 
                                        double *dx, double *dy, 
                                        double *xc, double *yc, 
                                        int* tag_threshold,
                                        int* init_flag,
                                        int* is_ghost)
{
#if 0

    int refinement_criteria = fclaw2d_clawpatch_get_refinement_criteria();

    int exceeds_th;

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt();
    clawpatch_fort_exceeds_threshold_t user_exceeds_threshold = 
                                clawpatch_vt->fort_user_exceeds_threshold;

    exceeds_th = 1;
    switch(refinement_criteria)
    {
        case FCLAW_REFINE_CRITERIA_VALUE:
            exceeds_th = FCLAW2D_CLAWPATCH_VALUE_EXCEEDS_TH(blockno,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_DIFFERENCE:
            exceeds_th = FCLAW2D_CLAWPATCH_DIFFERENCE_EXCEEDS_TH(blockno,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    tag_threshold,init_flag,is_ghost);
            break;
        case FCLAW_REFINE_CRITERIA_MINMAX:
            exceeds_th = FCLAW2D_CLAWPATCH_MINMAX_EXCEEDS_TH(blockno,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
                    tag_threshold,init_flag,is_ghost);

            break;
        case FCLAW_REFINE_CRITERIA_GRADIENT:
            exceeds_th = FCLAW2D_CLAWPATCH_GRADIENT_EXCEEDS_TH(blockno,
                    qval,qmin,qmax,quad, dx, dy, xc, yc, 
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
            exceeds_th = user_exceeds_threshold(blockno, qval,qmin,qmax,quad, 
                                                dx, dy, xc, yc,
                                                tag_threshold,init_flag,is_ghost);
            break;
        default:
            fclaw_global_essentialf("fclaw2d_clawpatch_exceeds_threshold.c): " \
                                    "No valid refinement criteria specified\n");
            exit(0);
            break;
    }
    return exceeds_th;
#endif
}
