/* # check to see if value exceeds threshold */

#include "fclaw_clawpatch.h"

int fclaw2d_clawpatch_exceeds_threshold(int* blockno,
                                        double qval[], 
                                        double* qmin, double *qmax,
                                        double quad[], 
                                        double *dx, double *dy, 
                                        double *xc, double *yc
                                        int* tag_threshold)
{

}
    int refine;

    int refine_criteria = fclaw2d_clawpatch_get_refinement_criteria();

    logical :: exceeds_th
    logical :: fclaw2d_clawpatch_gradient_exceeds_th
    logical :: fclaw2d_clawpatch_difference_exceeds_th
    logical :: fclaw2d_clawpatch_value_exceeds_th
    logical :: fclaw2d_clawpatch_minmax_exceeds_th

    if (refine_criteria .eq. 0) then 
        exceeds_th = fclaw2d_clawpatch_value_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 1) then
        exceeds_th = fclaw2d_clawpatch_minmax_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 2) then
        exceeds_th = fclaw2d_clawpatch_difference_exceeds_th(blockno, &
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    else if (refine_criteria .eq. 3) then
        exceeds_th = fclaw2d_clawpatch_gradient_exceeds_th(blockno, & 
                    qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                    tag_threshold)
    end if

    transport_exceeds_threshold = exceeds_th

end function transport_exceeds_threshold
