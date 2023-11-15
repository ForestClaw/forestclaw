/* Five bilinear patches */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
#define SQUARE_BASIS_COMPLETE FCLAW_F77_FUNC(square_basis_complete, \
                            SQUARE_BASIS_COMPLETE)

void SQUARE_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);
#endif


static int
fclaw2d_map_query_fivepatch(fclaw_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW_MAP_QUERY_IS_USED:
        return 1;
    case FCLAW_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW_MAP_QUERY_IS_AFFINE:
        return 1;
    case FCLAW_MAP_QUERY_IS_NONLINEAR:
        return 0;
    case FCLAW_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW_MAP_QUERY_IS_PLANAR:
        return 1;
    case FCLAW_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW_MAP_QUERY_IS_FLAT:
        return 1;
    case FCLAW_MAP_QUERY_IS_DISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_SPHERE:
        return 0;
    case FCLAW_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
    case FCLAW_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    case FCLAW_MAP_QUERY_IS_FIVEPATCH:
        return 1;
    case FCLAW_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_fivepatch (fclaw2d_map_query_defs.h) : " \
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}

#if 0
static void
fclaw2d_map_c2m_basis_fivepatch(fclaw_map_context_t * cont,
                               double xc, double yc, 
                               double *t, double *tinv, 
                               double *tderivs, int flag)
{
    /* Mapping must work for single Cartesian grid.   
       [xc,yc] \in [0,1]x[0,1] */
    SQUARE_BASIS_COMPLETE(&xc,&yc, t, tinv, tderivs, &flag);
}
#endif


static void
fclaw2d_map_c2m_fivepatch(fclaw_map_context_t* cont, int blockno,
                          double xc, double yc,
                          double *xp, double *yp, double *zp)
{
    double alpha = cont->user_double[0];

    /* Maps to [-1,1]x[-1,1] */
    MAPC2M_FIVEPATCH(&blockno,&xc,&yc,xp,yp,zp,&alpha);

    /* Scale to [-0.5,-0.5]x[0.5,0.5] */
    scale_map(cont, xp,yp,zp);

    /* Shift to [0,1]x[0,1] */
    shift_map(cont, xp,yp,zp);
}


fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha)
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_fivepatch;
    cont->mapc2m = fclaw2d_map_c2m_fivepatch;
    // cont->basis = fclaw2d_map_c2m_basis_fivepatch;

    set_scale(cont,scale);
    set_shift(cont,shift);

    cont->user_double[0] = alpha;

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
