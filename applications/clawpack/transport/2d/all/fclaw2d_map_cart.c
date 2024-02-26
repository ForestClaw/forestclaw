/* Cartesian grid, tranformed to Ax + b */

#include <fclaw_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define SQUARE_BASIS_COMPLETE FCLAW_F77_FUNC(square_basis_complete, \
                                             SQUARE_BASIS_COMPLETE)

void SQUARE_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);

static int
fclaw2d_map_query_cart (fclaw_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW_MAP_QUERY_IS_USED:
        return 1;
    case FCLAW_MAP_QUERY_IS_SCALEDSHIFT:
        return 1;
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
        return 0;
    case FCLAW_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_cart (fclaw2d_map_cart.c) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_cart.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_basis_cart(fclaw_map_context_t * cont,
                           double xc, double yc, 
                           double *t, double *tinv, 
                           double *tderivs, int flag)
{
    /* Mapping must work for single Cartesian grid.   
       [xc,yc] \in [0,1]x[0,1] */
    SQUARE_BASIS_COMPLETE(&xc, &yc, t, tinv, tderivs, &flag);
}


static void
fclaw2d_map_c2m_cart(fclaw_map_context_t * cont, int blockno,
                     double xc, double yc,
                     double *xp, double *yp, double *zp)
{
    /* Brick mapping to computational coordinates [0,1]x[0,1] */
    double xc1, yc1, zc1;
    FCLAW_MAP_2D_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    /* Unit square in [-1,1] x [-1,1] */
    FCLAW_MAP_2D_C2M_CART(&blockno,&xc1,&yc1,xp,yp,zp);

    scale_map(cont, xp,yp,zp);
    shift_map(cont, xp,yp,zp);
}


fclaw_map_context_t* fclaw2d_map_new_cart(fclaw_map_context_t *brick,
                                            const double scale[],
                                            const double shift[])
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_cart;
    cont->mapc2m = fclaw2d_map_c2m_cart;
    cont->basis = fclaw2d_map_c2m_basis_cart;
    cont->brick = brick;

    set_scale(cont,scale);
    set_shift(cont,shift);

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
