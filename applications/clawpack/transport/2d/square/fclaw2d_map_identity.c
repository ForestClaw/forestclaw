/* Cartesian grid, tranformed to Ax + b */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define CART_BASIS_COMPLETE FCLAW_F77_FUNC(cart_basis_complete, CART_BASIS_COMPLETE)

void CART_BASIS_COMPLETE(const int* blockno, const double* x, const double *y,
                          double t[], double tinv[], double uderivs[], 
                          const int* flag);

static int
fclaw2d_map_query_identity (fclaw2d_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW2D_MAP_QUERY_IS_USED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_AFFINE:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_NONLINEAR:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PLANAR:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FLAT:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_DISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FIVEPATCH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_identity (fclaw2d_map_identity.c) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_identity.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}

static void
fclaw2d_map_c2m_basis_identity(fclaw2d_map_context_t * cont, int blockno,
                               double xc, double yc, 
                               double *t, double *tinv, 
                               double *tderivs, int flag)
{
    /* Brick mapping to computational coordinates [0,1]x[0,1] */
    double xc1, yc1, zc1;
    FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    CART_BASIS_COMPLETE(&blockno, &xc1,&yc1, t, tinv, tderivs, &flag);
}


static void
    fclaw2d_map_c2m_identity(fclaw2d_map_context_t * cont, int blockno,
                             double xc, double yc,
                             double *xp, double *yp, double *zp)
{

    /* Brick mapping to computational coordinates [0,1]x[0,1] */
    double xc1, yc1, zc1;
    FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    *xp = xc1;
    *yp = yc1;
    *zp = 0;
}


fclaw2d_map_context_t* fclaw2d_map_new_identity(fclaw2d_map_context_t *brick)
{
    fclaw2d_map_context_t *cont;
    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_identity;
    cont->mapc2m = fclaw2d_map_c2m_identity;
    cont->basis = fclaw2d_map_c2m_basis_identity;

    cont->brick = brick;

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
