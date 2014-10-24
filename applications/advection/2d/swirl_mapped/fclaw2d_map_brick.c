/* Cartesian grid, tranformed to Ax + b */

#include <fclaw2d_map.h>

static int
fclaw2d_map_query_brick (fclaw2d_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW2D_MAP_QUERY_IS_USED:
        return 1;
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
    default:
        printf("\n");
        printf("fclaw2d_map_query_cart (fclaw2d_map_cart.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_cart.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_brick(fclaw2d_map_context_t * cont, int blockno,
                     double xc, double yc,
                     double *xp, double *yp, double *zp)
{
    /* brick in [0,mi]x[0,mj] */
    int mi, mj;
    mi = cont->user_int[0];
    mj = cont->user_int[1];
    MAPC2M_BRICK(&blockno,&xc,&yc,xp,yp,zp,&mi, &mj);

    /* Scale and shift to [0,1]x[0,1] */
    scale_map(cont, xp,yp,zp);
    shift_map(cont, xp,yp,zp);
}


fclaw2d_map_context_t* fclaw2d_map_new_brick(const double scale[],
                                             const double shift[],
                                             const double rotate[],
                                             int mi,
                                             int mj)
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_brick;
    cont->mapc2m = fclaw2d_map_c2m_brick;
    cont->user_int[0] = mi;
    cont->user_int[1] = mj;

    set_scale(cont,scale);
    set_shift(cont,shift);
    set_rotate(cont,rotate);

    return cont;
}
