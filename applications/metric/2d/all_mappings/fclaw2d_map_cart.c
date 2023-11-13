/* Cartesian grid, tranformed to Ax + b */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static int
fclaw2d_map_query_cart (fclaw_map_context_t * cont, int query_identifier)
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
        return 0;
#if 0
        printf("\n");
        printf("fclaw2d_map_query_cart (fclaw2d_map_cart.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_cart.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
#endif
    }
    return 0;
}


static void
fclaw2d_map_c2m_cart(fclaw_map_context_t * cont, int blockno,
                     double xc, double yc,
                     double *xp, double *yp, double *zp)
{
    /* Unit square in [-1,1]x[-1,1] */
    MAPC2M_CART(&blockno,&xc,&yc,xp,yp,zp);

    /*
       scale_map(cont,xp,yp,zp);
       shift_map(cont,xp,yp,zp);
       rotate_map(cont,xp,yp,zp);
    */
}

fclaw_map_context_t* fclaw2d_map_new_cart(const double scale[],
                                            const double shift[],
                                            const double rotate[])
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_cart;
    cont->mapc2m = fclaw2d_map_c2m_cart;

    /*
    set_scale(cont,scale);
    set_shift(cont,shift);
    set_rotate(cont,rotate);
    */


    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
