/* Spherical disk in xy plane.  Matches p4est_connectivity_new_disk (). */

#include <fclaw2d_map.h>



static int
fclaw2d_map_query_squareddisk (fclaw2d_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW2D_MAP_QUERY_IS_USED:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_AFFINE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_NONLINEAR:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PLANAR:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FLAT:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_DISK:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_SPHERE:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SQUAREDDISK:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_generic (fclaw2d_map_query_defs.h) : " \
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_squareddisk(fclaw2d_map_context_t * cont, int blockno,
                            double xc, double yc,
                            double *xp, double *yp, double *zp)
{
    double alpha = cont->user_double[0];
    MAPC2M_SQUAREDDISK(&blockno,&xc,&yc,xp,yp,zp,&alpha);

    /* These can probably be replaced by C functions at some point. */
    SCALE_MAP(xp,yp,zp);
    ROTATE_MAP(xp,yp,zp);
    SHIFT_MAP(xp,yp,zp);
}


fclaw2d_map_context_t* fclaw2d_map_new_squareddisk(const double rotate[],
                                                   const double scale,
                                                   const double alpha)
{
    fclaw2d_map_context_t *cont;
    double shift[3];

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_squareddisk;
    cont->mapc2m = fclaw2d_map_c2m_squareddisk;

    cont->user_double[0] = alpha;

    /* This stores rotate/scale parameters in common blocks for later
       retrieval by scale_map/rotate_map (called above).  These parameters
       can of course be stored as variables in a context field */
    SET_ROTATION(rotate);
    SET_SCALE(&scale);
    shift[0] = 1;
    shift[1] = 1;
    shift[2] = 0;
    SET_SHIFT(shift);

    return cont;
}
