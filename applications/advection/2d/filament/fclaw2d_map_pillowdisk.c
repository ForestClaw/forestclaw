/* Spherical disk in xy plane.  Matches p4est_connectivity_new_disk (). */

#include <fclaw2d_map.h>



static int
fclaw2d_map_query_pillowdisk (fclaw2d_map_context_t * cont, int query_identifier)
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
        return 1;
    case FCLAW2D_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_generic (fclaw2d_map_query_defs.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_pillowdisk(fclaw2d_map_context_t * cont, int blockno,
                      double xc, double yc,
                      double *xp, double *yp, double *zp)
{
    MAPC2M_PILLOWDISK(&xc,&yc,xp,yp,zp);

    /* These can probably be replaced by C functions at some point. */
    SCALE_MAP(xp,yp,zp);
    ROTATE_MAP(xp,yp,zp);
    SHIFT_MAP(xp,yp,zp);

}


fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk(const double rotate[],
                                                  const double scale)
{
    fclaw2d_map_context_t *cont;
    double shift[3];
    int m;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_pillowdisk;
    cont->mapc2m = fclaw2d_map_c2m_pillowdisk;

    /* This stores rotate/scale parameters in common blocks for later
       retrieval by scale_map/rotate_map (called above).  These parameters
       can of course be stored as variables in a context field */
    SET_ROTATION(rotate);
    SET_SCALE(&scale);
    for (m = 0; m < 2; m++)
    {
        shift[m] = 1.0;
    }
    SET_SHIFT(shift);

    return cont;
}
