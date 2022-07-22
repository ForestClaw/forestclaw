/* Spherical disk in xy plane.  Matches p4est_connectivity_new_disk (). */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


static int
fclaw2d_map_query_pillowdisk5 (fclaw2d_map_context_t * cont, int query_identifier)
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
    case FCLAW2D_MAP_QUERY_IS_BRICK:
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
fclaw2d_map_c2m_pillowdisk5(fclaw2d_map_context_t * cont, int blockno,
                                double xc, double yc,
                                double *xp, double *yp, double *zp)
{
<<<<<<< Updated upstream
    double alpha = cont->user_double[0];
=======

    /* Unit disk centered at (0,0) */
    double alpha = cont->user_double[0];    
>>>>>>> Stashed changes
    MAPC2M_PILLOWDISK5(&blockno,&xc,&yc,xp,yp,zp,&alpha);

    scale_map(cont, xp,yp,zp);
    rotate_map(cont,xp,yp,zp);
}


fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
<<<<<<< Updated upstream
=======
                                                   const double shift[],
>>>>>>> Stashed changes
                                                   const double rotate[],
                                                   const double alpha)
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_pillowdisk5;
    cont->mapc2m = fclaw2d_map_c2m_pillowdisk5;

    cont->user_double[0] = alpha;

    set_scale(cont, scale);
<<<<<<< Updated upstream
=======
    set_shift(cont, shift);
>>>>>>> Stashed changes
    set_rotate(cont, rotate);

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
