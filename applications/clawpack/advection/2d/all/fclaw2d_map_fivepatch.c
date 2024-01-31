/* Five bilinear patches */

#include <fclaw_map.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* fix syntax highlighting */    
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


static void
fclaw2d_map_c2m_fivepatch(fclaw_map_context_t* cont, int blockno,
                          double xc, double yc,
                          double *xp, double *yp, double *zp)
{
    /* five patch square in [-1,1]x[-1,1] */
    double alpha = cont->user_double[0];
    FCLAW_MAP_2D_C2M_FIVEPATCH(&blockno,&xc,&yc,xp,yp,zp,&alpha);

    if (cont->is_extruded == 0)
    {
        /* Shift [-1,1]x[-1,1] to [0,2]x[0,2] */
        scale_map(cont, xp, yp, zp);
        shift_map(cont, xp,yp,zp);        
    }

}

#if 0
static void
fclaw3dx_map_c2m_fivepatch(fclaw2d_map_context_t* cont, int blockno,
                           double xc, double yc,double zc,
                           double *xp, double *yp, double *zp)
{
    /* five patch square in [-1,1]x[-1,1] */
    double alpha = cont->user_double[0];
    MAPC2M_FIVEPATCH(&blockno,&xc,&yc,xp,yp,zp,&alpha);
    *zp = 8*zc;

    /* Shift [-1,1]x[-1,1] to [0,2]x[0,2] */
    scale_map(cont, xp, yp, zp);
    shift_map(cont, xp,yp,zp);

}
#endif


fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha)
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_fivepatch;
    cont->mapc2m = fclaw2d_map_c2m_fivepatch;

    //cont->mapc2m_3dx = fclaw3dx_map_c2m_fivepatch;

    set_scale(cont,scale);
    set_shift(cont,shift);

    cont->user_double[0] = alpha;

    cont->is_extruded = 0;

    return cont;
}

#ifdef __cplusplus
}
#endif
