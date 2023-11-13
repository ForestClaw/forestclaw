/* Cartesian grid, tranformed to Ax + b */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif


static int
fclaw2d_map_query_nomap_brick (fclaw_map_context_t * cont, int query_identifier)
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
        return 1;
    default:
        printf("\n");
        printf("fclaw2d_map_query_nomap_brick (fclaw2d_map_nomap_brick.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_nomap_brick.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


void
    fclaw2d_map_c2m_nomap_brick(fclaw_map_context_t * cont, int blockno,
                                double xc, double yc,
                                double *xp, double *yp, double *zp)
{
    /* Brick mapping to computational coordinates [0,1]x[0,1] */
    FCLAW_MAP_2D_BRICK2C(&cont,&blockno,&xc,&yc,xp,yp,zp);    
}


/* This shouldn't be called */
static void
fclaw3dx_map_c2m_nomap_brick(fclaw_map_context_t * cont, int blockno,
                             double xc, double yc, double zc,
                             double *xp, double *yp, double *zp)
{
    fclaw_global_essentialf("fclaw3dx_map_c2m_nomap_brick : Why is this being " \
                            "called? Maybe manifold=T.\n");
    exit(0);

    /* Call 2d mapping to get surface.  2d mapping is not scaled in the 
       extruded case. */
    cont->mapc2m(cont,blockno,xc,yc,xp,yp,zp);

    /* zc is already scaled into [az,bz] */
    *zp = zc;
}


fclaw_map_context_t* fclaw2d_map_new_nomap_brick(fclaw_map_context_t* brick)
{
    fclaw_map_context_t *cont;
    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_nomap_brick;
    cont->mapc2m = fclaw2d_map_c2m_nomap_brick;

    cont->mapc2m_3dx = fclaw3dx_map_c2m_nomap_brick;

    cont->brick = brick;

    /* this shouldn't really be referenced ... */
    cont->is_extruded = 0;

    return cont;
}

#ifdef __cplusplus
}
#endif
