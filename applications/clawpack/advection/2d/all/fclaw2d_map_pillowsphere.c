/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw_map.h>

#ifdef __cplusplus
extern "C"
{
#endif


static int
fclaw2d_map_query_pillowsphere (fclaw_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW_MAP_QUERY_IS_USED:
        return 1;
        /* no break necessary after return statement */
    case FCLAW_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW_MAP_QUERY_IS_AFFINE:
        return 0;
    case FCLAW_MAP_QUERY_IS_NONLINEAR:
        return 1;
    case FCLAW_MAP_QUERY_IS_CART:
        return 0;
    case FCLAW_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW_MAP_QUERY_IS_PLANAR:
        return 0;
    case FCLAW_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW_MAP_QUERY_IS_FLAT:
        return 0;
    case FCLAW_MAP_QUERY_IS_DISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_SPHERE:
        return 0;
    case FCLAW_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_PILLOWSPHERE:
        return 1;
    case FCLAW_MAP_QUERY_IS_CUBEDSPHERE:
        return 0;
    case FCLAW_MAP_QUERY_IS_FIVEPATCH:
        return 0;
    case FCLAW_MAP_QUERY_IS_HEMISPHERE:
        return 1;
    case FCLAW_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_pillowsphere (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_pillowsphere (fclaw_map_context_t * cont, int blockno,
                              double xc, double yc,
                              double *xp, double *yp, double *zp)
{
    FCLAW_MAP_2D_C2M_PILLOWSPHERE(&blockno,&xc,&yc,xp,yp,zp);

    if (cont->is_extruded == 0)
    {
        scale_map(cont,xp,yp,zp);
        rotate_map(cont,xp,yp,zp);        
    }
}

fclaw_map_context_t *
    fclaw2d_map_new_pillowsphere(const double scale[],
                                 const double rotate[])
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_pillowsphere;
    cont->mapc2m = fclaw2d_map_c2m_pillowsphere;
    
    set_scale(cont,scale); 
    set_rotate(cont, rotate);

    cont->is_extruded = 0;

    return cont;
}
#ifdef __cplusplus
}
#endif
