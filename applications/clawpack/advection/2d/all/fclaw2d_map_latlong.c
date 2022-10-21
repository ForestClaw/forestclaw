/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw2d_map.h>

#include "advection_user.h"

#if 0
/* Fix syntax highlighting */
#endif    


static int
fclaw2d_map_query_latlong (fclaw2d_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW2D_MAP_QUERY_IS_USED:
        return 1;
        /* no break necessary after return statement */
    case FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_AFFINE:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_NONLINEAR:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_CART:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PLANAR:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_FLAT:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_DISK:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_SPHERE:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_PILLOWDISK:
        return 0;
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
        printf("fclaw2d_map_query_latlong (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_latlong.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_latlong (fclaw2d_map_context_t * cont, int blockno,
                       double xc, double yc,
                       double *xp, double *yp, double *zp)
{
    double lat[2], longitude[2];
    double xc1,yc1,zc1,xc2,yc2;

    /* fclaw2d_map_context_t *brick_map = (fclaw2d_map_context_t*) cont->user_data; */

    /* Scale's brick mapping to [0,1]x[0,1] to create a single "logical" block */
    FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    /* Scale into lat/long box */
    lat[0] = cont->user_double[0];
    lat[1] = cont->user_double[1];
    longitude[0] = cont->user_double[2];
    longitude[1] = cont->user_double[3];
    xc2 = longitude[0] + (longitude[1]-longitude[0])*xc1;
    yc2 = lat[0] + (lat[1]-lat[0])*yc1;

    /* blockno is ignored in the current latlong mapping;  it just assumes
       a single "logical" block in [long0,long1]x[lat0,lat1] */
    MAPC2M_LATLONG(&blockno,&xc2,&yc2,xp,yp,zp);

    if (cont->is_extruded == 0)
    {        
        scale_map(cont,xp,yp,zp);
        rotate_map(cont,xp,yp,zp);
    }
}



fclaw2d_map_context_t *
    fclaw2d_map_new_latlong (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double rotate[],
                             const double lat[],
                             const double longitude[],
                             const int a, const int b)
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_latlong;
    cont->mapc2m = fclaw2d_map_c2m_latlong;    

    cont->user_double[0] = lat[0];
    cont->user_double[1] = lat[1];
    cont->user_double[2] = longitude[0];
    if (a == 1)
    {
        // Make sure periodic case as full 360 degrees.
        cont->user_double[3] = longitude[0] + 360.0;
    }
    else
    {
        cont->user_double[3] = longitude[1];
    }

    set_scale(cont,scale);
    set_rotate(cont,rotate);

    cont->brick = brick;

    cont->is_extruded = 0;

    return cont;
}
