/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw2d_map.h>

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

    int mi, mj;
    double lat[2], longitude[2];

    mi = cont->user_int[0];
    mj = cont->user_int[1];
    MAPC2M_BRICK(&blockno,&xc,&yc,xp,yp,zp,&mi, &mj);

    /* map back to [0,1]x[0,1] */
    xc = (double) *xp/mi;
    yc = (double) *yp/mj;

    /* Scale into lat/long box */
    lat[0] = cont->user_double[0];
    lat[1] = cont->user_double[1];
    longitude[0] = cont->user_double[2];
    longitude[1] = cont->user_double[3];
    xc = longitude[0] + (longitude[1]-longitude[0])*xc;
    yc = lat[0] + (lat[1]-lat[0])*yc;

    /* blockno is ignored in the current latlong mapping;  it just assumes
       a single "logical" block in [0,1]x[0,1] */
    MAPC2M_LATLONG(&blockno,&xc,&yc,xp,yp,zp);

    scale_map(cont,xp,yp,zp);
    rotate_map(cont,xp,yp,zp);
}

fclaw2d_map_context_t *
    fclaw2d_map_new_latlong (const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double lat[],
                             const double longitude[],
                             const int mi,
                             const int ni,
                             const int a,
                             const int b)
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
        cont->user_double[3] = longitude[0] + 360.0;
    }
    cont->user_int[0] = mi;
    cont->user_int[1] = ni;

    set_scale(cont,scale);
    set_shift(cont,shift);
    set_rotate(cont,rotate);

    return cont;
}
