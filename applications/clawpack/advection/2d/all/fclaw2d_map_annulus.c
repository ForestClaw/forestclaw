/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static int
fclaw2d_map_query_annulus (fclaw_map_context_t * cont, int query_identifier)
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
    case FCLAW_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_annulus (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_annulus.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_annulus (fclaw_map_context_t * cont, int blockno,
                       double xc, double yc,
                       double *xp, double *yp, double *zp)
{
    double xc1,yc1,zc1;

    /* Scale's brick mapping to [0,1]x[0,1] */
    /* fclaw_map_context_t *brick_map = (fclaw_map_context_t*) cont->user_data; */
    FCLAW_MAP_2D_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    /* blockno is ignored in the current annulus mapping;  it just assumes
       a single "logical" block in [0,1]x[0,1] */
    double beta = cont->user_double[0];
    double theta[2];
    theta[0] = cont->user_double[1];
    theta[1] = cont->user_double[2];
    MAPC2M_ANNULUS(&blockno,&xc1,&yc1,xp,yp,zp,&beta,theta);

    scale_map(cont, xp,yp,zp);
    rotate_map(cont, xp, yp, zp);
}

fclaw_map_context_t *
    fclaw2d_map_new_annulus (fclaw_map_context_t* brick,
                             const double scale[],
                             const double rotate[],
                             const double beta,
                             const double theta[])
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_annulus;
    cont->mapc2m = fclaw2d_map_c2m_annulus;

    cont->user_double[0] = beta;
    cont->user_double[1] = theta[0];
    cont->user_double[2] = theta[1];

    set_scale(cont,scale);
    set_rotate(cont,rotate);

    cont->brick = brick;

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
