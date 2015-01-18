/* Cartesian grid, tranformed to Ax + b */

#include <fclaw2d_map.h>
#include <fclaw2d_map_brick.h>

#define MAPC2M_MOUNTAIN FCLAW_F77_FUNC (mapc2m_mountain,MAPC2M_MOUNTAIN)
void MAPC2M_MOUNTAIN (int* blockno, double *xc, double *yc,
                      double *xp, double *yp, double *zp, double *scale);

#define MOUNTAIN_HEIGHT FCLAW_F77_FUNC (mountain_height,MOUNTAIN_HEIGHT)
double MOUNTAIN_HEIGHT (double *xp);


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static int
fclaw2d_map_query_mountain (fclaw2d_map_context_t * cont, int query_identifier)
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
        printf("\n");
        printf("fclaw2d_map_query_cart (fclaw2d_map_cart.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_cart.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_mountain(fclaw2d_map_context_t * cont, int blockno,
                         double xc, double yc,
                         double *xp, double *yp, double *zp)
{
    double xcorner[4],ycorner[4],xll,yll,h[4];
    double zbot, ztop,f;
    int i;

    fclaw2d_map_context_t* brick = cont->brick;
    fclaw2d_block_ll_t *bv = (fclaw2d_block_ll_t *) brick->user_data;

    /* Get xp */
    *xp = (double) cont->scale[0]*(bv->xv[blockno] + xc)/bv->mi;

    /* Get yp */
    xll = bv->xv[blockno];
    yll = bv->yv[blockno];
    xcorner[0] = xll;
    xcorner[1] = xll+1;
    xcorner[2] = xll+1;
    xcorner[3] = xll;

    ycorner[0] = yll;
    ycorner[1] = yll;
    ycorner[2] = yll+1;
    ycorner[3] = yll+1;
    for(i = 0; i < 4; i++)
    {
        xcorner[i] = cont->scale[0]*xcorner[i]/bv->mi;
        ycorner[i] = cont->scale[1]*ycorner[i]/bv->mj;
        h[i] = MOUNTAIN_HEIGHT(&xcorner[i]); /* Height of the four corners */
        f = ycorner[i]/cont->scale[1];
        ycorner[i] = h[i] + f*(cont->scale[1] - h[i]);
    }

    ztop = ycorner[3] + xc*(ycorner[2] - ycorner[3]);
    if (yll == 0)
    {
        zbot = MOUNTAIN_HEIGHT(xp);
    }
    else
    {
        zbot = ycorner[0] + xc*(ycorner[1] - ycorner[0]);
    }
    *yp = zbot + yc*(ztop - zbot);
    *zp = 0;
}


fclaw2d_map_context_t* fclaw2d_map_new_mountain(fclaw2d_map_context_t* brick,
                                                const double scale[],
                                                const double shift[],
                                                const double rotate[])
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_mountain;
    cont->mapc2m = fclaw2d_map_c2m_mountain;

    set_scale(cont,scale);
    set_shift(cont,shift);
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
