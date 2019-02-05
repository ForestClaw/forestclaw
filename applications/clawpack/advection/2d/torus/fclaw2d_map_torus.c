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
fclaw2d_map_query_torus (fclaw2d_map_context_t * cont, int query_identifier)
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
        return 0;
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
        printf("fclaw2d_map_query_torus (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_torus.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_torus (fclaw2d_map_context_t * cont, int blockno,
                       double xc, double yc,
                       double *xp, double *yp, double *zp)
{
    double xc1,yc1,zc1;
    double alpha, beta;
    double L[4];
    double x,y, r, R, pi, pi2;
    int i;

    pi = M_PI;
    pi2 = 2*pi;

    /* Scale brick mapping to [0,1]x[0,1] */
    if (blockno >= 0)
    {
        /* Data is not already in brick domain */
        FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);
    }
    else
    {
        xc1 = xc;
        yc1 = yc;
        zc1 = 0;
    }

    /* blockno is ignored in the current torus mapping;  it just assumes
       a single "logical" block in [0,1]x[0,1] */
    alpha = cont->user_double[0];
    beta = cont->user_double[1];

    for(i=0; i < 4; i++)
    {
        L[i] = cont->user_double[2+i];
    }

    /* Map from (a1,a2) back to (xi,eta) */
    x = L[0]*xc1 + L[1]*yc1;
    y = L[2]*xc1 + L[3]*yc1;

    r = alpha*(1 + beta*sin(pi2*x));
    R = 1 + r*cos(pi2*y);
    
    *xp = R*cos(pi2*x);
    *yp = R*sin(pi2*x);
    *zp = r*sin(pi2*y);

/*    
    int example = cont->user_int[0];
    if (example == 0)
    {
        MAPC2M_TORUS(&blockno,&xc1,&yc1,xp,yp,zp,&alpha);
    }
    else if (example == 1)
    {
        MAPC2M_TWISTED_TORUS(&blockno,&xc1,&yc1,xp,yp,zp,&alpha);
    }
    rotate_map(cont,xp,yp,zp);
*/    
}

fclaw2d_map_context_t *
    fclaw2d_map_new_torus (fclaw2d_map_context_t* brick,
                           const double scale[],
                           const double shift[],
                           const double rotate[],
                           const double alpha,
                           const double beta,
                           const int mapping)
{
    int i;
    double l0[4] = {1.,  0.,  0.,  1.};
    double l1[4] = {1.,  0.,  1.,  1.};
    double l2[4] = {1.,  1.,  1.,  0.};

    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_torus;
    cont->mapc2m = fclaw2d_map_c2m_torus;

    cont->user_double[0] = alpha;
    cont->user_double[1] = beta;

    for(i = 0; i < 4; i++)
    {
        if (mapping == 0)
        {
            /* Regular torus mapping.  L given rowwise*/
            cont->user_double[2+i] = l0[i];
        }
        else if (mapping == 1)
        {
            cont->user_double[2+i] = l1[i];
        }
        else if (mapping == 2)
        {
            cont->user_double[2+i] = l2[i];
        }
    }

    cont->user_int[0] = mapping;  /* Might not be used */

    /* Are we really using these? */
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
