/* Five bilinear patches */

#include <fclaw_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

static int
fclaw2d_map_query_pillowsphere5(fclaw_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW_MAP_QUERY_IS_USED:
        return 1;
    case FCLAW_MAP_QUERY_IS_SCALEDSHIFT:
        return 0;
    case FCLAW_MAP_QUERY_IS_AFFINE:
        return 0;
    case FCLAW_MAP_QUERY_IS_NONLINEAR:
        return 1;
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
        return 1;
    case FCLAW_MAP_QUERY_IS_PILLOWDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_SQUAREDDISK:
        return 0;
    case FCLAW_MAP_QUERY_IS_PILLOWSPHERE:
        return 0;
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
        printf("fclaw2d_map_query_fivepatch (fclaw2d_map_query_defs.h) : " \
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_pillowsphere5(fclaw_map_context_t * cont, int blockno,
                          double xc, double yc,
                          double *xp, double *yp, double *zp)
{
    double xp1, yp1, zp1;
    double alpha = cont->user_double[0];
    FCLAW_MAP_2D_C2M_FIVEPATCH(&blockno,&xc,&yc,&xp1,&yp1,&zp1,&alpha);
    xp1 = (xp1 + 1)/2.0;
    yp1 = (yp1 + 1)/2.0;
    FCLAW_MAP_2D_C2M_PILLOWSPHERE(&blockno,&xp1,&yp1,xp,yp,zp);

    /* These can probably be replaced by C functions at some point. */
    scale_map(cont,xp,yp,zp);
    rotate_map(cont,xp,yp,zp);

}

static void
fclaw3dx_map_c2m_pillowsphere5(fclaw_map_context_t * cont, int blockno,
                          double xc, double yc, double zc,
                          double *xp, double *yp, double *zp)
{

    /* map to fivepatch square in [-1,1]x[-1,1] */
    double xp1, yp1, zp1;
    double alpha = cont->user_double[0];
    FCLAW_MAP_2D_C2M_FIVEPATCH(&blockno,&xc,&yc,&xp1,&yp1,&zp1,&alpha);

    /* Map to unit square */
    xp1 = (xp1 + 1)/2.0;
    yp1 = (yp1 + 1)/2.0;

    /* Map to pillow sphere.  

    Note : Since we only have a single unit square, we will only have the 
    upper half of the sphere 
    */
    FCLAW_MAP_2D_C2M_PILLOWSPHERE(&blockno,&xp1,&yp1,xp,yp,zp);


    /* Map point on sphere to point in shell */
    double phi = asin(*zp);     // returns value in [-pi/2, pi/2]
    double theta = atan2(*yp,*xp);      // returns value in [-pi, pi]
    if (theta < 0)
        theta += 2*M_PI;

    double maxelev = cont->user_double[1];  /* should match value set in the options */

    /* Scaling depends on zc value, so we can't just scale using 'scale' below. */
    double R = maxelev*zc + 1;  
    *xp = R*cos(phi)*cos(theta);
    *yp = R*cos(phi)*sin(theta);
    *zp = R*sin(phi);

#if 0
    scale_map(cont,xp,yp,zp);
    rotate_map(cont,xp,yp,zp);
#endif    

}


fclaw_map_context_t* fclaw2d_map_new_pillowsphere5(const double scale[],
                                                     const double rotate[],
                                                     const double alpha)
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_pillowsphere5;
    cont->mapc2m = fclaw2d_map_c2m_pillowsphere5;
    cont->mapc2m_3dx = fclaw3dx_map_c2m_pillowsphere5;

    cont->user_double[0] = alpha;

    double maxelev = 0.5;
    cont->user_double[1] = maxelev;

    set_scale(cont,scale);
    set_rotate(cont,rotate);

    return cont;
}
#ifdef __cplusplus
#if 0
{
#endif
}
#endif
