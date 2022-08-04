/* Cubed sphere surface.  Matches p4est_connectivity_new_cubed (). */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#endif


static int
fclaw2d_map_query_cubedsphere (fclaw2d_map_context_t * cont, int query_identifier)
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
        return 1;
    case FCLAW2D_MAP_QUERY_IS_BRICK:
        return 0;
    default:
        printf("\n");
        printf("fclaw2d_map_query_cubedsphere (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_query_defs.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_cubedsphere (fclaw2d_map_context_t * cont, int blockno,
                         double xc, double yc,
                         double *xp, double *yp, double *zp)
{
    MAPC2M_CUBEDSPHERE(&blockno, &xc,&yc,xp,yp,zp);

    scale_map(cont,xp,yp,zp); 

    /* These can probably be replaced by C functions at some point. */
#if 0    
    /* These could be causing problems ... */
    rotate_map(cont,xp,yp,zp);
#endif    
}

static void
fclaw3dx_map_c2m_cubedsphere (fclaw2d_map_context_t * cont, int blockno,
                              double xc, double yc, double zc, 
                              double *xp, double *yp, double *zp)

{
    /* Maps point on blockno to sphere */
    MAPC2M_CUBEDSPHERE(&blockno, &xc,&yc,xp,yp,zp);

    /* Map point on sphere to point in shell */
    double phi = asin(*zp);     // returns value in [-pi/2, pi/2]
    double theta = atan2(*yp,*xp);      // returns value in [-pi, pi]
    if (theta < 0)
        theta += 2*M_PI;

    double maxelev = cont->user_double[0];  /* should match value set in the options */

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


fclaw2d_map_context_t *
    fclaw2d_map_new_cubedsphere(const double scale[],
                                const double rotate[])
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_cubedsphere;
    cont->mapc2m = fclaw2d_map_c2m_cubedsphere;
    cont->mapc2m_3dx = fclaw3dx_map_c2m_cubedsphere;

    set_rotate(cont,rotate);
    set_scale(cont,scale);

    /* 
       Should match value set in options. 

       Note : Does it make sense to pass this through as an argument 
       to the mapping, since this mapping is used both for 2d and 3d?
       Changing the argument list would break lots of examples ... 
    */
    cont->user_double[0] = 0.5;   



    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
