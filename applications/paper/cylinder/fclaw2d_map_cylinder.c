/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif



#define MAPC2M_CYLINDER FCLAW_F77_FUNC(mapc2m_cylinder,MAPC2M_CYLINDER)

void MAPC2M_CYLINDER(const double* xc1, const double *yc1, 
                   double* xp, double *yp, double *zp);


#define CYLINDER_BASIS_COMPLETE FCLAW_F77_FUNC(cylinder_basis_complete, \
                            CYLINDER_BASIS_COMPLETE)

void CYLINDER_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);


#define MAPC2M_LATLONG2 FCLAW_F77_FUNC(mapc2m_latlong2, MAPC2M_LATLONG2)

void MAPC2M_LATLONG2(const double* xc1, const double *yc1, 
                   double* xp, double *yp, double *zp);


#define LATLONG_BASIS_COMPLETE FCLAW_F77_FUNC(latlong_basis_complete, \
                            LATLONG_BASIS_COMPLETE)

void LATLONG_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);

static int
fclaw2d_map_query_cylinder (fclaw2d_map_context_t * cont, int query_identifier)
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
        printf("fclaw2d_map_query_cylinder (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_cylinder.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}


static void
fclaw2d_map_c2m_basis_cylinder(fclaw2d_map_context_t * cont,
                               double xc, double yc, 
                               double *t, double *tinv, 
                               double *tderivs, int flag)
{
    /* These coordinates are in [0,1]x[0,1] and are mapped to 
       [theta,phi] using mappings in sphere_basis.f */
    int mapping = cont->user_int[0];
    if (mapping == 0)
        CYLINDER_BASIS_COMPLETE(&xc, &yc, t, tinv, tderivs, &flag);
    else
        LATLONG_BASIS_COMPLETE(&xc, &yc, t, tinv, tderivs, &flag);
}


static void
fclaw2d_map_c2m_cylinder(fclaw2d_map_context_t * cont, int blockno,
                         double xc, double yc,
                         double *xp, double *yp, double *zp)
{
    /* Data is not already in brick domain */
    double xc1,yc1,zc1; /* We don't need zc1 - we are we computing it? */
    FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    int mapping = cont->user_int[0];
    if (mapping == 0)
        MAPC2M_CYLINDER(&xc1,&yc1,xp,yp,zp);
    else
        MAPC2M_LATLONG2(&xc1,&yc1,xp,yp,zp);
}

fclaw2d_map_context_t *
    fclaw2d_map_new_cylinder (fclaw2d_map_context_t* brick,
                              const double scale[], int mapping)
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_cylinder;
    cont->mapc2m = fclaw2d_map_c2m_cylinder;
    cont->basis = fclaw2d_map_c2m_basis_cylinder;

    cont->brick = brick;
    cont->user_int[0] = mapping;

    set_scale(cont,scale);

    return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
