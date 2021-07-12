/* Pillow grid surface.  Matches p4est_connectivity_new_pillow (). */

#include <fclaw2d_map.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


#define MAPC2M_ANNULUS2 FCLAW_F77_FUNC(mapc2m_annulus2,MAPC2M_ANNULUS2)

void MAPC2M_ANNULUS2(const double* xc1, const double *yc1, 
                   double* xp, double *yp, double *zp);


#define ANNULUS_BASIS_COMPLETE FCLAW_F77_FUNC(annulus_basis_complete, \
                            ANNULUS_BASIS_COMPLETE)

void ANNULUS_BASIS_COMPLETE(const double* x, const double *y,
                           double t[], double tinv[], double uderivs[], 
                           const int* flag);



static int
fclaw2d_map_query_annulus (fclaw2d_map_context_t * cont, int query_identifier)
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
    case FCLAW2D_MAP_QUERY_IS_BRICK:
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
fclaw2d_map_c2m_basis_annulus(fclaw2d_map_context_t * cont,
                            double xc, double yc, 
                            double *t, double *tinv, 
                            double *tderivs, int flag)
{
    /* These coordinates are in [0,1]x[0,1] and are mapped to 
       [theta,phi] using mappings in sphere_basis.f */
    ANNULUS_BASIS_COMPLETE(&xc, &yc, t, tinv, tderivs, &flag);
}


static void
fclaw2d_map_c2m_annulus (fclaw2d_map_context_t * cont, int blockno,
                       double xc, double yc,
                       double *xp, double *yp, double *zp)
{
#if 0
    double beta, theta[2];
    beta     = cont->user_double[0];
    theta[0] = cont->user_double[1];
    theta[1] = cont->user_double[2];
#endif    

    /* Scale's brick mapping to [0,1]x[0,1] */
    /* fclaw2d_map_context_t *brick_map = (fclaw2d_map_context_t*) cont->user_data; */
    double xc1,yc1,zc1;
    FCLAW2D_MAP_BRICK2C(&cont,&blockno,&xc,&yc,&xc1,&yc1,&zc1);

    /* blockno is ignored in the current annulus mapping;  it just assumes
       a single "logical" block in [0,1]x[0,1] */
    MAPC2M_ANNULUS2(&xc1,&yc1,xp,yp,zp);
}

fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double shift[],
                             const double rotate[],
                             const double beta, const double theta[])
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_annulus;
    cont->mapc2m = fclaw2d_map_c2m_annulus;
    cont->basis = fclaw2d_map_c2m_basis_annulus;

    cont->user_double[0] = beta;
    cont->user_double[1] = theta[0];
    cont->user_double[2] = theta[1];

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
