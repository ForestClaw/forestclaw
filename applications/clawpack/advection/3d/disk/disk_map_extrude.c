#include "disk_user.h"

/* For mappings headers */
#include "../all/advection_user.h"

static void
disk_map_3dx(fclaw2d_map_context_t * cont, 
             int blockno,
             double xc, double yc, double zc,
             double *xp, double *yp, double *zp)
{
    /* Get surface map for one of two disk mappings. Because this 
       surface will be extruded, the 2d mapping won't be scaled or 
       shifted or rotated. */
    cont->mapc2m(cont,blockno,xc,yc,xp,yp,zp);

    double maxelev = cont->user_double_3dx[0];  

    /* Extrude surface map in z direction */
    *zp = maxelev*zc;

    scale_map(cont,xp,yp,zp);
    shift_map(cont,xp,yp,zp);
}

void disk_map_extrude(fclaw2d_map_context_t *cont, 
                      const double maxelev)
{
    /* Modify the 2d mapping to included extruded details */
    cont->mapc2m_3dx = disk_map_3dx;
    cont->user_double_3dx[0] = maxelev;

    cont->is_extruded = 1;
}

