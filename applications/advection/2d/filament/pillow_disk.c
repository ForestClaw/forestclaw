/* Spherical disk in xy plane.  Matches p4est_connectivity_new_disk (). */

void mapc2m_disk_(double *xc, double *yc, double *xp, double *yp, double *zp);

static int
fclaw2d_map_query_pillowdisk (fclaw2d_map_context_t * cont, int query_identifier)
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
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_PLANAR:
        return 1;
    case FCLAW2D_MAP_QUERY_IS_ALIGNED:
        return 0;
    case FCLAW2D_MAP_QUERY_IS_PILLOWDISK:   // testing out this idea
        return 1;
    }
    return 0;
}


static void
fclaw2d_map_c2m_pillowdisk(fclaw2d_map_context_t * cont, int blockno,
                      double xc, double yc,
                      double *xp, double *yp, double *zp)
{
    FCLAW_ASSERT (0. <= xc && xc <= 1.);
    FCLAW_ASSERT (0. <= yc && yc <= 1.);
    set_block_(blockno);
    mapc2m_disk_(&xc,&yc,xp,yp,zp);
}

fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk()
{
    fclaw2d_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_pillowdisk;
    cont->mapc2m = fclaw2d_map_c2m_pillowdisk;

    return cont;
}
