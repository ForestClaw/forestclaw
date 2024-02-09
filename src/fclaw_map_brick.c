/* Cartesian grid, tranformed to Ax + b */

#include <fclaw_map_brick.h>

#include <fclaw_map.h>
#include <fclaw_convenience.h>

static int
fclaw2d_map_query_brick (fclaw_map_context_t * cont, int query_identifier)
{
    switch (query_identifier)
    {
    case FCLAW_MAP_QUERY_IS_USED:
        return 1;
    case FCLAW_MAP_QUERY_IS_SCALEDSHIFT:
        return 1;
    case FCLAW_MAP_QUERY_IS_AFFINE:
        return 1;
    case FCLAW_MAP_QUERY_IS_NONLINEAR:
        return 0;
    case FCLAW_MAP_QUERY_IS_GRAPH:
        return 0;
    case FCLAW_MAP_QUERY_IS_PLANAR:
        return 1;
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
    case FCLAW_MAP_QUERY_IS_FIVEPATCH:
        return 0;
    case FCLAW_MAP_QUERY_IS_CART:
        return 1;
    case FCLAW_MAP_QUERY_IS_BRICK:
        return 1;
    default:
        printf("\n");
        printf("fclaw2d_map_query_brick (fclaw2d_map_brick.h) : "\
               "Query id not identified;  Maybe the query is not up to "\
               "date?\nSee fclaw2d_map_brick.c.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}

void
FCLAW_MAP_BRICK_GET_2D_DIM(fclaw_map_context_t **cont,
                          int *mi, int* mj)
{
    fclaw_block_ll_t *bv = (fclaw_block_ll_t *) (*cont)->brick->user_data;
    *mi = bv->mi;
    *mj = bv->mj;
}


static void
fclaw_map_c2m_2d_brick(fclaw_map_context_t * cont, int blockno,
                     double xc, double yc,
                     double *xp, double *yp, double *zp)
{
    /* Map the brick coordinates to global [0,1] coordinates. */
    fclaw_block_ll_t *bv = (fclaw_block_ll_t *) cont->user_data;
    *xp = (double) (bv->xv[blockno] + xc)/bv->mi;
    *yp = (double) (bv->yv[blockno] + yc)/bv->mj;
    *zp = 0;
}
static void
fclaw_map_c2m_3d_brick(fclaw_map_context_t * cont, int blockno,
                     double xc, double yc, double zc,
                     double *xp, double *yp, double *zp)
{
    /* Map the brick coordinates to global [0,1] coordinates. */
    fclaw_block_ll_t *bv = (fclaw_block_ll_t *) cont->user_data;
    *xp = (double) (bv->xv[blockno] + xc)/bv->mi;
    *yp = (double) (bv->yv[blockno] + yc)/bv->mj;
    *zp = (double) (bv->zv[blockno] + zc)/bv->mk;
}

void fclaw_map_destroy_brick(fclaw_map_context_t *cont)
{
    fclaw_block_ll_t *bv = (fclaw_block_ll_t *) cont->user_data;
    FCLAW_FREE(bv->xv);
    FCLAW_FREE(bv->yv);
    FCLAW_FREE(bv->zv);

    FCLAW_FREE (cont->user_data);
    FCLAW_FREE (cont);
}

static
fclaw_block_ll_t* new_block_data(fclaw_domain_t *domain,
                                 int mi, int mj, int mk,
                                 int periodic_i, int periodic_j, int periodic_k)
{
    int nb = domain->num_blocks;

    fclaw_block_ll_t *bv = FCLAW_ALLOC_ZERO(fclaw_block_ll_t,1);

    /* We don't store this in user_double[], since that only has limited
       storage (16 doubles) */
    bv->xv = FCLAW_ALLOC_ZERO(double,nb);
    bv->yv = FCLAW_ALLOC_ZERO(double,nb);
    bv->zv = FCLAW_ALLOC_ZERO(double,nb);

    bv->nb = nb;  /* These integer values could also be stored in user_int[] data */
    bv->mi = mi;
    bv->mj = mj;
    bv->mk = mk;

    for (int i = 0; i < nb; i++)
    {
        fclaw_block_t *block = &domain->blocks[i];

        /* (x,y) reference coordinates of lower-left block corner */
        bv->xv[i] = block->vertices[0];
        bv->yv[i] = block->vertices[1];
        bv->zv[i] = block->vertices[2];
    }
    return bv;
}

fclaw_map_context_t*
fclaw_map_new_2d_brick (fclaw_domain_t *domain,
                       int mi, int mj, int periodic_i, int periodic_j)
{
    fclaw_map_context_t* cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_brick;
    cont->mapc2m = fclaw_map_c2m_2d_brick;
    cont->destroy = fclaw_map_destroy_brick;

    fclaw_block_ll_t *bv = new_block_data(domain, mi, mj, 1, periodic_i, periodic_j, 0);

    cont->user_data = (void*) bv;

    /* Write data for Matlab plotting */
    FILE *fid;
    fid = fopen("brick.dat","w");
    fprintf(fid,"%8d %8d\n",mi, mj);
    for(int i = 0; i < bv->nb; i++)
    {
        fprintf(fid,"%12g %12g\n",bv->xv[i], bv->yv[i]);
    }
    fclose(fid);
    return cont;
}

fclaw_map_context_t*
fclaw_map_new_3d_brick (fclaw_domain_t *domain,
                       int mi, int mj, int mk, int periodic_i, int periodic_j, int periodic_k)
{
    fclaw_map_context_t* cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_brick;
    cont->mapc2m_3d = fclaw_map_c2m_3d_brick;
    cont->destroy = fclaw_map_destroy_brick;

    fclaw_block_ll_t *bv = new_block_data(domain, mi, mj, mk, periodic_i, periodic_j, periodic_k);

    cont->user_data = (void*) bv;

    /* Write data for Matlab plotting */
    FILE *fid;
    fid = fopen("brick.dat","w");
    fprintf(fid,"%8d %8d\n",mi, mj);
    for(int i = 0; i < bv->nb; i++)
    {
        fprintf(fid,"%12g %12g %12g\n",bv->xv[i], bv->yv[i], bv->zv[i]);
    }
    fclose(fid);
    return cont;
}
