/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "fclaw2d_p4est.h"
#include <p4est.h>

static const double fclaw2d_p4est_root_len = (double) P4EST_ROOT_LEN;
/* static const double fclaw2d_p4est_smallest_h = 1. / (double) P4EST_ROOT_LEN; */

static unsigned fclaw2d_p4est_magic;
#define FCLAW2D_P4EST_MAGIC \
  (fclaw2d_p4est_magic ?: (fclaw2d_p4est_magic = \
   sc_hash_function_string ("fclaw2d_p4est_magic", NULL)))

static int
fclaw2d_p4est_map_query (fclaw_map_context_t * cont, int query_identifier)
{
    P4EST_ASSERT (cont->magic == FCLAW2D_P4EST_MAGIC);

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
        return 0;
    }
    return 0;
}

static void
fclaw2d_p4est_map_c2m (fclaw_map_context_t * cont, int blockno,
                       double cx, double cy,
                       double *mx, double *my, double *mz)
{
    P4EST_ASSERT (cont->magic == FCLAW2D_P4EST_MAGIC);

    p4est_connectivity_t *conn = (p4est_connectivity_t *) cont->user_data;
    double xyzp[3];

    p4est_qcoord_to_vertex (conn, (p4est_topidx_t) blockno,
                            (p4est_qcoord_t) (cx * fclaw2d_p4est_root_len),
                            (p4est_qcoord_t) (cy * fclaw2d_p4est_root_len),
                            xyzp);
    *mx = xyzp[0];
    *my = xyzp[1];
    *mz = xyzp[2];
}

fclaw_map_context_t *
fclaw2d_p4est_map_new (const p4est_connectivity_t * conn)
{
    fclaw_map_context_t *cont;

    P4EST_ASSERT (conn->vertices != NULL);

    cont = P4EST_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->magic = FCLAW2D_P4EST_MAGIC;
    cont->query = fclaw2d_p4est_map_query;
    cont->mapc2m = fclaw2d_p4est_map_c2m;
    cont->user_data = (void *) conn;

    return cont;
}

void
fclaw2d_p4est_map_destroy (fclaw_map_context_t * cont)
{
    P4EST_ASSERT (cont->magic == FCLAW2D_P4EST_MAGIC);
    P4EST_FREE (cont);
}
