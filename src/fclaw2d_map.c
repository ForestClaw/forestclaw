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

#include "fclaw2d_map.h"
#include <p4est_base.h>

/* This prototype is declared in src/clawpack_fort.H.  Clean up later. */
void set_block_ (const int *a_blockno);

/* This function can be called from Fortran inside of ClawPatch. */
void
fclaw2d_map_query_ (fclaw2d_map_context_t * cont,
                    const int *query_identifier, int *iresult)
{
    *iresult = cont->query (cont, *query_identifier);
}

/* This function can be called from Fortran inside of ClawPatch. */
void
fclaw2d_map_c2m_ (fclaw2d_map_context_t * cont, int *blockno,
                  const double *cx, const double *cy,
                  double *mx, double *my, double *mz)
{
    cont->mapc2m (cont, *blockno, *cx, *cy, mx, my, mz);
}

/* Torus.  Uses user_double[0,1] for R1 and R2, respectively. */

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
    case FCLAW2D_MAP_QUERY_IS_GRAPH:
        return 0;
    }
    return 0;
}

static void
fclaw2d_map_c2m_torus (fclaw2d_map_context_t * cont, int blockno,
                       double cx, double cy,
                       double *mx, double *my, double *mz)
{
    const double R2 = cont->user_double[1];
    const double L = cont->user_double[0] + R2 * cos (2. * M_PI * cy);

    *mx = L * cos (2. * M_PI * cx);
    *my = L * sin (2. * M_PI * cx);
    *mz = R2 * sin (2. * M_PI * cy);
}

fclaw2d_map_context_t *
fclaw2d_map_new_torus (double R1, double R2)
{
    fclaw2d_map_context_t *cont;

    P4EST_ASSERT (0. <= R1);
    P4EST_ASSERT (0. <= R2 && R2 <= R1);

    cont = P4EST_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_torus;
    cont->mapc2m = fclaw2d_map_c2m_torus;
    cont->user_double[0] = R1;
    cont->user_double[1] = R2;

    return cont;
}

void
fclaw2d_map_destroy_torus (fclaw2d_map_context_t * cont)
{
    P4EST_FREE (cont);
}

/* Use an existing Fortran mapc2m routine.
 * The answers to the queries are expected in user_int[0] through [4].
 * The pointer to the Fortran mapping function is stored in user_data.
 */

static int
fclaw2d_map_query_fortran (fclaw2d_map_context_t * cont, int query_identifier)
{
    return 0 <= query_identifier
        && query_identifier <
        FCLAW2D_MAP_QUERY_LAST ? cont->user_int[query_identifier] : 0;
}

static void
fclaw2d_map_c2m_fortran (fclaw2d_map_context_t * cont, int blockno,
                         double cx, double cy,
                         double *mx, double *my, double *mz)
{
    /* call Fortran functions */
    set_block_ (&blockno);
    (*(fclaw2d_map_c2m_fortran_t) cont->user_data) (&cx, &cy, mx, my, mz);
}

fclaw2d_map_context_t *
fclaw2d_map_new_fortran (fclaw2d_map_c2m_fortran_t mapc2m,
                         const int query_results[FCLAW2D_MAP_QUERY_LAST])
{
    fclaw2d_map_context_t *cont;

    cont = P4EST_ALLOC_ZERO (fclaw2d_map_context_t, 1);
    cont->query = fclaw2d_map_query_fortran;
    cont->mapc2m = fclaw2d_map_c2m_fortran;
    memcpy (cont->user_int, query_results,
            FCLAW2D_MAP_QUERY_LAST * sizeof (int));
    cont->user_data = (void *) mapc2m;

    return cont;
}

void
fclaw2d_map_destroy_fortran (fclaw2d_map_context_t * cont)
{
    P4EST_FREE (cont);
}
