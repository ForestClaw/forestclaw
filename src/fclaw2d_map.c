/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef P4_TO_P8

#include <fclaw_global.h>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>  /* Needed for pillowsphere query */

#define MAP_KEY "map_2d"

#else

#include <fclaw3d_map.h>

#define MAP_KEY "map_3d"

#endif

#ifndef P4_TO_P8

void
fclaw_map_store (fclaw_global_t* glob,
                          fclaw_map_context_t * map)
{
    fclaw_global_attribute_store(glob, MAP_KEY, map, NULL);
}

fclaw_map_context_t*
fclaw_map_get(fclaw_global_t* glob)
{
    return (fclaw_map_context_t*) fclaw_global_get_attribute(glob,MAP_KEY);
}

/* This function can be called from Fortran inside of ClawPatch. */
void
FCLAW_MAP_QUERY (fclaw_map_context_t ** pcont,
                 const int *query_identifier, int *iresult)
{
    fclaw_map_context_t *cont = *pcont;
    *iresult = cont->query (cont, *query_identifier);
}

/* This function can be called from Fortran inside of ClawPatch. */

void
FCLAW_MAP_2D_C2M (fclaw_map_context_t ** pcont, int *blockno,
                 const double *xc, const double *yc,
                 double *xp, double *yp, double *zp)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->mapc2m != NULL);
    cont->mapc2m (cont, *blockno, *xc, *yc, xp, yp, zp);
}


/* This function can be called from Fortran inside of ClawPatch. */

void
FCLAW2D_MAP_C2M_BASIS (fclaw_map_context_t ** pcont,
                       const double *xc, const double *yc,
                       double *t, double *tinv, double *tderivs,
                       int *flag)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->basis != NULL);

    cont->basis (cont, *xc, *yc, t, tinv, tderivs, *flag);
}


/* This pragma requires a fclaw3dx_map.{c,h} file, so I"ll hold off for 
   now on adding these.   Ideally, the patch metric won't know whether
   the mapping is for an extruded mesh or for a full octree mesh. 
*/

//#if PATCH_DIM == 3 && REFINE_DIM == 2
void
FCLAW3D_MAP_C2M (fclaw_map_context_t ** pcont, int *blockno,
                 const double *xc, const double *yc, const double *zc,
                 double *xp, double *yp, double *zp)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->mapc2m_3dx != NULL);
    cont->mapc2m_3dx (cont, *blockno, *xc, *yc, *zc, xp, yp, zp);
}



void
FCLAW3D_MAP_C2M_BASIS (fclaw_map_context_t ** pcont,
                       const double *xc, const double *yc, const double *zc, 
                       double *t, double *tinv, double *tderivs,
                       int *flag)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->basis_3dx != NULL);

    /* DAC : Still need to decide what should go here */
    cont->basis_3dx (cont, *xc, *yc, *zc, t, tinv, tderivs, *flag);
}
//#endif

void FCLAW2D_MAP_BRICK2C (fclaw_map_context_t ** pcont, int *blockno,
                          const double *xc, const double *yc,
                          double *xp, double *yp, double *zp)
{
    fclaw_map_context_t *cont = *pcont;

    if (cont->brick != NULL)
    {
        fclaw_map_context_t *brick = cont->brick;  
        brick->mapc2m (brick, *blockno, *xc, *yc, xp, yp, zp);
    }
    else
    {
        /* We only have one tree */
        FCLAW_ASSERT(blockno == 0);
        *xp = *xc;
        *yp = *yc;
        *zp = 0;
    }
}


/* This function is expected to be called from C or C++. */
void
fclaw2d_map_destroy (fclaw_map_context_t * cont)
{
#ifndef P4_TO_P8
    if (cont->brick != NULL)
    {
        fclaw2d_map_destroy(cont->brick);  /* recursive call */
    }
#endif
    if (cont->destroy == NULL)
    {
        FCLAW_FREE (cont);
    }
    else
    {
        cont->destroy (cont);
    }
}

#endif /* !P4_TO_P8 */

#ifndef P4_TO_P8

#if 0
/* Cubed sphere surface.  Matches p4est_connectivity_new_cubed (). */

static int
fclaw2d_map_query_csphere (fclaw_map_context_t * cont, int query_identifier)
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
    default:
        printf("\n");
        printf("fclaw2d_map_query_csphere (fclaw2d_map.c) : Query id not "\
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_query.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}

static inline void
fclaw2d_map_c2m_csphere_help (double R, double xi, double eta,
                             double *x, double *y, double *z)
{
    const double tan_xi = tan (.5 * M_PI * (xi - .5));
    const double tan_eta = tan (.5 * M_PI * (eta - .5));

    /*
      This won't be true for ghost cell values.
    FCLAW_ASSERT (0. <= xi && xi <= 1.);
    FCLAW_ASSERT (0. <= eta && eta <= 1.);
    */

    *z = R / sqrt (SC_SQR (tan_xi) + SC_SQR (tan_eta) + 1.);
    *x = *z * tan_xi;
    *y = *z * tan_eta;
}

static void
fclaw2d_map_c2m_csphere (fclaw_map_context_t * cont, int blockno,
                         double xc, double yc,
                         double *xp, double *yp, double *zp)
{
    const double R = cont->user_double[0];

    switch (blockno)
    {
    case 0:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, yp, xp, zp);
        *zp = -*zp;
        break;
    case 1:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, zp, xp, yp);
        break;
    case 2:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, zp, yp, xp);
        *xp = -*xp;
        break;
    case 3:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, xp, yp, zp);
        break;
    case 4:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, xp, zp, yp);
        *yp = -*yp;
        break;
    case 5:
        fclaw2d_map_c2m_csphere_help (R, xc, yc, yp, zp, xp);
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }
}

fclaw_map_context_t *
fclaw2d_map_new_csphere (double R)
{
    fclaw_map_context_t *cont;

    FCLAW_ASSERT (0. <= R);

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_csphere;
    cont->mapc2m = fclaw2d_map_c2m_csphere;
    cont->user_double[0] = R;

    return cont;
}

/* Spherical disk in xy plane.  Matches p4est_connectivity_new_disk (). */

static int
fclaw2d_map_query_disk (fclaw_map_context_t * cont, int query_identifier)
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
    default:
        printf("\n");
        printf("fclaw2d_map_query_disk (fclaw2d_map.c) : Query id not " \
               "identified;  Maybe the query is not up to date?\nSee "  \
               "fclaw2d_map_query.h.\n");
        printf("Requested query id : %d\n",query_identifier);
        SC_ABORT_NOT_REACHED ();
    }
    return 0;
}

static inline void
fclaw2d_map_c2m_disk_help (double R2sqrbyR1, double R1byR2,
                           double xi, double eta, double *x, double *y)
{
    const double R = R2sqrbyR1 * pow (R1byR2, 1. + eta);
    const double tan_xi = tan (.5 * M_PI * (xi - .5));
    const double xi_prime = (1. - eta) * 2. * (xi - .5) + eta * tan_xi;

    *y = R / sqrt (1. + eta * SC_SQR (tan_xi) + (1. - eta));
    *x = *y * xi_prime;
}

static void
fclaw2d_map_c2m_disk (fclaw_map_context_t * cont, int blockno,
                      double xc, double yc,
                      double *xp, double *yp, double *zp)
{
    FCLAW_ASSERT (0. <= xc && xc <= 1.);
    FCLAW_ASSERT (0. <= yc && yc <= 1.);

    *zp = 0.;
    if (blockno == 2) {
        /* center square */
        const double half_length = cont->user_double[2];

        *xp = (2. * xc - 1.) * half_length;
        *yp = (2. * yc - 1.) * half_length;
    }
    else {
        /* outer patches */
        const double R2sqrbyR1 = cont->user_double[0];
        const double R1byR2 = cont->user_double[1];

        switch (blockno)
        {
        case 0:
            fclaw2d_map_c2m_disk_help (R2sqrbyR1, R1byR2,
                                       xc, 1. - yc, xp, yp);
            *yp = -*yp;
            break;
        case 1:
            fclaw2d_map_c2m_disk_help (R2sqrbyR1, R1byR2,
                                       yc, 1. - xc, yp, xp);
            *xp = -*xp;
            break;
        case 3:
            fclaw2d_map_c2m_disk_help (R2sqrbyR1, R1byR2, yc, xc, yp, xp);
            break;
        case 4:
            fclaw2d_map_c2m_disk_help (R2sqrbyR1, R1byR2, xc, yc, xp, yp);
            break;
        default:
            SC_ABORT_NOT_REACHED ();
        }
    }
}

fclaw_map_context_t *
fclaw2d_map_new_disk (double R1, double R2)
{
    fclaw_map_context_t *cont;

    FCLAW_ASSERT (0. < R2 && R2 <= R1);

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_disk;
    cont->mapc2m = fclaw2d_map_c2m_disk;
    cont->user_double[0] = R2 * R2 / R1;
    cont->user_double[1] = R1 / R2;
    cont->user_double[2] = R2 / M_SQRT2;        /* half length of square */

    return cont;
}
#endif

#if 0
/* Use an existing Fortran mapc2m routine.
 * The answers to the queries are expected in user_int[0] through [4].
 * The pointer to the Fortran mapping function is stored in user_data.
 */

static int
fclaw2d_map_query_fortran (fclaw_map_context_t * cont, int query_identifier)
{
    return 0 <= query_identifier && query_identifier <
        FCLAW2D_MAP_QUERY_LAST ? cont->user_int[query_identifier] : 0;
}

static void
fclaw2d_map_c2m_fortran (fclaw_map_context_t * cont, int blockno,
                         double xc, double yc,
                         double *xp, double *yp, double *zp)
{
    /* call Fortran functions */
    set_block_ (&blockno);
    (*(fclaw2d_map_c2m_fortran_t) cont->user_data) (&xc, &yc, xp, yp, zp);
}

fclaw_map_context_t *
fclaw2d_map_new_fortran (fclaw2d_map_c2m_fortran_t mapc2m,
                         const int query_results[FCLAW2D_MAP_QUERY_LAST])
{
    fclaw_map_context_t *cont;

    cont = FCLAW_ALLOC_ZERO (fclaw_map_context_t, 1);
    cont->query = fclaw2d_map_query_fortran;
    cont->mapc2m = fclaw2d_map_c2m_fortran;
    memcpy (cont->user_int, query_results,
            FCLAW2D_MAP_QUERY_LAST * sizeof (int));
    cont->user_data = (void *) mapc2m;

    return cont;
}
#endif

#if 0
void fclaw2d_register_map_data(sc_options_t* options, fclaw2d_map_data_t* map_data)
{
    sc_options_add_int (options, 0, "mapping:mi", &map_data->mi, 1,
                        "[mapping] Number of blocks in x direction [1]");

    sc_options_add_int (options, 0, "mapping:mj", &map_data->mj, 1,
                        "[mapping] Number of blocks in y direction  [1]");

    sc_options_add_int (options, 0, "mapping:periodic_x", &map_data->periodic_x, 0,
                        "[mapping] Periodic in x direction [0]");

    sc_options_add_int (options, 0, "mapping:periodic_y", &map_data->periodic_y, 0,
                        "[mapping] Periodic in y direction  [0]");

    /* --------------------------------------------------------------------
       Scale
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (options,0, "mapping:scale",
                                    &map_data->scale_string, "1 1 1",
                                    &map_data->scale, 3,
                                    "[mapping] Scale factor [1 1 1]");

    /* --------------------------------------------------------------------
       Shift
       --------------------------------------------------------------------*/
    fclaw_options_add_double_array (options,0, "mapping:shift",
                                    &map_data->shift_string, "0 0 0",
                                    &map_data->shift, 3,
                                    "[mapping] Shift array [0 0 0]");

    /* --------------------------------------------------------------------
       Rotate
       --------------------------------------------------------------------*/
    sc_options_add_double (options, 0, "mapping:phi", &map_data->phi, 0,
                           "[mapping] Rotation angle about x axis (degrees) [0]");

    sc_options_add_double (options, 0, "mapping:theta", &map_data->theta, 0,
                           "[mapping] Rotation angle about z axis (degrees) [0]");

}
#endif

#if 0
void fclaw2d_map_destroy_arrays(fclaw2d_map_data_t* map_data)
{
    FCLAW_FREE (map_data->scale);
    FCLAW_FREE (map_data->shift);
}
#endif

#if 0
void
fclaw2d_options_postprocess_map_data(fclaw2d_map_data_t * map_data)
{
    fclaw_options_convert_double_array (map_data->scale_string, &map_data->scale,3);
    fclaw_options_convert_double_array (map_data->shift_string, &map_data->shift,3);
}
#endif

int fclaw2d_map_pillowsphere(fclaw_global_t* glob)
{
    fclaw_map_context_t *cont = fclaw_map_get(glob);
    return FCLAW2D_MAP_IS_PILLOWSPHERE(&cont) != 0;    
}

void set_scale(fclaw_map_context_t* cont, const double scale[])
{
    memcpy(cont->scale,scale,3*sizeof(double));
}

void set_shift(fclaw_map_context_t* cont, const double shift[])
{
    memcpy(cont->shift,shift,3*sizeof(double));
}

void set_default_transform(double scale[],double shift[],double rotate[])
{
  shift[0] = 0;
  shift[1] = 0;
  shift[2] = 0;
  scale[0] = 1;
  scale[1] = 1;
  scale[2] = 1;
  rotate[0] = 0;
  rotate[1] = 0;
}

void set_rotate(fclaw_map_context_t* cont, const double rotate[])
{
    double rotate_mat[9];
    SET_ROTATION_MATRIX(rotate,rotate_mat);
    memcpy(cont->rotate,rotate_mat,9*sizeof(double));
}

void scale_map(fclaw_map_context_t* cont, double *xp, double *yp, double *zp)
{
    *xp *= cont->scale[0];
    *yp *= cont->scale[1];
    *zp *= cont->scale[2];
}

void shift_map(fclaw_map_context_t* cont, double *xp, double *yp, double *zp)
{
    *xp += cont->shift[0];
    *yp += cont->shift[1];
    *zp += cont->shift[2];
}

void rotate_map(fclaw_map_context_t* cont, double *xp, double *yp, double *zp)
{
    double v[3], vrot[3];
    int i,j;

    v[0] = *xp;
    v[1] = *yp;
    v[2] = *zp;

    for(i = 0; i < 3; i++)
    {
        vrot[i] = 0;
        for(j = 0; j < 3; j++)
        {
            vrot[i] += cont->rotate[3*j + i]*v[j];
        }
    }
    *xp = vrot[0];
    *yp = vrot[1];
    *zp = vrot[2];
}

#endif /* !P4_TO_P8 */
