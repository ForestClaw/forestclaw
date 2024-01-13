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
#include <fclaw_map.h>
#include <fclaw_map_query.h>  /* Needed for pillowsphere query */

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
FCLAW_MAP_2D_C2M_BASIS (fclaw_map_context_t ** pcont,
                        const double *xc, const double *yc,
                        double *t, double *tinv, double *tderivs,
                        int *flag)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->basis != NULL);

    cont->basis (cont, *xc, *yc, t, tinv, tderivs, *flag);
}


void
FCLAW_MAP_3D_C2M (fclaw_map_context_t ** pcont, int *blockno,
                  const double *xc, const double *yc, const double *zc,
                  double *xp, double *yp, double *zp)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->mapc2m_3dx != NULL);
    cont->mapc2m_3dx (cont, *blockno, *xc, *yc, *zc, xp, yp, zp);
}



void
FCLAW_MAP_3D_C2M_BASIS (fclaw_map_context_t ** pcont,
                        const double *xc, const double *yc, const double *zc, 
                        double *t, double *tinv, double *tderivs,
                        int *flag)
{
    fclaw_map_context_t *cont = *pcont;
    FCLAW_ASSERT(cont->basis_3dx != NULL);

    /* DAC : Still need to decide what should go here */
    cont->basis_3dx (cont, *xc, *yc, *zc, t, tinv, tderivs, *flag);
}

void FCLAW_MAP_2D_BRICK2C (fclaw_map_context_t ** pcont, int *blockno,
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
fclaw_map_destroy (fclaw_map_context_t * cont)
{
#ifndef P4_TO_P8
    if (cont->brick != NULL)
    {
        fclaw_map_destroy(cont->brick);  /* recursive call */
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

int fclaw_map_pillowsphere(fclaw_global_t* glob)
{
    fclaw_map_context_t *cont = fclaw_map_get(glob);
    return FCLAW_MAP_IS_PILLOWSPHERE(&cont) != 0;    
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
