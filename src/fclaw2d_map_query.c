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

#ifndef REFINE_DIM
#define REFINE_DIM 2
#endif

#ifndef PATCH_DIM
#define PATCH_DIM 2
#endif

#if PATCH_DIM == 2

#include <fclaw2d_map_query.h>

#elif PATCH_DIM == 3

#include <fclaw3d_map_query.h>
#include <_fclaw2d_to_fclaw3d.h>

#endif

int FCLAW2D_MAP_IS_USED(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_USED;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_DISK(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_DISK;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_CART(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_CART;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_AFFINE(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_AFFINE;
    iresult = cont->query(cont,id);
    return iresult != 0;
}


int FCLAW2D_MAP_IS_PILLOWDISK(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_PILLOWDISK;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_SQUAREDDISK(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_SQUAREDDISK;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int  FCLAW2D_MAP_IS_PILLOWSPHERE(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_PILLOWSPHERE;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_CUBEDSPHERE(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_CUBEDSPHERE;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_FLAT(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_FLAT;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_SPHERE(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_SPHERE;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_HEMISPHERE(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_HEMISPHERE;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_TORUS(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_TORUS;
    iresult = cont->query(cont,id);
    return iresult != 0;
}

int FCLAW2D_MAP_IS_BRICK(fclaw2d_map_context_t** pcont)
{
    int iresult;
    fclaw2d_map_context_t *cont = *pcont;
    int id = FCLAW2D_MAP_QUERY_IS_BRICK;
    iresult = cont->query(cont,id);
    return iresult != 0;
}
