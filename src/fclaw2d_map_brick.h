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

#ifndef FCLAW2D_MAP_BRICK_H
#define FCLAW2D_MAP_BRICK_H

#include <fclaw_base.h>
#include <fclaw2d_domain.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw2d_map_context;



#define FCLAW2D_MAP_BRICK_GET_DIM FCLAW_F77_FUNC (fclaw2d_map_brick_get_dim, \
                                                  FCLAW2D_MAP_BRICK_GET_DIM)

void FCLAW2D_MAP_BRICK_GET_DIM(struct fclaw2d_map_context **cont,
                               int *mi, int* mj);



/* This is used for communicating with Matlab about the mapping */
#define WRITE_BRICK_DATA FCLAW_F77_FUNC (write_brick_data,WRITE_BRICK_DATA)
void WRITE_BRICK_DATA(int* n,
                      int* mi,
                      int* mj,
                      double xv[],
                      double yv[]);

typedef struct fclaw2d_block_ll
{
    int nb;
    int mi, mj;
    double *xv;
    double *yv;
}
fclaw2d_block_ll_t;

struct fclaw2d_map_context*
fclaw2d_map_new_brick (fclaw_domain_t *domain,
                       int mi, int mj, int periodic_i, int periodic_j);

void fclaw2d_map_destroy_brick (struct fclaw2d_map_context *cont);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
