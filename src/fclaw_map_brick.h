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

#ifndef FCLAW_MAP_BRICK_H
#define FCLAW_MAP_BRICK_H

#include <fclaw_base.h>
#include <fclaw_domain.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_map_context;



#define FCLAW_MAP_BRICK_GET_2D_DIM FCLAW_F77_FUNC (fclaw_map_brick_get_2d_dim, \
                                                  FCLAW_MAP_BRICK_GET_2D_DIM)

void FCLAW_MAP_BRICK_GET_2D_DIM(struct fclaw_map_context **cont,
                               int *mi, int* mj);



/* This is used for communicating with Matlab about the mapping */
#define WRITE_BRICK_DATA FCLAW_F77_FUNC (write_brick_data,WRITE_BRICK_DATA)
void WRITE_BRICK_DATA(int* n,
                      int* mi,
                      int* mj,
                      double xv[],
                      double yv[]);

typedef struct fclaw_block_ll
{
    int nb;
    int mi, mj, mk; /**< Number of blocks in each direction */
    double *xv;
    double *yv;
    double *zv;
}
fclaw_block_ll_t;

/**
 * @brief Creates a new 2D brick map context.
 *
 * @param domain The domain to which the map context belongs.
 * @param mi The number of cells in the x-direction of the brick.
 * @param mj The number of cells in the y-direction of the brick.
 * @param periodic_i Flag indicating whether the x-direction is periodic.
 * @param periodic_j Flag indicating whether the y-direction is periodic.
 * @return A pointer to the newly created map context.
 */
struct fclaw_map_context*
fclaw_map_new_2d_brick (fclaw_domain_t *domain,
                       int mi, int mj, int periodic_i, int periodic_j);

/**
 * @brief Creates a new 2D brick map context.
 *
 * @param domain The domain to which the map context belongs.
 * @param mi The number of cells in the x-direction of the brick.
 * @param mj The number of cells in the y-direction of the brick.
 * @param periodic_i Flag indicating whether the x-direction is periodic.
 * @param periodic_j Flag indicating whether the y-direction is periodic.
 * @return A pointer to the newly created map context.
 */
struct fclaw_map_context*
fclaw_map_new_3d_brick (fclaw_domain_t *domain,
                       int mi, int mj, int mk,
                       int periodic_i, int periodic_j, int periodic_k);


/**
 * @brief Destroys a brick mapping context.
 *
 * @param cont The fclaw_map_context to destroy
 */
void fclaw_map_destroy_brick (struct fclaw_map_context *cont);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
