/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton,
Hannes Brandt
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

#ifndef OVERLAP_H
#define OVERLAP_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <fclaw2d_options.h>

#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>    
#include <fclaw2d_global.h>


typedef struct overlap_prodata
{
  double myvalue[7];
  int    isset;
}
overlap_prodata_t;

typedef struct overlap_point
{
  size_t  lnum;
  double  xy[3];
  overlap_prodata_t   prodata;
}
overlap_point_t;

typedef struct overlap_consumer
{
  fclaw2d_global_t   *glob;
  fclaw2d_domain_t   *domain;
  sc_array_t         *query_points;
  size_t              cell_idx;
  int                 num_cells_in_patch;
}
overlap_consumer_t;


typedef struct overlap_geometry
{
    fclaw_options_t *fclaw_opt;
    fclaw2d_block_t *blocks;
}
overlap_geometry_t;


void apply_consumer_mapping (overlap_point_t * op);

void create_query_points (overlap_consumer_t * c);

int apply_inverse_producer_mapping (overlap_point_t * op, double xy[3],
                                    int blockno, overlap_geometry_t * geo);

int overlap_interpolate (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
                         int blockno, int patchno, void *point, void *user);

void output_query_points (overlap_consumer_t * c);


void add_cell_centers (fclaw2d_domain_t * domain, fclaw2d_patch_t * patch,
                       int blockno, int patchno, void *user);


#ifdef __cplusplus
}
#endif

#endif

