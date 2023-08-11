/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_INTEGRATE_H
#define FCLAW2D_INTEGRATE_H

#include "fclaw2d_convenience.h"  /* Needed for def. of fclaw2d_integrate_t */

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* fix syntax highlighting */
#endif    

struct fclaw2d_global;
struct fclaw_domain;



typedef struct fclaw2d_ray_vtable fclaw2d_ray_vtable_t;

struct fclaw2d_global;

/* This is a polynmorphic type */
typedef struct fclaw2d_ray
{
    int num;                 /* User defined ID */
    void* ray_data;          /* User defined */
    double integral;         /**< scalar integral value lives here */
    int untrustworthy;       /**< set to nonzero if integration
                                  result may be inaccurate */
} fclaw2d_ray_t;


typedef void (*fclaw2d_ray_allocate_and_define_t)(struct fclaw2d_global *glob, 
                                                  fclaw2d_ray_t **rays, 
                                                  int *num);

typedef void (*fclaw2d_ray_deallocate_t)(struct fclaw2d_global *glob, 
                                         fclaw2d_ray_t **rays, 
                                         int *num);



#if 0
typedef void (*fclaw2d_ray_create_files_t)(struct fclaw2d_global *glob, 
                                           struct fclaw2d_ray *rays, 
                                           int num_rays);

typedef void (*fclaw2d_ray_normalize_t)(struct fclaw2d_global *glob, 
                                       struct fclaw2d_block *block,
                                       int blockno, 
                                       struct fclaw2d_ray *g,
                                       double *xc, double *yc);


typedef void (*fclaw2d_ray_update_t)(struct fclaw2d_global* glob, 
                                     struct fclaw2d_block* block,
                                     struct fclaw2d_patch* patch, 
                                     int blockno, int patchno,
                                     double tcurr, struct fclaw2d_ray *g);

typedef void (*fclaw2d_ray_print_t)(struct fclaw2d_global *glob, 
                                    struct fclaw2d_ray *ray);

typedef void (*fclaw2d_ray_destroy_buffer_data_t)(struct fclaw2d_global *glob, 
                                                  void* gdata);

#endif


struct fclaw2d_ray_vtable
{
    fclaw2d_ray_allocate_and_define_t   allocate_and_define;
    fclaw2d_ray_deallocate_t            deallocate;

    fclaw2d_integrate_ray_t             integrate;   /* Function that does the integration */
#if 0
    fclaw2d_ray_create_files_t  create_ray_files;
    fclaw2d_ray_update_t        update_ray;
    fclaw2d_ray_print_t         print_ray_buffer;
    fclaw2d_ray_normalize_t     normalize_coordinates;
#endif    

    int is_set;
};


void fclaw2d_ray_allocate_and_define(struct fclaw2d_global* glob, 
                                     fclaw2d_ray_t **rays, 
                                     int *num_rays);

void fclaw2d_ray_deallocate(struct fclaw2d_global* glob, 
                            fclaw2d_ray_t **rays, 
                            int *num_rays);

void fclaw2d_ray_set_ray(fclaw2d_ray_t *r, 
                         int id, 
                         void* ray_data);

void* fclaw2d_ray_get_ray(fclaw2d_ray_t *r, 
                          int *id);

fclaw2d_ray_vtable_t* fclaw2d_ray_vt(struct fclaw2d_global *glob);

void fclaw2d_ray_vtable_initialize(struct fclaw2d_global *glob);

fclaw2d_ray_t* fclaw2d_ray_allocate_rays(int num_rays);

int fclaw2d_ray_deallocate_rays(fclaw2d_ray_t **rays);


#endif

#ifdef __cplusplus
}
#endif
