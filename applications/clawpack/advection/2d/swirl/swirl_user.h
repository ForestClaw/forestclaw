/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef SWIRL_USER_H
#define SWIRL_USER_H

#include "../all/advection_user.h"

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* fix syntax highlighting */
#endif

typedef struct user_options
{
    double period;
    int claw_version;
    int is_registered;

} user_options_t;

void swirl_link_solvers(fclaw_global_t *glob);

/* ------------------------------------- Options ---------------------------------------*/
user_options_t* swirl_options_register (fclaw_app_t * app,
                                        const char *configfile);

void swirl_options_store (fclaw_global_t* glob, user_options_t* user);

const user_options_t* swirl_get_options(fclaw_global_t* glob);


/* ------------------------------------- Rays ---------------------------------------*/

#if 0
typedef enum 
{
  SWIRL_RAY_LINE,
  SWIRL_RAY_CIRCLE,
  SWIRL_RAY_TYPE_LAST
} swirl_ray_type_t;


typedef struct swirl_ray
{
  swirl_ray_type_t    rtype;
  double              xy[2];
  union {
    struct {
      double              vec[2];
    } line;
    struct {
      double              radius;
    } circle;
  } r;
} swirl_ray_t;
#endif

int swirl_intersect_ray(fclaw_domain_t *domain, 
                         fclaw_patch_t * patch,
                         int blockno, 
                         int patchno, 
                         void *ray, 
                         double *integral,
                         void* user);

#if 0
void swirl_define_ray(fclaw_global_t *glob, 
                      fclaw2d_ray_t** rays, 
                      int* num_rays);

void swirl_destroy_ray(fclaw_global_t *glob,void* ray);
#endif

void swirl_initialize_rays(fclaw_global_t* glob);


sc_array_t * swirl_integrals_new(void);
sc_array_t * swirl_rays_new (void);



#ifdef __cplusplus
}
#endif

#endif
