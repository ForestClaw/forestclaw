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

#include "swirl_user.h"

#include "../all/advection_user.h"

#include <fclaw2d_rays.h>

typedef enum 
{
    SWIRL_RAY_LINE,
    SWIRL_RAY_CIRCLE,
    SWIRL_RAY_TYPE_LAST
} swirl_ray_type_t;

typedef struct swirl_ray
{
    swirl_ray_type_t rtype;
    double xy[2];
    union 
    {
        struct 
        {
            double vec[2];
        } line;
        struct 
        {
            double radius;
        } circle;
    } r;
} swirl_ray_t;

int swirl_intersect_ray (fclaw2d_domain_t *domain, 
                         fclaw2d_patch_t * patch,
                         int blockno, 
                         int patchno, 
                         void *ray, 
                         double *integral)
{
  /* assert that ray is a valid swirl_ray_t */
  swirl_ray_t *swirl_ray = (swirl_ray_t *) ray;
  FCLAW_ASSERT(swirl_ray != NULL);
  FCLAW_ASSERT(swirl_ray->rtype == SWIRL_RAY_LINE); /* Circles not there yet. */

  if (patchno >= 0) 
  {
     /* We are at a leaf and the patch is a valid patch of the domain.
     * Based on the patch, the domain, the blockno and the information stored
     * in the swirl_ray_t we defined, we now have to set *integral to be the
     * contribution of this ray-patch combination to the ray integral.
     * We should return 1 (even though a leaf return value is ignored). */

    /* This is a dummy example.  We add the ray's x component for each patch.
       Truly, this example must be updated to compute the exact ray integral. */
    *integral = swirl_ray->r.line.vec[0];
    return 1;
  } 
  else 
  {
    /* We are not at a leaf and the patch is an artificial patch containing all
     * standard patch information except for the pointer to the next patch and
     * user-data of any kind. Only the FCLAW2D_PATCH_CHILDID and the
     * FCLAW2D_PATCH_ON_BLOCK_FACE_* flags are set.
     * Based on this, we now can run a test to check if the ray and the patch
     * intersect.
     * We return 0 if we are certain that the ray does not intersect any
     * descendant of this patch.
     * We return 1 if the test concludes that the ray may intersect the patch.
     * This test may be overinclusive / false positive to optimize for speed.
     *
     * The purpose of this test is to remove irrelevant ancestor
     * patch-ray-combinations early on to avoid unnecessary computations.
     * We do not need to assign to the integral value for ancestor patches. */

    /* This is a dummy example.  Truly, implement a fast and over-inclusive test
     * to see if this ray may possibly intersect the patch and return the answer. */
    return 1;
  }
}


static int nlines = 3;

sc_array_t * swirl_rays_new (void)
{
    swirl_ray_t *ray;
    sc_array_t  *a = sc_array_new (sizeof (swirl_ray_t));

    /* add a couple straight rays */
    for (int i = 0; i < nlines; ++i) 
    {
        ray = (swirl_ray_t *) sc_array_push (a);
        ray->rtype = SWIRL_RAY_LINE;
        ray->xy[0] = 0.;
        ray->xy[1] = 0.;
        ray->r.line.vec[0] = cos (i * M_PI / nlines);
        ray->r.line.vec[1] = sin (i * M_PI / nlines);
    }

    /* add no circles yet */
    return a;
}




/* Virtual function for setting rays */
static
void swirl_allocate_and_define_rays(fclaw2d_global_t *glob, 
                                    fclaw2d_ray_t** rays, 
                                    int* num_rays)
{
    *num_rays = nlines;

    /* We let the user allocate an array of rays, although what is inside the 
       generic ray type is left opaque. This is destroy in matching FREE,
       below. */
    *rays = (fclaw2d_ray_t*) FCLAW_ALLOC(fclaw2d_ray_t,*num_rays);
    for (int i = 0; i < nlines; ++i) 
    {
        swirl_ray_t *sr = (swirl_ray_t*) FCLAW_ALLOC(swirl_ray_t,1);
        sr->rtype = SWIRL_RAY_LINE;
        sr->xy[0] = 0.;
        sr->xy[1] = 0.;
        sr->r.line.vec[0] = cos (i * M_PI / nlines);
        sr->r.line.vec[1] = sin (i * M_PI / nlines);
        fclaw2d_ray_set_ray(glob,&(*rays)[i],i, sr);
    }
}


static
void swirl_deallocate_rays(fclaw2d_global_t *glob, 
                           fclaw2d_ray_t** rays, 
                           int* num_rays)
{
    for(int i = 0; i < *num_rays; i++)
    {
        /* Retrieve rays set above and deallocate them */
        int id;
        swirl_ray_t *rs = (swirl_ray_t*) fclaw2d_ray_get_ray(glob,&(*rays)[i],&id);
        FCLAW_ASSERT(rs != NULL);
        FCLAW_FREE(rs);
        rs = NULL;
    }
    /* Match FCLAW_ALLOC, above */
    FCLAW_FREE(*rays);  
    *num_rays = 0;  
    *rays = NULL;
}

void swirl_initialize_rays(fclaw2d_global_t* glob)
{
    /* Set up rays */
    fclaw2d_ray_vtable_t* rays_vt = fclaw2d_ray_vt(); 

    rays_vt->allocate_and_define = swirl_allocate_and_define_rays;
    rays_vt->deallocate = swirl_deallocate_rays;
}


sc_array_t * swirl_integrals_new(void)
{
    return sc_array_new_count (sizeof (double), nlines);
}