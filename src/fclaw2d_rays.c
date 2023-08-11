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

#include "fclaw2d_rays.h"
#include "fclaw_global.h"
#include "fclaw2d_options.h"
#include "fclaw2d_diagnostics.h"

#include "fclaw2d_convenience.h"

#include <fclaw_pointer_map.h>


//static fclaw2d_ray_vtable_t s_ray_vt;


typedef struct fclaw2d_ray_acc
{
    int num_rays;
    //sc_array_t *rays;
    //sc_array_t *integrals;
    fclaw2d_ray_t *rays;
} fclaw2d_ray_acc_t;



/* ------------------------------ Virtualized functions ------------------------------- */

/* Do I need this level of indirection? */
static
void ray_allocate_and_define(fclaw_global_t* glob, 
                             fclaw2d_ray_t **rays, 
                             int *num_rays)
{
    const fclaw2d_ray_vtable_t* ray_vt = fclaw2d_ray_vt(glob);

    FCLAW_ASSERT(ray_vt->allocate_and_define != NULL);
    ray_vt->allocate_and_define(glob, rays, num_rays);    
}

static
void ray_deallocate(fclaw_global_t* glob, 
                    fclaw2d_ray_t **rays, 
                    int *num_rays)
{
    const fclaw2d_ray_vtable_t* ray_vt = fclaw2d_ray_vt(glob);

    FCLAW_ASSERT(ray_vt->deallocate != NULL);
    ray_vt->deallocate(glob, rays, num_rays);    
}


static
void ray_initialize(fclaw_global_t* glob, void** acc)
{
    fclaw2d_ray_acc_t *ray_acc = FCLAW_ALLOC(fclaw2d_ray_acc_t,1);

    /* Check to see if user wants ray output */
    int i, num_rays;
    const fclaw_options_t * fclaw_opt = fclaw2d_get_options(glob);
    if (!fclaw_opt->output_rays)
    {
        num_rays = 0;
        ray_acc->rays = NULL;
    }
    else
    {
        /* This allocates memory for pointers to fclaw2d_ray_t types, 
           User defined rays are also allocated and stored. 
           We could allocate ray_acc->rays here, but we don't know ahead of time 
           how many rays to create.  Only the user knows this */
        ray_allocate_and_define(glob, &ray_acc->rays, &num_rays);

        /* We initialize untrustworthy to 0 */
        fclaw2d_ray_t *rays = ray_acc->rays;
        for (i =0; i < num_rays; i++)
        {
            fclaw2d_ray_t *ray = &rays[i];
            ray->untrustworthy = 0;
        }
    }
    *acc = ray_acc;
    ray_acc->num_rays = num_rays;

#if 0
    if (num_rays > 0)
    {
        acc->rays = sc_array_new (sizeof (fclaw2d_ray_t));        
        for(int i = 0; i < num_rays; i++)
        {
            fclaw2d_ray_t *fclaw_ray = (fclaw2d_ray_t *) sc_array_push (sc_rays);        
            fclaw_ray->num = ray_acc->rays[i].num;
            fclaw_ray->ray_data = (void*) ray_acc->rays[i].ray_data;
        }

        acc->integrals = sc_array_new_count (sizeof (double), num_rays);

        /* Things that could be done here : 
           1. Open ray files to write out time series data for each ray 
           2. Set integral to 0? 
        */
    }
#endif
}

static
void ray_reset(fclaw_global_t *glob, void*acc)
{
    fclaw2d_ray_acc_t* ray_acc = (fclaw2d_ray_acc_t*) acc;

    if (ray_acc->num_rays == 0)
    {
        return;
    }
    FCLAW_ASSERT(ray_acc->rays != NULL);

    /* We reset untrustworthy to 0 */
    int i, num_rays;
    fclaw2d_ray_t *ray;
    fclaw2d_ray_t *rays = ray_acc->rays;
    num_rays = ray_acc->num_rays;
    for (i =0; i < num_rays; i++)
    {
        ray = &rays[i];
        ray->untrustworthy = 0;
    }
}


static
void ray_integrate(fclaw_global_t *glob, void *acc)
{
    fclaw2d_ray_acc_t* ray_acc = (fclaw2d_ray_acc_t*) acc;

    if (ray_acc->num_rays == 0)
    {
        return;
    }
    FCLAW_ASSERT(ray_acc->rays != NULL);

    /* Copy arrays stored in accumulator to an sc_array */
    int i, num_rays;
    double intval;
    fclaw2d_ray_t *ray, *fclaw_ray;
    sc_array_t  *sc_rays = sc_array_new (sizeof (fclaw2d_ray_t));
    num_rays = ray_acc->num_rays;
    for(i = 0; i < num_rays; i++)
    {
        /* Add one ray to sc_rays;  return newly pushed array */
        fclaw_ray = (fclaw2d_ray_t *) sc_array_push (sc_rays);

        /* Set data for ray */
        fclaw_ray->num = ray_acc->rays[i].num;
        fclaw_ray->ray_data = (void*) ray_acc->rays[i].ray_data;
        fclaw_ray->untrustworthy = 0;
    }

    sc_array_t *integrals = sc_array_new_count (sizeof (double), num_rays);

    const fclaw2d_ray_vtable_t* ray_vt = fclaw2d_ray_vt(glob);

    /* This does a compute and a gather. */
    fclaw2d_domain_integrate_rays(glob->domain, ray_vt->integrate, 
                                  sc_rays, integrals, glob);

    /* Copy integral value back to fclaw2d_ray_t */
    fclaw2d_ray_t *rays = ray_acc->rays;
    for(i = 0; i < num_rays; i++)
    {
        /* Return values from integrals */
        intval = *((double*) sc_array_index_int(integrals,i));

        fclaw_ray = (fclaw2d_ray_t*)sc_array_index_int(sc_rays, i);

        /* Update rays stored in accumulator so we can report values later */
        ray = &rays[i];
        ray->integral = intval;
        ray->untrustworthy = fclaw_ray->untrustworthy;
    }

    sc_array_destroy (sc_rays);
    sc_array_destroy (integrals);
}

static
void ray_gather(fclaw_global_t *glob, void* acc, int init_flag)
{
    /* Here is where we would do an all-reduce, but again this
      is handled by fclaw2d_convenience routines */
    fclaw2d_ray_acc_t* ray_acc = (fclaw2d_ray_acc_t*) acc;
    if (ray_acc->num_rays == 0)
    {
        return;
    }

    FCLAW_ASSERT(ray_acc->rays != NULL);

    /* we check if the rays yielded trustworthy results on all processes */
    int i, *local_untrustworthy, *global_untrustworthy;
    fclaw2d_ray_t *ray;
    int num_rays = ray_acc->num_rays;
    fclaw2d_ray_t *rays = ray_acc->rays;

    local_untrustworthy = FCLAW_ALLOC (int, num_rays);
    for (i = 0; i < num_rays; i++)
    {
        ray = &rays[i];
        local_untrustworthy[i] = ray->untrustworthy;
    }

    /* aggregate the boolean untrustworthy condition for each ray */
    global_untrustworthy = FCLAW_ALLOC (int, num_rays);
    sc_MPI_Allreduce(local_untrustworthy, global_untrustworthy, num_rays,
                     sc_MPI_INT, sc_MPI_LOR, glob->mpicomm);

    for(int i = 0; i < num_rays; i++)
    {
        ray = &rays[i];
        int id = ray->num;
        double intval = ray->integral;
        if (global_untrustworthy[i])
        {
            fclaw_global_essentialf("The integral of ray %2d is not trustworthy.\n",id);
        }
        fclaw_global_essentialf("ray %2d; integral = %24.16e\n",id,intval);
    }

    FCLAW_FREE (local_untrustworthy);
    FCLAW_FREE (global_untrustworthy);
}


static
void ray_finalize(fclaw_global_t *glob, void** acc)
{
    fclaw2d_ray_acc_t* ray_acc = *((fclaw2d_ray_acc_t**) acc);
    if (ray_acc->rays != NULL)
    {
        ray_deallocate(glob,&ray_acc->rays,&ray_acc->num_rays);
    }
#if 0
    sc_array_destroy (ray_acc->rays);
    sc_array_destroy (ray_accc->integrals);
#endif    

    /* Matches allocation in ray_initialize */
    FCLAW_FREE(*acc);
    *acc = 0;
}


/* ---------------------------------- Virtual table  ---------------------------------- */

static
fclaw2d_ray_vtable_t* fclaw2d_ray_vt_new()
{
    return (fclaw2d_ray_vtable_t*) FCLAW_ALLOC_ZERO (fclaw2d_ray_vtable_t, 1);
}

static
void fclaw2d_ray_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

#if 0
static
fclaw2d_ray_vtable_t* fclaw2d_ray_vt_new()
{
    //FCLAW_ASSERT(s_rays_vt.is_set == 0);
    return &s_ray_vt;
}

fclaw2d_ray_vtable_t* fclaw2d_ray_vt()
{
    FCLAW_ASSERT(s_ray_vt.is_set != 0);
    return &s_ray_vt;
}
#endif

void fclaw2d_ray_vtable_initialize(fclaw_global_t *glob)
{
    fclaw2d_diagnostics_vtable_t * diag_vt = fclaw2d_diagnostics_vt(glob);
    diag_vt->ray_init_diagnostics     = ray_initialize;    
    diag_vt->ray_compute_diagnostics  = ray_integrate;
    diag_vt->ray_gather_diagnostics   = ray_gather;
    diag_vt->ray_reset_diagnostics    = ray_reset;
    diag_vt->ray_finalize_diagnostics = ray_finalize;

    fclaw2d_ray_vtable_t* rays_vt = fclaw2d_ray_vt_new(); 

    FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fclaw2d_rays") == NULL);
    fclaw_pointer_map_insert(glob->vtables, "fclaw2d_rays", rays_vt, fclaw2d_ray_vt_destroy);

    rays_vt->is_set = 1;
}


/* ---------------------------- Get Access Functions ---------------------------------- */


fclaw2d_ray_vtable_t* fclaw2d_ray_vt(fclaw_global_t* glob)
{
    fclaw2d_ray_vtable_t* ray_vt = (fclaw2d_ray_vtable_t*) 
                                fclaw_pointer_map_get(glob->vtables, "fclaw2d_rays");
    FCLAW_ASSERT(ray_vt != NULL);
    FCLAW_ASSERT(ray_vt->is_set != 0);
    return ray_vt;
}


fclaw2d_ray_t* fclaw2d_ray_allocate_rays(int num_rays)
{
    fclaw2d_ray_t* rays = FCLAW_ALLOC(fclaw2d_ray_t,num_rays);
    return rays;
}

int fclaw2d_ray_deallocate_rays(fclaw2d_ray_t **rays)
{
    FCLAW_ASSERT(*rays != NULL);
    FCLAW_FREE(*rays);
    *rays = NULL;
    return 0;  /* Number of rays */
}


/* Routines that operate on a single array */
void fclaw2d_ray_set_ray(fclaw2d_ray_t *r, 
                         int id, 
                         void* user_ray){
    /* User calls this to set ray data */
    r->num = id;
    r->ray_data = user_ray;
}



void* fclaw2d_ray_get_ray(fclaw2d_ray_t *r, 
                          int *id)
{
    *id = r->num;
    return r->ray_data;
}

