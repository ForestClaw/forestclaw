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
#include "fclaw2d_global.h"
#include "fclaw2d_options.h"
#include "fclaw2d_diagnostics.h"

static fclaw2d_ray_vtable_t s_ray_vt;


typedef struct fclaw2d_ray_acc
{
    int num_rays;
    int is_latest_domain;
    fclaw2d_ray_t *rays;
} fclaw2d_ray_acc_t;



/* ------------------------------ Virtualized functions ------------------------------- */

/* Do I need this level of indirection? */
static
void ray_allocate_and_define(fclaw2d_global_t* glob, 
                             fclaw2d_ray_t **rays, 
                             int *num_rays)
{
    const fclaw2d_ray_vtable_t* ray_vt = fclaw2d_ray_vt();

    FCLAW_ASSERT(ray_vt->allocate_and_define != NULL);
    ray_vt->allocate_and_define(glob, rays, num_rays);    
}

static
void ray_deallocate(fclaw2d_global_t* glob, 
                    fclaw2d_ray_t **rays, 
                    int *num_rays)
{
    const fclaw2d_ray_vtable_t* ray_vt = fclaw2d_ray_vt();

    FCLAW_ASSERT(ray_vt->deallocate != NULL);
    ray_vt->deallocate(glob, rays, num_rays);    
}




static
void ray_initialize(fclaw2d_global_t* glob, void** acc)
{
    fclaw2d_ray_acc_t *ray_acc = (fclaw2d_ray_acc_t*) FCLAW_ALLOC(fclaw2d_ray_acc_t,1);

    /* Check to see if user wants ray output */
    int num_rays;
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
    }
    *acc = ray_acc;
    ray_acc->num_rays = num_rays;


    if (num_rays > 0)
    {
        /* Things that could be done here : 
           1. Open ray files to write out time series data for each ray 
           2. Set integral to 0? 
        */
        //fclaw2d_ray_t *rays = ray_acc->rays;
        for(int i = 0; i < num_rays; i++)
        {
            //fclaw2d_ray_t *r = &rays[i];
            fclaw_global_essentialf("ray_initialize : Setting up ray %d : \n",i);
            /* do something? */
        }
    }
}

#if 0
static
void rays_reset(fclaw2d_global_t *glob, void**acc)
{
    /* Do nothing for now */
}
#endif


static
void ray_integrate(fclaw2d_global_t *glob, void *acc)
{
    /* Do nothing for now */
}

static
void ray_finalize(fclaw2d_global_t *glob, void** acc)
{
    fclaw2d_ray_acc_t* ray_acc = *((fclaw2d_ray_acc_t**) acc);
    FCLAW_ASSERT(ray_acc->rays != NULL);

    ray_deallocate(glob,&ray_acc->rays,&ray_acc->num_rays);

    FCLAW_FREE(*acc);
    *acc = 0;
}


/* ---------------------------------- Virtual table  ---------------------------------- */

static
fclaw2d_ray_vtable_t* fclaw2d_ray_vt_init()
{
    //FCLAW_ASSERT(s_rays_vt.is_set == 0);
    return &s_ray_vt;
}

fclaw2d_ray_vtable_t* fclaw2d_ray_vt()
{
    FCLAW_ASSERT(s_ray_vt.is_set != 0);
    return &s_ray_vt;
}

void fclaw2d_ray_vtable_initialize()
{

    fclaw2d_diagnostics_vtable_t * diag_vt = fclaw2d_diagnostics_vt();
    fclaw2d_ray_vtable_t* rays_vt = fclaw2d_ray_vt_init(); 

    diag_vt->ray_init_diagnostics     = ray_initialize;    
    diag_vt->ray_compute_diagnostics  = ray_integrate;
    diag_vt->ray_gather_diagnostics  = NULL; /* Gather is handled in compute step */
    diag_vt->ray_finalize_diagnostics = ray_finalize;

    rays_vt->is_set = 1;
}


/* ---------------------------- Get Access Functions ---------------------------------- */


/* Routines that operate on a single array */
void fclaw2d_ray_set_ray(fclaw2d_global_t *glob, 
                         fclaw2d_ray_t *r,
                         int id, 
                         void* ray_data)
{
    /* User calls this to set ray data */
    r->num = id;
    r->ray_data = ray_data;
}



void* fclaw2d_ray_get_ray(fclaw2d_global_t *glob, 
                          fclaw2d_ray_t *r,
                          int *id)
{
    *id = r->num;
    return r->ray_data;
}


#ifdef __cplusplus
}
#endif



