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

#include <fclaw2d_domain.h>
#include <fclaw2d_patch.hpp>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


void fclaw2d_domain_data_new(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t* ddata = (fclaw2d_domain_data_t*) domain->user;
    ddata = FCLAW2D_ALLOC_ZERO(fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;

    ddata->count_set_clawpatch = ddata->count_delete_clawpatch = 0;
    ddata->count_amr_advance = 0;
    ddata->count_ghost_exchange = 0;
    ddata->count_amr_regrid = 0;
    ddata->is_latest_domain = 0;        /* set 1 by amrinit or rebuild_domain */

    ddata->domain_exchange = NULL;

    ddata->curr_time = 0;
}

void fclaw2d_domain_data_delete(fclaw2d_domain_t* domain)
{
    fclaw2d_domain_data_t* ddata = (fclaw2d_domain_data_t*) domain->user;

    FCLAW2D_FREE (ddata);
    domain->user = NULL;
}

fclaw2d_domain_data_t *fclaw2d_domain_get_data(fclaw2d_domain_t *domain)
{
    return (fclaw2d_domain_data_t *) domain->user;
}


void fclaw2d_domain_data_copy(fclaw2d_domain_t *old_domain, fclaw2d_domain_t *new_domain)
{
    fclaw2d_domain_data_t *ddata_old = fclaw2d_domain_get_data(old_domain);

    /* Has the data already been allocated? */
    fclaw2d_domain_data_t *ddata_new = fclaw2d_domain_get_data(new_domain);


    /* Move timers over to the new domain */
    ddata_old->is_latest_domain = 0;
    memcpy (ddata_new->timers, ddata_old->timers,
            sizeof (fclaw2d_timer_t) * FCLAW2D_TIMER_COUNT);
    ddata_new->is_latest_domain = 1;
    ddata_new->count_amr_advance = ddata_old->count_amr_advance;
    ddata_new->count_ghost_exchange = ddata_old->count_ghost_exchange;
    ddata_new->count_amr_regrid = ddata_old->count_amr_regrid;


    ddata_new->curr_time = ddata_old->curr_time;

}

void fclaw2d_domain_reset(fclaw2d_domain_t** domain)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (*domain);

    for(int i = 0; i < (*domain)->num_blocks; i++)
    {
        fclaw2d_block_t *block = (*domain)->blocks + i;
        fclaw2d_block_data_t *bd = (fclaw2d_block_data_t *) block->user;

        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_delete(*domain,patch);
        }

        FCLAW2D_FREE (bd);
        block->user = NULL;

    }

    fclaw2d_partition_delete(domain);

    /* Output memory discrepancy for the ClawPatch */
    if (ddata->count_set_clawpatch != ddata->count_delete_clawpatch)
    {
        printf ("[%d] This domain had Clawpatch set %d and deleted %d times\n",
                (*domain)->mpirank,
                ddata->count_set_clawpatch, ddata->count_delete_clawpatch);
    }

#if 0
    /* TODO : Do we need this here?  Will we ever want to report stats during
       a run? */
    if (ddata->is_latest_domain)  /* ddata->is_last_domain?  */
    {
        fclaw2d_timer_report(*domain);
    }
#endif

    fclaw2d_domain_data_delete(*domain);  // Delete allocated pointers to set of functions.

    fclaw2d_domain_destroy(*domain);
    *domain = NULL;
}


fclaw_app_t* fclaw2d_domain_get_app(fclaw2d_domain_t* domain)
{
    fclaw_app_t *app;

    app = (fclaw_app_t*)
          fclaw2d_domain_attribute_access(domain,"fclaw_app",NULL);

    FCLAW_ASSERT(app != NULL);
    return app;
}

void fclaw2d_domain_set_app(fclaw2d_domain_t* domain,fclaw_app_t* app)
{
    FCLAW_ASSERT(app != NULL);
    fclaw2d_domain_attribute_add (domain,"fclaw_app",app);
}

int fclaw2d_domain_get_num_patches(fclaw2d_domain_t* domain)
{
    return domain->global_num_patches;
}


void fclaw2d_domain_set_time(fclaw2d_domain_t *domain, double time)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    ddata->curr_time = time;
}

double fclaw2d_domain_get_time(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata = fclaw2d_domain_get_data (domain);
    return ddata->curr_time;
}

const amr_options_t* fclaw2d_forestclaw_get_options(fclaw2d_domain_t *domain)
{
    const amr_options_t *gparms;
    fclaw_app_t *app;

    app = fclaw2d_domain_get_app(domain);
    gparms = fclaw_forestclaw_get_options(app);
    return gparms;
}

void* fclaw2d_domain_get_user_options(fclaw2d_domain_t* domain)
{
    fclaw_app_t *app;

    app = fclaw2d_domain_get_app(domain);
    return fclaw_app_get_user(app);
}


const amr_options_t* get_domain_parms(fclaw2d_domain_t *domain)
{
    return fclaw2d_forestclaw_get_options(domain);
}

fclaw2d_map_context_t* fclaw2d_domain_get_map_context(fclaw2d_domain_t* domain)
{
    fclaw2d_map_context_t* cont;
  cont = (fclaw2d_map_context_t*)
         fclaw2d_domain_attribute_access (domain, "fclaw_map_context", NULL);
  FCLAW_ASSERT (cont != NULL);
  return cont;
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
