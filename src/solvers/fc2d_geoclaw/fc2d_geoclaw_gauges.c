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

#include "fc2d_geoclaw_gauges.h"
#include "fc2d_geoclaw_gauges_default.h"

#include "fc2d_geoclaw_options.h"

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_fort.h"

#include "fclaw2d_options.h"
#include <fclaw2d_global.h>
#include <fclaw2d_convenience.h>  /* Needed to get search function for gauges */
#include <fclaw2d_diagnostics.h>

/* Some mapping functions */
#include <fclaw2d_map_brick.h>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* -------------------------------------------------------------------------------------*/

static fc2d_geoclaw_gauges_vtable_t s_geoclaw_gauges_vt;

typedef struct fc2d_geoclaw_gauge_acc
{
    int num_gauges;
    int is_latest_domain;
    struct fc2d_geoclaw_gauge *gauges;
} fc2d_geoclaw_gauge_acc_t;


typedef struct fc2d_geoclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fc2d_geoclaw_gauge_info_t;


static fc2d_geoclaw_gauge_info_t gauge_info;

static
void gauge_initialize(fclaw2d_global_t* glob, void** acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc;
    fc2d_geoclaw_gauge_t *gauges;
    int i, num_gauges;

    const fc2d_geoclaw_gauges_vtable_t* gauge_vt = fc2d_geoclaw_gauges_vt();

    const fclaw_options_t * gparms = fclaw2d_get_options(glob);
    const fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);

    /* ------------------------------------------------------------------
       These two calls are the only calls that should worry about the format
       GeoClaw of the files gauges.data (created with make_data.py) and 
       gauge output files (e.g. gauge00123.txt)
       ---------------------------------------------------------------- */
    
    gauge_vt->read_gauges_data(glob, &gauges, &num_gauges);

    /* Set diagnostic accumulutor info  */
    gauge_acc = FCLAW_ALLOC(fc2d_geoclaw_gauge_acc_t,1);
    *acc = gauge_acc;
    gauge_acc->num_gauges = num_gauges;
    gauge_acc->gauges = gauges;  /* Might be NULL */

    gauge_vt->create_gauge_files(glob, gauges, num_gauges);    

    /* ------------------------------------------------------------------
        Finish setting gauges with ForestClaw specific info 
        For  q_gauges, users must still allocate space for variables to be
        stored in the print_buffer.
        ---------------------------------------------------------------- */
    int bl = geo_opt->gauge_print_buffer_length;   
    for(i = 0; i < num_gauges; i++)
    {
        gauges[i].last_time = gauges[i].t1;
        gauges[i].patchno = -1;
        gauges[i].blockno = -1;
        gauges[i].location_in_results = -1;
        gauges[i].buffer = FCLAW_ALLOC(void*,bl);
        gauges[i].buffer_index = 0;
    }       

    /* -----------------------------------------------------
       Set up block offsets and coordinate list for p4est
       search function
       ----------------------------------------------------- */

    if (num_gauges == 0)
    {
        gauge_info.block_offsets = NULL;
        gauge_info.coordinates = NULL;
        return;
    }

    fclaw2d_map_context_t* cont = glob->cont;

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    gauge_info.block_offsets = sc_array_new_size(sizeof(int), 
                                                 glob->domain->num_blocks+1);
    gauge_info.coordinates = sc_array_new_size(2*sizeof(double), num_gauges);

    int *block_offsets = (int*) sc_array_index_int(gauge_info.block_offsets, 0);
    double *coordinates = (double*) sc_array_index_int(gauge_info.coordinates, 0);

    if (is_brick)
    {
        /* We don't know how the blocks are arranged in the brick domain
           so we reverse engineer this information
        */
        int nb,mi,mj;
        double x,y;
        double z;
        double xll,yll;
        double xur,yur;
        int number_of_gauges_set;


        double x0,y0,x1,y1;
        x0 = 0;
        y0 = 0;
        x1 = 1;
        y1 = 1;

        FCLAW2D_MAP_BRICK_GET_DIM(&cont,&mi,&mj);

        number_of_gauges_set = 0;
        block_offsets[0] = 0;
        for (nb = 0; nb < glob->domain->num_blocks; ++nb)
        {
            /* Scale to [0,1]x[0,1], based on blockno */
            fclaw2d_map_c2m_nomap_brick(cont,nb,x0,y0,&xll,&yll,&z);
            fclaw2d_map_c2m_nomap_brick(cont,nb,x1,y1,&xur,&yur,&z);
            for(int i = 0; i < num_gauges; i++)
            {
                /* Map gauge to global [0,1]x[0,1] space */
                x = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
                y = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
                if (xll <= x && x <= xur && yll <= y && y <= yur)
                {
                    int ng = number_of_gauges_set;
                    gauges[i].blockno = nb;
                    coordinates[2*ng] = mi*(x - xll);
                    coordinates[2*ng+1] = mj*(y - yll);
                    gauges[i].location_in_results = ng;
                    number_of_gauges_set++;
                }
            }
            block_offsets[nb+1] = number_of_gauges_set;
        }
    }
    else
    {
        for(int i = 0; i < num_gauges; i++)
        {
            gauges[i].blockno = 0;
            gauges[i].location_in_results = i;
        }

        block_offsets[0] = 0;
        block_offsets[1] = num_gauges;

        for (int i = 0; i < num_gauges; ++i)
        {
            coordinates[2*i] = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
            coordinates[2*i+1] = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
        }
    }
}


static
void gauge_update(fclaw2d_global_t *glob, void* acc)
{
    double tcurr;
    int num_gauges;

    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;
    fc2d_geoclaw_gauge_t *g;

    const fc2d_geoclaw_gauges_vtable_t* gauge_vt = fc2d_geoclaw_gauges_vt();

    fc2d_geoclaw_gauge_acc_t* gauge_acc = (fc2d_geoclaw_gauge_acc_t*) acc;
    fc2d_geoclaw_gauge_t *gauges = gauge_acc->gauges;

    const fc2d_geoclaw_options_t *geo_opt = fc2d_geoclaw_get_options(glob);

    tcurr = glob->curr_time;
    num_gauges = gauge_acc->num_gauges;

    for (int i = 0; i < num_gauges; i++)
    {
        g = &gauges[i];
        block = &glob->domain->blocks[g->blockno];
        if (g->patchno >= 0)
        {
            patch = &block->patches[g->patchno];
            FCLAW_ASSERT(patch != NULL);

            if (tcurr >= g->t1 && tcurr <= g->t2 &&
                tcurr - g->last_time >= g->min_time_increment)
            {
                gauge_vt->update_gauge(glob,block,patch,g->blockno,g->patchno,tcurr,g);
                g->buffer_index++;
                g->last_time = tcurr;
                
                if (g->buffer_index == geo_opt->gauge_print_buffer_length)
                {
                    /* This prints buffers and deletes buffer storage */
                    gauge_vt->print_gauge_buffer(glob,g);
                    g->buffer_index = 0;
                }  
            }
        }
    }
}


void fc2d_geoclaw_locate_gauges(fclaw2d_global_t *glob)
{
    int i,index;
    fc2d_geoclaw_gauge_acc_t* gauge_acc = 
              (fc2d_geoclaw_gauge_acc_t*) glob->acc->solver_accumulator;

    /* Locate each gauge in the new mesh */
    int num = gauge_acc->num_gauges;

    if (num == 0)
    {
        return;
    }

    sc_array_t *results = sc_array_new_size(sizeof(int), num);
    fclaw2d_domain_search_points(glob->domain, gauge_info.block_offsets,
                                 gauge_info.coordinates, results);

    for (i = 0; i < gauge_acc->num_gauges; ++i)
    {
        index = gauge_acc->gauges[i].location_in_results;
        FCLAW_ASSERT(index >= 0 && index < num);

        /* patchno == -1  : Patch is not on this processor
           patchno >= 0   : Patch number is local patch list.
        */

        int patchno = *((int *) sc_array_index_int(results, index));
        gauge_acc->gauges[i].patchno = patchno;
    }
    sc_array_destroy(results);
}

static
void gauge_finalize(fclaw2d_global_t *glob, void** acc)
{
    int i;
    fc2d_geoclaw_gauge_t *g;

    /* Clean up gauges and print anything left over in buffers */
    fc2d_geoclaw_gauge_acc_t* gauge_acc = *((fc2d_geoclaw_gauge_acc_t**) acc);
    const fc2d_geoclaw_gauges_vtable_t* gauge_vt = fc2d_geoclaw_gauges_vt();

    fc2d_geoclaw_gauge_t *gauges = gauge_acc->gauges;

    for(i = 0; i < gauge_acc->num_gauges; i++)
    {
        g = &gauges[i];
        gauge_vt->print_gauge_buffer(glob,g);  
        gauge_vt->destroy_gauge(glob,g);   
    }

    if (gauge_acc->gauges != NULL)
    {
        FCLAW_FREE(gauge_acc->gauges); 
    }       

    if (gauge_info.block_offsets != NULL)
    {
        sc_array_destroy(gauge_info.block_offsets);
    }

    if (gauge_info.coordinates != NULL)
    {
        sc_array_destroy(gauge_info.coordinates);
    }
    
    FCLAW_FREE(gauge_acc);
    *acc = NULL;    
}

/* -------------------------- Virtual table  ---------------------------- */

static
fc2d_geoclaw_gauges_vtable_t* fc2d_geoclaw_gauges_vt_init()
{
    FCLAW_ASSERT(s_geoclaw_gauges_vt.is_set == 0);
    return &s_geoclaw_gauges_vt;
}

fc2d_geoclaw_gauges_vtable_t* fc2d_geoclaw_gauges_vt()
{
    FCLAW_ASSERT(s_geoclaw_gauges_vt.is_set != 0);
    return &s_geoclaw_gauges_vt;
}

void fc2d_geoclaw_gauges_vtable_set()
{
    fclaw2d_diagnostics_vtable_t * diag_vt = fclaw2d_diagnostics_vt();
    fc2d_geoclaw_gauges_vtable_t* gauges_vt = fc2d_geoclaw_gauges_vt_init();

    gauges_vt->read_gauges_data   = geoclaw_read_gauges_data_default;
    gauges_vt->create_gauge_files = geoclaw_create_gauge_files_default; 
    gauges_vt->destroy_gauge      = geoclaw_gauge_destroy_default; 

    gauges_vt->update_gauge       = geoclaw_gauge_update_default;
    gauges_vt->print_gauge_buffer = geoclaw_print_gauges_default;

    diag_vt->solver_init_diagnostics     = gauge_initialize;
    diag_vt->solver_compute_diagnostics  = gauge_update;
    diag_vt->solver_finalize_diagnostics = gauge_finalize;

    gauges_vt->is_set = 1;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
