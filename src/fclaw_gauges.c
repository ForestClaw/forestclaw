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

#include <fclaw_gauges.h>

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

static fclaw_gauges_vtable_t s_gauges_vt;

typedef struct fclaw_gauge_acc
{
    int num_gauges;
    int is_latest_domain;
    struct fclaw_gauge *gauges;
} fclaw_gauge_acc_t;


typedef struct fclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fclaw_gauge_info_t;


static fclaw_gauge_info_t gauge_info;

static
void gauge_initialize(fclaw2d_global_t* glob, void** acc)
{
    fclaw_gauge_acc_t* gauge_acc;
    fclaw_gauge_t *gauges;
    int i, num_gauges;

    const fclaw_options_t * fclaw_opt = fclaw2d_get_options(glob);

    /* ------------------------------------------------------------------
       These two calls are the only calls that should worry about the format
       GeoClaw of the files gauges.data (created with make_data.py) and 
       gauge output files (e.g. gauge00123.txt)
       ---------------------------------------------------------------- */
    
    fclaw_set_gauge_data(glob, &gauges, &num_gauges);

    /* Set diagnostic accumulutor info  */
    gauge_acc = FCLAW_ALLOC(fclaw_gauge_acc_t,1);
    *acc = gauge_acc;
    gauge_acc->num_gauges = num_gauges;
    gauge_acc->gauges = gauges;  /* Might be NULL */

    if (num_gauges == 0)
    {
        // gauge_info.block_offsets = NULL;
        // gauge_info.coordinates = NULL;
        return;
    }

    fclaw_create_gauge_files(glob,gauges,num_gauges);    

    /* ------------------------------------------------------------------
        Finish setting gauges with ForestClaw specific info 
        For  q_gauges, users must still allocate space for variables to be
        stored in the print buffer.
        ---------------------------------------------------------------- */
    int buffer_len = fclaw_opt->gauge_buffer_length;
    for(i = 0; i < num_gauges; i++)
    {
        gauges[i].last_time = gauges[i].t1;
        gauges[i].patchno = -1;
        gauges[i].blockno = -1;
        gauges[i].location_in_results = -1;
        gauges[i].buffer = FCLAW_ALLOC(void*,buffer_len);  /* Array of generic ptrs */
        gauges[i].next_buffer_location = 0;
    }       

    /* -----------------------------------------------------
       Set up block offsets and coordinate list for p4est
       search function
       ----------------------------------------------------- */

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
            for(i = 0; i < num_gauges; i++)
            {
                /* Map gauge to global [0,1]x[0,1] space */
                x = (gauges[i].xc - fclaw_opt->ax)/(fclaw_opt->bx-fclaw_opt->ax);
                y = (gauges[i].yc - fclaw_opt->ay)/(fclaw_opt->by-fclaw_opt->ay);
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
        for(i = 0; i < num_gauges; i++)
        {
            gauges[i].blockno = 0;
            gauges[i].location_in_results = i;
        }

        block_offsets[0] = 0;
        block_offsets[1] = num_gauges;

        for (i = 0; i < num_gauges; ++i)
        {
            coordinates[2*i] = (gauges[i].xc - fclaw_opt->ax)/(fclaw_opt->bx-fclaw_opt->ax);
            coordinates[2*i+1] = (gauges[i].yc - fclaw_opt->ay)/(fclaw_opt->by-fclaw_opt->ay);
        }
    }
}


static
void gauge_update(fclaw2d_global_t *glob, void* acc)
{
    double tcurr;
    int i, num_gauges;

    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;
    fclaw_gauge_t *g;

    const fclaw_options_t * fclaw_opt = fclaw2d_get_options(glob);

    fclaw_gauge_acc_t* gauge_acc = (fclaw_gauge_acc_t*) acc;
    fclaw_gauge_t *gauges = gauge_acc->gauges;

    int buffer_len = fclaw_opt->gauge_buffer_length;
    tcurr = glob->curr_time;
    num_gauges = gauge_acc->num_gauges;

    for (i = 0; i < num_gauges; i++)
    {
        g = &gauges[i];
        if (tcurr >= g->t1 && tcurr <= g->t2 &&
            tcurr - g->last_time >= g->min_time_increment)
        {
            /* Update the last time, even though this gauge may not be local to
               this processor. This keeps the time consistent across all processors, 
               so that when this gauge is local to this processor, it knows when 
               it was last updated (even if it was updated on another processor). */
            g->last_time = tcurr;            

            if (g->is_local)
            {
                block = &glob->domain->blocks[g->blockno];
                patch = &block->patches[g->patchno]; 
                fclaw_update_gauge(glob,block,patch,
                                          g->blockno,g->patchno,
                                          tcurr,g);

                g->next_buffer_location++;
                
                if (g->next_buffer_location == buffer_len)
                {
                    fclaw_print_gauge_buffer(glob,g);
                    g->next_buffer_location = 0;
                }  
            }
            else
            {
                /* If this gauge is not local, then it should not have anything in the 
                   buffer */
                FCLAW_ASSERT(g->next_buffer_location == 0);
            }
        }
    }
}


void fclaw_locate_gauges(fclaw2d_global_t *glob)
{
    int i,index,num;
    fclaw_gauge_t *g;

    fclaw_gauge_acc_t* gauge_acc = 
              (fclaw_gauge_acc_t*) glob->acc->solver_accumulator;

    /* Locate each gauge in the new mesh */
    num = gauge_acc->num_gauges;

    if (num == 0)
    {
        return;
    }

    sc_array_t *results = sc_array_new_size(sizeof(int), num);

    fclaw2d_domain_search_points(glob->domain, 
                                 gauge_info.block_offsets,
                                 gauge_info.coordinates, results);

    for (i = 0; i < gauge_acc->num_gauges; ++i)
    {
        g = &gauge_acc->gauges[i];

        index = g->location_in_results;
        FCLAW_ASSERT(index >= 0 && index < num);

        /* patchno == -1  : Patch is not on this processor
           patchno >= 0   : Patch number is in local patch list.
        */

        /* Current patch no (patches can move under gauges, but blocks 
           remain fixed. */
        g->patchno = *((int *) sc_array_index_int(results, index));
        g->is_local = (g->patchno >= 0);  /* Local to this processor */
        if (!g->is_local && g->next_buffer_location > 0)
        {
            /* Patch moved off of processor, but the buffer is not empty. */
            fclaw_print_gauge_buffer(glob,g);
            g->next_buffer_location = 0;
        }
    }
    sc_array_destroy(results);
}

static
void gauge_finalize(fclaw2d_global_t *glob, void** acc)
{
    int i;
    fclaw_gauge_t *g;

    /* Clean up gauges and print anything left over in buffers */
    fclaw_gauge_acc_t* gauge_acc = *((fclaw_gauge_acc_t**) acc);
    fclaw_gauge_t *gauges = gauge_acc->gauges;

    for(i = 0; i < gauge_acc->num_gauges; i++)
    {
        g = &gauges[i];

        /* Every processor owns every gauge (which will scale up to a few 
        hundred gauges).  But we only want to print those gauge buffers that 
        for gauges that are on the local processor */        
        if (g->is_local)
        {
            fclaw_print_gauge_buffer(glob,g);
        }
        FCLAW_FREE(g->buffer);               
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

/* ---------------------------------- Virtual table  ---------------------------------- */

static
fclaw_gauges_vtable_t* fclaw_gauges_vt_init()
{
    FCLAW_ASSERT(s_gauges_vt.is_set == 0);
    return &s_gauges_vt;
}

fclaw_gauges_vtable_t* fclaw_gauges_vt()
{
    FCLAW_ASSERT(s_gauges_vt.is_set != 0);
    return &s_gauges_vt;
}

void fclaw_gauges_vtable_initialize()
{
    fclaw2d_diagnostics_vtable_t * diag_vt = fclaw2d_diagnostics_vt();

    fclaw_gauges_vtable_t* gauges_vt = fclaw_gauges_vt_init();

    diag_vt->solver_init_diagnostics     = gauge_initialize;
    diag_vt->solver_compute_diagnostics  = gauge_update;
    diag_vt->solver_finalize_diagnostics = gauge_finalize;

    gauges_vt->is_set = 1;
}
/* ---------------------------- Virtualized Functions --------------------------------- */

void fclaw_set_gauge_data(fclaw2d_global_t* glob, 
                          fclaw_gauge_t **gauges, 
                          int *num_gauges)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt();
    if (gauge_vt->set_gauge_data == NULL)
    {
        *gauges = NULL;
        *num_gauges = 0;
    }
    else
    {
        gauge_vt->set_gauge_data(glob, gauges, num_gauges);    
    }
}

void fclaw_create_gauge_files(fclaw2d_global_t* glob, 
                              fclaw_gauge_t *gauges, 
                              int num_gauges)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt();
    FCLAW_ASSERT(gauge_vt->create_gauge_files != NULL);
    gauge_vt->create_gauge_files(glob, gauges, num_gauges);    

}

void  fclaw_update_gauge(fclaw2d_global_t* glob, 
                         fclaw2d_block_t *block,
                         fclaw2d_patch_t *patch,
                         int blockno, int patchno,
                         double tcurr, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt();
    FCLAW_ASSERT(gauge_vt->update_gauge != NULL);

    gauge_vt->update_gauge(glob,block,patch,blockno,patchno,tcurr,g);
}

void fclaw_print_gauge_buffer(fclaw2d_global_t* glob, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt();
    FCLAW_ASSERT(gauge_vt->print_gauge_buffer != NULL);

    gauge_vt->print_gauge_buffer(glob,g);
}

/* ---------------------------- Get Access Functions ---------------------------------- */

void fclaw_gauge_allocate(fclaw2d_global_t *glob, int num_gauges,
                          fclaw_gauge_t **g)
{
    *g = (fclaw_gauge_t*) FCLAW_ALLOC(fclaw_gauge_t,num_gauges);
}

void fclaw_gauge_set_data(fclaw2d_global_t *glob, 
                          fclaw_gauge_t *g,
                          int num, double xc, double yc, 
                          double  t1, double t2, 
                          double min_time_increment)
{
    g->num = num;
    g->xc = xc;
    g->yc = yc;
    g->t1 = t1;
    g->t2 = t2;
    g->min_time_increment = min_time_increment;
}



void fclaw_gauge_get_data(fclaw2d_global_t *glob, 
                          fclaw_gauge_t *g,
                          int *num, 
                          double *xc, double *yc, 
                          double  *t1, double *t2)
{
    *num = g->num;
    *xc = g->xc;
    *yc = g->yc;
    *t1 = g->t1;
    *t2 = g->t2;
}

void fclaw_gauge_get_buffer(fclaw2d_global_t *glob,
                            fclaw_gauge_t *g,
                            int *kmax, void*** gauge_buffer)
{
    *kmax = g->next_buffer_location;
    *gauge_buffer = g->buffer;
}

void fclaw_gauge_set_buffer_entry(fclaw2d_global_t *glob,
                                  fclaw_gauge_t* g,
                                  void* guser)
{
    int k = g->next_buffer_location;
    g->buffer[k] = guser;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
