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

#include "fclaw_base.h"
#include <fclaw_gauges.h>

#include <fclaw_pointer_map.h>

#include <fclaw2d_options.h>
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
#endif

/* -------------------------------------------------------------------------------------*/

typedef struct fclaw_gauge_acc
{
    int num_gauges;
    int is_latest_domain;
    struct fclaw_gauge *gauges;
} fclaw_gauge_acc_t;


/* This is only used for blocks */
typedef struct fclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fclaw_gauge_info_t;

static void gauge_info_destroy(void* gauge_info)
{
    FCLAW_FREE(gauge_info);
}

static
void gauge_initialize(fclaw2d_global_t* glob, void** acc)
{
    const fclaw_options_t * fclaw_opt = fclaw2d_get_options(glob);

    /* ------------------------------------------------------------------
       These two calls are the only calls that should worry about the format
       GeoClaw of the files gauges.data (created with make_data.py) and 
       gauge output files (e.g. gauge00123.txt)
       ---------------------------------------------------------------- */

    fclaw_gauge_acc_t* gauge_acc = FCLAW_ALLOC(fclaw_gauge_acc_t,1);

    int num_gauges;
    fclaw_gauge_t *gauges;
    if (!fclaw_opt->output_gauges)
    {
        /* User does not want any gauge output, so no point in creating gauges */
        num_gauges = 0;
        gauge_acc->gauges = NULL;
    }
    else
    {
        fclaw_set_gauge_data(glob, &gauge_acc->gauges, &num_gauges);

        /* Set diagnostic accumulutor info  */
        //gauge_acc->gauges = gauges;  /* Might be NULL */    
    }
    *acc = gauge_acc;
    gauge_acc->num_gauges = num_gauges;


    fclaw_gauge_info_t* gauge_info = FCLAW_ALLOC_ZERO(fclaw_gauge_info_t,1);
    fclaw2d_global_attribute_store(glob, "gauge_info", gauge_info, gauge_info_destroy);

    if (num_gauges > 0)
    {
        gauges = gauge_acc->gauges;
        fclaw_create_gauge_files(glob,gauges,num_gauges);    

        /* ------------------------------------------------------------------
           Finish setting gauges with ForestClaw specific info 
           For  q_gauges, users must still allocate space for variables to be
           stored in the print buffer.
           ---------------------------------------------------------------- */
        int buffer_len = fclaw_opt->gauge_buffer_length;
        for(int i = 0; i < num_gauges; i++)
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

           Basic idea : Create a list of gauges for each block.

          -- For each block j, loop over all gauges.  If 
             gauge[i] is in block[j], add gauge[i] to 
             block[j] list of gauges. 

          ----------------------------------------------------- */

        fclaw2d_map_context_t* cont = glob->cont;

        int num_blocks = glob->domain->num_blocks;

        gauge_info->block_offsets = sc_array_new_count(sizeof(int), 
                                                      num_blocks+1);
        gauge_info->coordinates = sc_array_new_count(2*sizeof(double), num_gauges);

        /* We don't know how the blocks are arranged in the brick domain
           so we reverse engineer this information */
        double x0 = 0;
        double y0 = 0;
        double x1 = 1;
        double y1 = 1;

        int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);
        int mi = fclaw_opt->mi;
        int mj = fclaw_opt->mj;

        fclaw2d_block_t *blocks = glob->domain->blocks;

        int number_of_gauges_set = 0;
        int *bo = (int*) sc_array_index_int(gauge_info->block_offsets,0);
        bo[0] = 0;

        double xll,yll, xur, yur;
        for (int nb = 0; nb < num_blocks; nb++)
        {
            fclaw2d_block_t *block = &blocks[nb];
            if (is_brick)
            {
                /* Scale local block coordinates to global 
                [0,1]x[0,1] coordinates.

                Example : 4 x 1 arrangement of blocks 
                The point [0.5,0.5] in block 1  will map to
                the point (0.125,0.5) in [0,1]x[0,1].  

                We can't do this directly, because we don't know 
                (easily) how the blocks are numbered in the 
                brick grid */ 

                double z;
                fclaw2d_map_c2m_nomap_brick(cont,nb,x0,y0,&xll,&yll,&z);
                fclaw2d_map_c2m_nomap_brick(cont,nb,x1,y1,&xur,&yur,&z);                
            }
            else
            {
                /* eventually, we should set up the brick domain so that it 
                   defaults to these values */
                xll = x0;
                yll = y0;
                xur = x1;
                yur = y1;
            }
            for(int i = 0; i < num_gauges; i++)
            {
                /* Map gauge to global [0,1]x[0,1] space. This works for the brick
                   but not clear what happens for the cubed sphere */
                double x, y;
                fclaw_gauge_normalize_coordinates(glob,block,nb,&gauges[i],&x,&y);

                int gauge_in_block = (xll <= x && x <= xur && yll <= y && y <= yur);
                if (gauge_in_block)
                {
                    int ng = number_of_gauges_set;
                    gauges[i].blockno = nb;     /* gauge[i] is in block nb */

                    /* Store gauge coordinates in global brick domain.
                       This assumes that coordinates are in a global brick
                       domain in [0,mi]x[0,mj].  This may not work for the 
                       cubed sphere case */

                    double *c = (double*) sc_array_index_int(gauge_info->coordinates, ng);
                    c[0] = mi*(x - xll);
                    c[1] = mj*(y - yll);

                    /* Store location of this gauge in a global list */
                    gauges[i].location_in_results = ng;
                    number_of_gauges_set++;
                }
            }
            bo = (int*) sc_array_index_int(gauge_info->block_offsets, nb+1);
            bo[0] = number_of_gauges_set;        
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

    fclaw2d_diagnostics_accumulator_t* acc = fclaw2d_global_get_attribute(glob, "acc");
    fclaw_gauge_acc_t* gauge_acc = 
              (fclaw_gauge_acc_t*) acc->gauge_accumulator;
    fclaw_gauge_info_t* gauge_info = (fclaw_gauge_info_t*) 
              fclaw2d_global_get_attribute(glob, "gauge_info");

    /* Locate each gauge in the new mesh */
    num = gauge_acc->num_gauges;

    if (num == 0)
    {
        return;
    }

    sc_array_t *results = sc_array_new_size(sizeof(int), num);

    fclaw2d_domain_search_points(glob->domain, 
                                 gauge_info->block_offsets,
                                 gauge_info->coordinates, results);

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
    /* Clean up gauges and print anything left over in buffers */
    fclaw_gauge_acc_t* gauge_acc = *((fclaw_gauge_acc_t**) acc);
    fclaw_gauge_t *gauges = gauge_acc->gauges;
    fclaw_gauge_info_t* gauge_info = (fclaw_gauge_info_t*) 
              fclaw2d_global_get_attribute(glob, "gauge_info");

    for(int i = 0; i < gauge_acc->num_gauges; i++)
    {
        fclaw_gauge_t *g = &gauges[i];

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

    if (gauge_info->block_offsets != NULL)
    {
        sc_array_destroy(gauge_info->block_offsets);
    }

    if (gauge_info->coordinates != NULL)
    {
        sc_array_destroy(gauge_info->coordinates);
    }
    
    FCLAW_FREE(gauge_acc);
    FCLAW_FREE(gauge_info);
    *acc = NULL;    
}

/* ---------------------------------- Virtual table  ---------------------------------- */
static
fclaw_gauges_vtable_t* fclaw_gauges_vt_new()
{
    return (fclaw_gauges_vtable_t*) FCLAW_ALLOC_ZERO (fclaw_gauges_vtable_t, 1);
}

static
void fclaw_gauges_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fclaw_gauges_vtable_t* fclaw_gauges_vt(fclaw2d_global_t* glob)
{
	fclaw_gauges_vtable_t* gauges_vt = (fclaw_gauges_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fclaw_gauges");
	FCLAW_ASSERT(gauges_vt != NULL);
	FCLAW_ASSERT(gauges_vt->is_set != 0);

    return gauges_vt;
}

void fclaw_gauges_vtable_initialize(fclaw2d_global_t* glob)
{
    fclaw2d_diagnostics_vtable_t * diag_vt = fclaw2d_diagnostics_vt(glob);

    fclaw_gauges_vtable_t* gauges_vt = fclaw_gauges_vt_new();

    diag_vt->gauges_init_diagnostics     = gauge_initialize;
    diag_vt->gauges_compute_diagnostics  = gauge_update;
    diag_vt->gauges_finalize_diagnostics = gauge_finalize;

    gauges_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fclaw_gauges") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fclaw_gauges", gauges_vt, fclaw_gauges_vt_destroy);
}
/* ---------------------------- Virtualized Functions --------------------------------- */

void fclaw_set_gauge_data(fclaw2d_global_t* glob, 
                          fclaw_gauge_t **gauges, 
                          int *num_gauges)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
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
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->create_gauge_files != NULL);
    gauge_vt->create_gauge_files(glob, gauges, num_gauges);    

}

void fclaw_gauge_normalize_coordinates(fclaw2d_global_t *glob, 
                                      fclaw2d_block_t *block,
                                      int blockno, 
                                      fclaw_gauge_t *g,
                                      double *xc, double *yc)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->normalize_coordinates != NULL);
    gauge_vt->normalize_coordinates(glob, block,blockno,g,xc,yc);    
}


void  fclaw_update_gauge(fclaw2d_global_t* glob, 
                         fclaw2d_block_t *block,
                         fclaw2d_patch_t *patch,
                         int blockno, int patchno,
                         double tcurr, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->update_gauge != NULL);

    gauge_vt->update_gauge(glob,block,patch,blockno,patchno,tcurr,g);
}

void fclaw_print_gauge_buffer(fclaw2d_global_t* glob, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
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

int fclaw_gauge_get_id(fclaw2d_global_t *glob, 
                          fclaw_gauge_t *g)
{
    return g->num;
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

void fclaw_gauge_set_user_data(fclaw2d_global_t *glob,
                               fclaw_gauge_t* g,
                               void* user)
{
    g->user_data = user;
}

void* fclaw_gauge_get_user_data(fclaw2d_global_t *glob,
                                  fclaw_gauge_t* g)
{
    return g->user_data;
}


#ifdef __cplusplus
}
#endif
