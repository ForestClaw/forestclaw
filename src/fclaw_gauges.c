/*
Copyright (c) 2012-2025 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_options.h>
#include <fclaw_global.h>
#include <fclaw_convenience.h>  /* Needed to get search function for gauges */
#include <fclaw_diagnostics.h>

/* Some mapping functions */
#include <fclaw_map_brick.h>
#include <fclaw_map.h>
#include <fclaw_map_query.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif    

/* -------------------------------------------------------------------------------------*/

typedef struct fclaw_gauge_acc
{
    int dim;
    int num_gauges;
    int num_gauges_set;  /* In case some gauges are not in the domain */    
    int is_latest_domain;
    fclaw_gauge_t *gauges;
} fclaw_gauge_acc_t;


/* This is only used for blocks */
typedef struct fclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
    int is_initialized;
    int moving_gauges;
} fclaw_gauge_info_t;

/* ----------------------------------------------------------------------------
    These five routines call virtualized gauge functions that are all 
    defined elsewhere (defaults are supplied in `clawpatch` code).  

    The five routines are 

        -- void gauges_read_data(...)
            -- Reads data from a "gauge.data" file. 

        -- void gauges_create_files(...)
            -- Creates gauge output files

        -- void gauges_normalize_coordinates(...)
            -- this is used to locate the gauge within a block.  This may 
               be replaced by a higher level "gauge_in_block" routine which can 
               be used for cubed-sphere and other non-Cartesian mappings. 

        -- void gauges_update(...)
            -- Update the gauge value.  This requires interpolation from the 
               mesh to the gauge point.

        -- void gauges_print_buffer(...)
            -- Print out the gauge buffer to a gauge file.

    These all call virtualized functions which are set in fclaw_clawpatch.cpp

        // In fc2d_clawpatch.cpp (vtable_initialize)
        fclaw_gauges_vtable_t*  gauges_vt = fclaw_gauges_vt(glob);
        gauges_vt->read_data              = fclaw_clawpatch_gauges_read_data;
        gauges_vt->create_files           = fclaw_clawpatch_gauges_create_files; 
        gauges_vt->normalize_coordinates  = fclaw_clawpatch_gauges_normalize_coordinates;
        gauges_vt->update                 = fclaw_clawpatch_gauges_update;
        gauges_vt->print_buffer           = fclaw_clawpatch_gauges_print;

*/ 


static
void gauges_read_data(fclaw_global_t* glob, 
                     fclaw_gauge_t **gauges, 
                     int *num_gauges,
                     int *dim)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    if (gauge_vt->read_data == NULL)
    {
        *gauges = NULL;
        *num_gauges = 0;
        *dim = 0;
    }
    else
    {
        gauge_vt->read_data(glob, gauges, num_gauges,dim);  
    }
}


static 
void gauges_create_files(fclaw_global_t* glob, 
                        fclaw_gauge_t *gauges, 
                        int num_gauges)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->create_files != NULL);
    gauge_vt->create_files(glob, gauges, num_gauges);    

}

static
void gauges_normalize_coordinates(fclaw_global_t *glob, 
                                      fclaw_block_t *block,
                                      int blockno, 
                                      fclaw_gauge_t *g,
                                      double *xc, double *yc, double *zc)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->normalize_coordinates != NULL);
    gauge_vt->normalize_coordinates(glob, block,blockno,g,xc,yc,zc);    
}


static
void  gauges_update(fclaw_global_t* glob, 
                    fclaw_block_t *block,
                    fclaw_patch_t *patch,
                    int blockno, int patchno,
                    double tcurr, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->update != NULL);

    gauge_vt->update(glob,block,patch,blockno,patchno,tcurr,g);
}

static
void gauges_print_buffer(fclaw_global_t* glob, fclaw_gauge_t *g)
{
    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->print_buffer != NULL);

    gauge_vt->print_buffer(glob,g);
}


/* ---------------------------------------------------------------------- 
    The routines below are used to set diagnostics routines

    Virtualized diagnostic routines : 

        -- void fclaw_gauges_initialize(fclaw_global_t *glob, void** acc)

        -- void fclaw_gauges_update(glob, acc)

        -- void fclaw_gauges_finalize(glob, acc)

    These are set in the gauges vtable, below.


*/    


/* Function needed for diagnostics */
static
void gauges_initialize(fclaw_global_t* glob, void** acc)
{
    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);

    /* ------------------------------------------------------------------
       These two calls are the only calls that should worry about the format
       GeoClaw of the files gauges.data (created with make_data.py) and 
       gauge output files (e.g. gauge00123.txt)
       ---------------------------------------------------------------- */
 
    fclaw_gauge_acc_t* gauge_acc = FCLAW_ALLOC(fclaw_gauge_acc_t,1);

    int num_gauges, gauge_dim;
    if (!fclaw_opt->output_gauges)
    {
        /* User does not want any gauge output, so no point in creating gauges */
        num_gauges = 0;
        gauge_acc->gauges = NULL;
        gauge_dim = -1;
    }
    else
    {
        /* 
            Read custom gauges file (e.g. gauges.data).  After this call, all
            gauges have current data 
        */
        gauges_read_data(glob, &gauge_acc->gauges, &num_gauges, &gauge_dim);
    }
    *acc = gauge_acc;
    gauge_acc->num_gauges = num_gauges;
    gauge_acc->dim = gauge_dim;


    fclaw_gauge_info_t *gauge_info = FCLAW_ALLOC_ZERO(fclaw_gauge_info_t,1);
    gauge_info->moving_gauges = fclaw_opt->moving_gauges;
    gauge_info->is_initialized = 0;

    fclaw_global_attribute_store(glob,
                                 "gauge_info",
                                 gauge_info,
                                 NULL,
                                 NULL);

    if (num_gauges > 0)
    {
        fclaw_gauge_t *gauges = gauge_acc->gauges;
        gauges_create_files(glob,gauges,num_gauges);    

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
            gauges[i].in_domain = 0;
            gauges[i].is_moving = gauge_info->moving_gauges; /* Assume static gauge */
            gauges[i].buffer = FCLAW_ALLOC(void*,buffer_len);  /* Array of generic ptrs */
            gauges[i].next_buffer_location = 0;
        }       
    }

    /* Create block lists and locate gauges in patches. 
       For static gauges, we only need to create the block list once. 
       For moving gauges, the block lists must be recreated each time 
       we want to locate gauges in patches. */

    fclaw_gauges_locate_patches(glob);
}    

static
int fclaw_gauges_setup_block_lists(fclaw_global_t* glob, 
                                    fclaw_gauge_acc_t *gauge_acc,
                                    sc_array_t **block_offsets_ptr,
                                    sc_array_t **coordinates_ptr)
{
    fclaw_gauge_t *gauges = gauge_acc->gauges;

    int num_gauges = gauge_acc->num_gauges;
    int gauge_dim = gauge_acc->dim;

    if (num_gauges > 0)
    {
        /* -----------------------------------------------------
           Set up block offsets and coordinate list for p4est
           search function

           Basic idea : Create a list of gauges for each block.

          -- For each block j, loop over all gauges.  If 
             gauge[i] is in block[j], add gauge[i] to 
             block[j] list of gauges. 

          ----------------------------------------------------- */

        fclaw_map_context_t* cont = glob->cont;

        int num_blocks = glob->domain->num_blocks;

        sc_array_t *block_offsets = sc_array_new_count(sizeof(int), 
                                                      num_blocks+1);

        sc_array_t *coordinates = sc_array_new_count(gauge_dim*sizeof(double), 
                                                     num_gauges);

        /* We don't know how the blocks are arranged in the brick domain
           so we reverse engineer this information */
        int is_brick = FCLAW_MAP_IS_BRICK(&cont);
        //int brick_dim[gauge_dim];

        fclaw_block_t *blocks = glob->domain->blocks;

        int num_gauges_set = 0;
        int *bo = (int*) sc_array_index_int(block_offsets,0);
        bo[0] = 0;

        double p0[gauge_dim];
        double p1[gauge_dim];

        /* 
            Loop over each block;  search all gauges to build list of gauges for 
            each block 
        */

        /* needed to compute coordinates of gauges */
        const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);
        int mi = fclaw_opt->mi;
        int mj = fclaw_opt->mj;

        for (int nb = 0; nb < num_blocks; nb++)
        {
            double xll,yll,zll,xur,yur,zur;

            fclaw_block_t *block = &blocks[nb];
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

                /* 
                    These mappings map local block coordinates (in block 'nb') to 
                    global coordinates.
                */

                p0[0] = 0;
                p0[1] = 0;
                p1[0] = 1;
                p1[1] = 1;                
                if (gauge_dim == 2)
                {
                    fclaw_map_2d_c2m_nomap_brick(cont,nb,p0[0],p0[1],&xll,&yll,&zll);
                    fclaw_map_2d_c2m_nomap_brick(cont,nb,p1[0],p1[1],&xur,&yur,&zur);
                }
                else if (gauge_dim == 3)
                {
                    p0[2] = 0;
                    p1[2] = 1;
                    fclaw_map_3d_c2m_nomap_brick(cont,nb,p0[0],p0[1],p0[2],&xll,&yll,&zll);
                    fclaw_map_3d_c2m_nomap_brick(cont,nb,p1[0],p1[1],p1[2],&xur,&yur,&zur); 
                }
                else
                {
                    fclaw_global_essentialf("gauge_initialize (fclaw_gauges.c) : " \
                                            "dim not set properly; dim=%d\n",gauge_dim);
                    exit(0);
                }
            }
            else
            {
                /* eventually, we should set up the brick domain so that it 
                   defaults to these values */
                xll = 0;
                yll = 0;
                zll = 0;
                xur = 1;
                yur = 1;
                zur = 1;
            }

            for(int i = 0; i < num_gauges; i++)
            {

                fclaw_gauge_t *g = &gauges[i];

                /* Map gauge to global [0,1]x[0,1] space. This works for the brick
                   but not clear what happens for the cubed sphere */
                double gx,gy,gz;
                gauges_normalize_coordinates(glob,block,nb,g,&gx,&gy,&gz);                

                int gauge_in_block = (xll <= gx && gx < xur) && 
                                     (yll <= gy && gy < yur);
                if (gauge_dim == 3)
                {                    
                    gauge_in_block  &= (zll <= gz && gz < zur);
                }

                if (gauge_in_block)
                {
                    int ng = num_gauges_set;
                    g->blockno = nb;     /* gauge[i] is in block nb */

                    /* Store gauge coordinates in global brick domain.
                       This assumes that coordinates are in a global brick
                       domain in [0,mi]x[0,mj].  This may not work for the 
                       cubed sphere case */

                    double *c = (double*) sc_array_index_int(coordinates, ng);

                    c[0] = mi*(gx - xll);
                    c[1] = mj*(gy - yll);
                    if (gauge_dim == 3)
                    {
                        c[2] = fclaw_opt->mk*(gz - zll);
                    }

                    /* Store location of this gauge in a global list */
                    g->location_in_results = ng;
                    g->in_domain = 1;
                    num_gauges_set++;
                }
            }  /* end of num_gauges */

            /* Set number of gauges for block nb */
            bo = (int*) sc_array_index_int(block_offsets, nb+1);

            bo[0] = num_gauges_set;  
        }
        gauge_acc->num_gauges_set = num_gauges_set;

        /* Each processor checks all gauges and all blocks, so no need for a 
           collective call here. */        

        if (num_gauges_set < num_gauges)            
        {            
            for(int i = 0; i < num_gauges; i++)
                if (!gauges[i].in_domain)
                    fclaw_global_infof("Gauge %d is not in the domain.\n",
                                            gauges[i].num);

            /* Ignore any gauges that are outside of the domain */
            sc_array_resize (coordinates, num_gauges_set);
        }
        *block_offsets_ptr = block_offsets;
        *coordinates_ptr = coordinates;
        return num_gauges_set;
    } 
    else
    {
        *block_offsets_ptr = NULL;
        *coordinates_ptr = NULL;   
        return 0;     
    } /* num_gauges > 0 */
}

/* This is called after regrid and when gauge the positions are updated */
void fclaw_gauges_locate_patches(fclaw_global_t *glob)
{
    fclaw_gauge_acc_t* gauge_acc = 
        (fclaw_gauge_acc_t*) fclaw_diagnostics_get_acc(glob)->gauge_accumulator;

    /* Set block offsets for all gauges;  return number of gauges
       that are actually in the domain */
    fclaw_gauge_info_t* gauge_info = 
        (fclaw_gauge_info_t *) fclaw_global_get_attribute(glob,"gauge_info");

    int num_gauges_set;
    if (!gauge_info->is_initialized || gauge_info->moving_gauges)
    {
        num_gauges_set = fclaw_gauges_setup_block_lists(glob,
                                                        gauge_acc,
                                                        &gauge_info->block_offsets,
                                                        &gauge_info->coordinates); 
        gauge_info->is_initialized = 1;
    }
    else
    {
        /* Gauges are not moving but have been initialized */
        num_gauges_set = gauge_acc->num_gauges_set;
    }

    if (num_gauges_set == 0)
        return;

    /* This calls the main p4est search routine */
    sc_array_t *results = sc_array_new_size(sizeof(int), num_gauges_set);
    fclaw_domain_search_points(glob->domain, 
                               gauge_info->block_offsets,
                               gauge_info->coordinates, results);

    fclaw_gauge_t *gauges = gauge_acc->gauges;
    for (int i = 0; i < num_gauges_set; i++)
    {
        fclaw_gauge_t *g = &gauges[i];

        if (!g->in_domain)
            continue;

        int index = g->location_in_results;
        FCLAW_ASSERT(index >= 0 && index < num_gauges_set);

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
            gauges_print_buffer(glob,g);
            g->next_buffer_location = 0;
        }
    }
    sc_array_destroy(results);

    if (gauge_info->moving_gauges)
    {
        /* We need to relocate gauges within patches each time */
        sc_array_destroy(gauge_info->block_offsets);
        sc_array_destroy(gauge_info->coordinates);
    }
}

/* This will be called after an update on a level */
void  fclaw_gauges_update_positions(fclaw_global_t* glob,
                                    double t, double dt)
{
    fclaw_gauge_acc_t* gauge_acc = 
        (fclaw_gauge_acc_t*) fclaw_diagnostics_get_acc(glob)->gauge_accumulator;

    fclaw_gauge_t *gauges = gauge_acc->gauges;

    const fclaw_gauges_vtable_t* gauge_vt = fclaw_gauges_vt(glob);
    FCLAW_ASSERT(gauge_vt->move != NULL);

    int num_gauges_set = gauge_acc->num_gauges_set;    
    for(int i = 0; i < num_gauges_set; i++)
    {
        fclaw_gauge_t* g = &gauges[i];

        if (!g->in_domain)
            continue;

        if (g->is_moving) 
            gauge_vt->move(glob,g, t, dt);
    }
}


/* Needed for diagnostics */
static
void gauges_compute(fclaw_global_t *glob, void* acc)
{
    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);

    fclaw_gauge_acc_t* gauge_acc = (fclaw_gauge_acc_t*) acc;
    fclaw_gauge_t *gauges = gauge_acc->gauges;

    int buffer_len = fclaw_opt->gauge_buffer_length;
    double tcurr = glob->curr_time;
    //int num_gauges = gauge_acc->num_gauges;
    int num_gauges_set = gauge_acc->num_gauges_set;

    /* Update location;  setup gauges;  determine if they are local or remote */

    /* Only compute those gauges that are actually set */
    for (int i = 0; i < num_gauges_set; i++)
    {
        fclaw_gauge_t *g = &gauges[i];

        /* Ignore any gauge that is not in the domain. */
        if (!g->in_domain)
            continue;

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
                fclaw_block_t *block = &glob->domain->blocks[g->blockno];
                fclaw_patch_t *patch = &block->patches[g->patchno]; 
                gauges_update(glob,block,patch,
                              g->blockno,g->patchno,
                              tcurr,g);

                g->next_buffer_location++;
                
                if (g->next_buffer_location == buffer_len)
                {
                    gauges_print_buffer(glob,g);
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


static
void gauges_finalize(fclaw_global_t *glob, void** acc)
{
    fclaw_gauge_info_t* gauge_info = 
        (fclaw_gauge_info_t *) fclaw_global_get_attribute(glob,"gauge_info");

    /* Clean up gauges and print anything left over in buffers */
    fclaw_gauge_acc_t* gauge_acc = *((fclaw_gauge_acc_t**) acc);


    fclaw_gauge_t *gauges = gauge_acc->gauges;
    int num_gauges = gauge_acc->num_gauges;
    for(int i = 0; i < num_gauges; i++)
    {
        fclaw_gauge_t *g = &gauges[i];

        /* Every processor owns every gauge (which will scale up to a few 
        hundred gauges).  But we only want to print those gauge buffers that 
        for gauges that are on the local processor */        
        if (g->is_local)
            gauges_print_buffer(glob,g);

        FCLAW_FREE(g->buffer);               
    }

    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);
    if (!fclaw_opt->moving_gauges)
    {
        /* info for moving gauges is destroyed at the end of each 
        "locate patches" call.  Only static gauge info needs 
        to be destroyed here. */
        sc_array_destroy(gauge_info->block_offsets);
        sc_array_destroy(gauge_info->coordinates);            
    }

    if (gauge_acc->gauges != NULL)
        FCLAW_FREE(gauge_acc->gauges); 

    /* What about the user gauge ? */
    
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

fclaw_gauges_vtable_t* fclaw_gauges_vt(fclaw_global_t* glob)
{
	fclaw_gauges_vtable_t* gauges_vt = (fclaw_gauges_vtable_t*) 
	   							fclaw_global_get_vtable(glob, "fclaw_gauges");
	FCLAW_ASSERT(gauges_vt != NULL);
	FCLAW_ASSERT(gauges_vt->is_set != 0);

    return gauges_vt;
}

void fclaw_gauges_vtable_initialize(fclaw_global_t* glob)
{
    /* All gauges functions are set in clawpatch routines, since this 
      is where all the information needed to interpolate, etc resides. 
      Also, the data written to gauges is very solver dependent, e.g.
      qvar, aux_var, etc. 
    */
    fclaw_gauges_vtable_t* gauges_vt = fclaw_gauges_vt_new();

    /* This functions will be called if the option `output-gauges` is set 
       to true
    */

    fclaw_diagnostics_vtable_t *diag_vt  = fclaw_diagnostics_vt(glob);
    diag_vt->gauges_init_diagnostics     = gauges_initialize;
    diag_vt->gauges_compute_diagnostics  = gauges_compute;
    diag_vt->gauges_finalize_diagnostics = gauges_finalize;

    gauges_vt->is_set = 1;

	fclaw_global_vtable_store(glob, "fclaw_gauges", gauges_vt, fclaw_gauges_vt_destroy);
}

/* ---------------------------- Access Functions ---------------------------------- */

/* Functions called from outside */

/* 
    This routine identifies the patchnos of all gauges and whether they are local
    or remote.  If a patch is now remote, any current buffer is written out. 

    This routine is called from fclaw_initialize.c and fclaw_regrid.c 
    (if we have a new mesh).  
*/




void fclaw_gauges_allocate(fclaw_global_t *glob, int num_gauges,
                          fclaw_gauge_t **g)
{
    *g = (fclaw_gauge_t*) FCLAW_ALLOC(fclaw_gauge_t,num_gauges);
}

void fclaw_gauges_set_data(fclaw_global_t *glob, 
                             fclaw_gauge_t *g,
                             int num, int dim,
                             double xc, double yc, double zc,
                             double  t1, double t2, 
                             double min_time_increment)
{
    g->num = num;
    g->dim = dim;
    g->xc = xc;
    g->yc = yc;
    g->zc = zc;
    g->t1 = t1;
    g->t2 = t2;
    g->min_time_increment = min_time_increment;
}

void fclaw_gauges_get_data(fclaw_global_t *glob, 
                          fclaw_gauge_t *g,
                          int *num, int *dim,
                          double *xc, double *yc, double *zc,
                          double  *t1, double *t2)
{
    *num = g->num;
    *dim = g->dim;
    *xc = g->xc;
    *yc = g->yc;
    *zc = g->zc;
    *t1 = g->t1;
    *t2 = g->t2;
}

void fclaw_gauges_set_position(fclaw_global_t *glob, 
                               fclaw_gauge_t *g,
                               double xc, double yc, double zc)
{
    g->xc = xc;
    g->yc = yc;
    g->zc = zc;
}

int fclaw_gauges_get_id(fclaw_global_t *glob, 
                          fclaw_gauge_t *g)
{
    return g->num;
}


void fclaw_gauges_get_buffer(fclaw_global_t *glob,
                            fclaw_gauge_t *g,
                            int *kmax, void*** gauge_buffer)
{
    *kmax = g->next_buffer_location;
    *gauge_buffer = g->buffer;
}

void fclaw_gauges_set_buffer_entry(fclaw_global_t *glob,
                                  fclaw_gauge_t* g,
                                  void* guser)
{
    int k = g->next_buffer_location;
    g->buffer[k] = guser;
}

void fclaw_gauges_set_user_data(fclaw_global_t *glob,
                               fclaw_gauge_t* g,
                               void* user)
{
    g->user_data = user;
}

void* fclaw_gauges_get_user_data(fclaw_global_t *glob,
                                  fclaw_gauge_t* g)
{
    return g->user_data;
}


#ifdef __cplusplus
}
#endif
