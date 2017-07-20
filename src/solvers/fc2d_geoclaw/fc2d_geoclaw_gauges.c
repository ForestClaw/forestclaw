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

#include "fc2d_geoclaw_options.h"

#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_fort.h"

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

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

typedef struct geoclaw_gauge
{
    int blockno;
    int patchno;
    int location_in_results;

    double xc;
    double yc;
    double t1;
    double t2;
    int num;
    /* double* buffer; */  /* Not yet used */

} geoclaw_gauge_t;

typedef struct fc2d_geoclaw_gauge_acc
{
    int num_gauges;
    int is_latest_domain;
    struct geoclaw_gauge *gauges;
} fc2d_geoclaw_gauge_acc_t;


typedef struct fc2d_geoclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fc2d_geoclaw_gauge_info_t;


static fc2d_geoclaw_gauge_info_t gauge_info;

static
void geoclaw_gauge_initialize(fclaw2d_global_t* glob, void** acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc;
    gauge_acc = FCLAW_ALLOC(fc2d_geoclaw_gauge_acc_t,1);
    *acc = gauge_acc;

    const fclaw_options_t * gparms = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    /* --------------------------------------------------------
       Read gauges files 'gauges.data' to get number of gauges
       -------------------------------------------------------- */
    char fname[] = "gauges.data";
    int num = FC2D_GEOCLAW_GAUGES_GETNUM(fname);
    int restart = 0;

    gauge_acc->num_gauges = num;
    gauge_acc->gauges = FCLAW_ALLOC(geoclaw_gauge_t,num);
    
    if (num == 0)
    {
        return;
    }

    /* Read gauges file for the locations, etc. of all gauges */
    FC2D_GEOCLAW_GAUGES_INIT(&restart, &clawpatch_opt->meqn, &num,  gauge_acc->gauges, fname);

    /* -----------------------------------------------------
       Open gauge files and add header information
       ----------------------------------------------------- */
    char filename[15];    /* gaugexxxxx.txt */
    FILE *fp;

    for (int i = 0; i < gauge_acc->num_gauges; ++i)
    {
        geoclaw_gauge_t g = gauge_acc->gauges[i];
        sprintf(filename,"gauge%05d.txt",g.num);
        fp = fopen(filename, "w");
        fprintf(fp, "# gauge_id= %5d location=( %15.7e %15.7e ) num_eqn= %2d\n",
                g.num, g.xc, g.yc, clawpatch_opt->meqn+1);
        fprintf(fp, "# Columns: level time h    hu    hv    eta\n");
        fclose(fp);
    }

    /* -----------------------------------------------------
       Set up block offsets and coordinate list for p4est
       search function
       ----------------------------------------------------- */

    fclaw2d_map_context_t* cont = glob->cont;

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    gauge_info.block_offsets = sc_array_new_size(sizeof(int), glob->domain->num_blocks+1);
    gauge_info.coordinates = sc_array_new_size(2*sizeof(double), num);

    int *block_offsets = (int*) sc_array_index_int(gauge_info.block_offsets, 0);
    double *coordinates = (double*) sc_array_index_int(gauge_info.coordinates, 0);

    geoclaw_gauge_t *gauges = gauge_acc->gauges;

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
            for(int i = 0; i < num; i++)
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
        for(int i = 0; i < num; i++)
        {
            gauges[i].blockno = 0;
            gauges[i].location_in_results = i;
        }

        block_offsets[0] = 0;
        block_offsets[1] = num;

        for (int i = 0; i < num; ++i)
        {
            coordinates[2*i] = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
            coordinates[2*i+1] = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
        }
    }
}



static
void geoclaw_gauge_update(fclaw2d_global_t *glob, void* solver_acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc = (fc2d_geoclaw_gauge_acc_t*) solver_acc;

    const double tcurr = glob->curr_time;
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower,eta;
    double *q, *aux;
    double *var;
    char filename[15];  /* gaugexxxxx.txt */
    FILE *fp;

    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    if (gauge_acc->num_gauges == 0)
    {
        return;
    }

    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    var = FCLAW_ALLOC(double, clawpatch_opt->meqn);
    geoclaw_gauge_t *gauges = gauge_acc->gauges;
    geoclaw_gauge_t g;
    for (int i = 0; i < gauge_acc->num_gauges; ++i)
    {
        g = gauges[i];
        block = &glob->domain->blocks[g.blockno];
        if (g.patchno >= 0)
        {
            patch = &block->patches[g.patchno];
            FCLAW_ASSERT(patch != NULL);
            fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                        &xlower,&ylower,&dx,&dy);

            fclaw2d_clawpatch_soln_data(glob,patch,&q,&meqn);
            fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);

            FCLAW_ASSERT(g.xc >= xlower && g.xc <= xlower+mx*dx);
            FCLAW_ASSERT(g.yc >= ylower && g.yc <= ylower+my*dy);
            if (tcurr >= g.t1 && tcurr <= g.t2)
            {
                FC2D_GEOCLAW_UPDATE_GAUGE(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                      q,&maux,aux,&g.xc,&g.yc,var,&eta);
                sprintf(filename,"gauge%05d.txt",g.num);
                fp = fopen(filename, "a");
                fprintf(fp, "%5d %15.7e %15.7e %15.7e %15.7e %15.7e\n",
                        patch->level,tcurr,var[0],var[1],var[2],eta);
                fclose(fp);
            }
        }
    }
    FCLAW_FREE(var);
}


void fc2d_geoclaw_locate_gauges(fclaw2d_global_t *glob)
{
    int i,index;
    fc2d_geoclaw_gauge_acc_t* gauge_acc = (fc2d_geoclaw_gauge_acc_t*) glob->acc->solver_accumulator;

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
void geoclaw_gauge_finalize(fclaw2d_global_t *glob, void** acc)
{
    fc2d_geoclaw_gauge_acc_t* gauge_acc = *((fc2d_geoclaw_gauge_acc_t**) acc);
    FCLAW_FREE(gauge_acc->gauges);
    FCLAW_FREE(gauge_acc);
    *acc = NULL;

    if (gauge_info.block_offsets != NULL)
    {
        sc_array_destroy(gauge_info.block_offsets);
    }
    if (gauge_info.coordinates != NULL)
    {
        sc_array_destroy(gauge_info.coordinates);
    }    
}


/* --------------------------- Virtual table entries ---------------------------- */

void fc2d_geoclaw_gauges_vtable_set()
{

    fclaw2d_diagnostics_vtable_t *    diag_vt = fclaw2d_diagnostics_vt();

    diag_vt->solver_init_diagnostics     = geoclaw_gauge_initialize;
    diag_vt->solver_compute_diagnostics  = geoclaw_gauge_update;
    diag_vt->solver_finalize_diagnostics = geoclaw_gauge_finalize;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
