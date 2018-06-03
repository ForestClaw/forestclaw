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

#ifndef FC2D_GEOCLAW_GAUGES_H
#define FC2D_GEOCLAW_GAUGES_H

#include <fclaw_base.h>  /* Needed for fclaw_app_t */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fc2d_geoclaw_gauges_vtable fc2d_geoclaw_gauges_vtable_t;

struct fclaw2d_patch;

typedef struct fc2d_geoclaw_gauge
{
    int blockno;
    int location_in_results;

    /* Some data needed to get around fact that in parallel, we don't communicate
       gauge information */
    int is_local;
    int patchno;

    double xc;
    double yc;
    double t1;   /* Tstart */
    double t2;   /* Tend */
    int num;     /* Gauge number 1001, 1002, 1003, ...*/

    /* Control output times for gauges */
    double min_time_increment;  /* Output gauges this often */
    double last_time;           /* Last time we output gauge */

    /* Store data in buffer before outputting gauges */
    int next_buffer_location;     /* Where are we in the gauge output */
    void **buffer;
} fc2d_geoclaw_gauge_t;

struct fclaw2d_global;
struct fc2d_geoclaw_gauge;

typedef void (*fc2d_geoclaw_read_gauges_data_t)(struct fclaw2d_global *glob, 
                                                struct fc2d_geoclaw_gauge **gauges, 
                                                int *num);

typedef void (*fc2d_geoclaw_create_gauge_files_t)(struct fclaw2d_global *glob, 
                                                  struct fc2d_geoclaw_gauge *gauges, 
                                                  int num_gauges);

typedef void (*fc2d_geoclaw_gauge_update_t)(struct fclaw2d_global* glob, 
                                            struct fclaw2d_block* block,
                                            struct fclaw2d_patch* patch, 
                                            int blockno, int patchno,
                                            double tcurr, struct fc2d_geoclaw_gauge *g);

typedef void (*fc2d_geoclaw_gauge_print_t)(struct fclaw2d_global *glob, 
                                           struct fc2d_geoclaw_gauge *gauge);

typedef void (*fc2d_geoclaw_gauge_destroy_buffer_data_t)(struct fclaw2d_global *glob, 
                                                         void* gdata);




struct fc2d_geoclaw_gauges_vtable
{
    fc2d_geoclaw_read_gauges_data_t    read_gauges_data;
    fc2d_geoclaw_create_gauge_files_t  create_gauge_files;
    fc2d_geoclaw_gauge_update_t        update_gauge;
    fc2d_geoclaw_gauge_print_t        print_gauge_buffer;

    int is_set;
};

void fc2d_geoclaw_locate_gauges(struct fclaw2d_global *glob);

void fc2d_geoclaw_gauges_vtable_set();

fc2d_geoclaw_gauges_vtable_t* fc2d_geoclaw_gauges_vt();

/* ---------------------------------- Gauges ------------------------------------------ */

void fc2d_geoclaw_gauge_allocate(struct fclaw2d_global *glob, int num_gauges,
                                 struct fc2d_geoclaw_gauge **g);

void fc2d_geoclaw_gauge_set_data(struct fclaw2d_global *glob, 
                                 struct fc2d_geoclaw_gauge *g,
                                 int num, 
                                 double xc, double yc, 
                                 double  t1, double t2,
                                 double min_time_increment);

void fc2d_geoclaw_gauge_get_data(struct fclaw2d_global *glob, 
                                 struct fc2d_geoclaw_gauge *g,                             
                                 int *num, 
                                 double *xc, double *yc, 
                                 double  *t1, double *t2);

void fc2d_geoclaw_gauge_set_buffer_entry(struct fclaw2d_global *glob,
                                           struct fc2d_geoclaw_gauge* g,
                                           void* guser);

void fc2d_geoclaw_gauge_get_buffer(struct fclaw2d_global *glob,
                                   struct fc2d_geoclaw_gauge *g,
                                   int *kmax, void*** gauge_buffer);

void fc2d_geoclaw_gauge_buffer_reset(struct fclaw2d_global *glob, 
                                     struct fc2d_geoclaw_gauge *g);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
