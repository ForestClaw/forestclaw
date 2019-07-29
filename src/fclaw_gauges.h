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

#ifndef FCLAW_GAUGES_H
#define FCLAW_GAUGES_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw_gauges_vtable fclaw_gauges_vtable_t;

typedef struct fclaw_gauge
{
    int blockno;
    int location_in_results;

    /* Some data needed to get around fact that in parallel, we don't communicate
       gauge information */
    int is_local;
    int patchno;

    /* Relative to [ax,ay]x[bx,by] set in fclaw2d_options */
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

    void* user_data;  /* Data about each gauge that doesn't change */

} fclaw_gauge_t;

struct fclaw2d_global;
struct fclaw2d_patch;
struct fclaw2d_block;
struct fclaw_gauge;


typedef void (*fclaw_gauge_set_data_t)(struct fclaw2d_global *glob, 
                                       struct fclaw_gauge **gauges, 
                                       int *num);

typedef void (*fclaw_gauge_create_files_t)(struct fclaw2d_global *glob, 
                                           struct fclaw_gauge *gauges, 
                                           int num_gauges);

typedef void (*fclaw_gauge_normalize_t)(struct fclaw2d_global *glob, 
                                       struct fclaw2d_block *block,
                                       int blockno, 
                                       struct fclaw_gauge *g,
                                       double *xc, double *yc);


typedef void (*fclaw_gauge_update_t)(struct fclaw2d_global* glob, 
                                     struct fclaw2d_block* block,
                                     struct fclaw2d_patch* patch, 
                                     int blockno, int patchno,
                                     double tcurr, struct fclaw_gauge *g);

typedef void (*fclaw_gauge_print_t)(struct fclaw2d_global *glob, 
                                    struct fclaw_gauge *gauge);

typedef void (*fclaw_gauge_destroy_buffer_data_t)(struct fclaw2d_global *glob, 
                                                  void* gdata);


struct fclaw_gauges_vtable
{
    fclaw_gauge_set_data_t      set_gauge_data;
    fclaw_gauge_create_files_t  create_gauge_files;
    fclaw_gauge_update_t        update_gauge;
    fclaw_gauge_print_t         print_gauge_buffer;
    fclaw_gauge_normalize_t     normalize_coordinates;

    int is_set;
};

void fclaw_locate_gauges(struct fclaw2d_global *glob);

void fclaw_gauges_vtable_initialize();

fclaw_gauges_vtable_t* fclaw_gauges_vt();



/* ------------------------ Virtualized gauge functions ------------------------------- */

void fclaw_set_gauge_data(struct fclaw2d_global* glob, 
                          struct fclaw_gauge **gauges, 
                          int *num_gauges);

void fclaw_create_gauge_files(struct fclaw2d_global* glob, 
                              struct fclaw_gauge *gauges, 
                              int num_gauges);

void fclaw_gauge_normalize_coordinates(struct fclaw2d_global *glob, 
                                      struct fclaw2d_block *block,
                                      int blockno, 
                                      struct fclaw_gauge *g,
                                      double *xc, double *yc);


void  fclaw_update_gauge(struct fclaw2d_global* glob, 
                         struct fclaw2d_block *block,
                         struct fclaw2d_patch *patch,
                         int blockno, int patchno,
                         double tcurr, fclaw_gauge_t *g);

void fclaw_print_gauge_buffer(struct fclaw2d_global* glob, 
                              struct fclaw_gauge *g);


/* ---------------------------------- Gauges ------------------------------------------ */

void fclaw_gauge_allocate(struct fclaw2d_global *glob, int num_gauges,
                          struct fclaw_gauge **g);

void fclaw_gauge_set_data(struct fclaw2d_global *glob, 
                          struct fclaw_gauge *g,
                          int num, 
                          double xc, double yc, 
                          double  t1, double t2,
                          double min_time_increment);

void fclaw_gauge_get_data(struct fclaw2d_global *glob, 
                          struct fclaw_gauge *g,                             
                          int *num, 
                          double *xc, double *yc, 
                          double  *t1, double *t2);

int fclaw_gauge_get_id(struct fclaw2d_global *glob, 
                       struct fclaw_gauge *g);
    

void fclaw_gauge_set_buffer_entry(struct fclaw2d_global *glob,
                                  struct fclaw_gauge* g,
                                  void* guser);

void fclaw_gauge_get_buffer(struct fclaw2d_global *glob,
                            struct fclaw_gauge *g,
                            int *kmax, void*** gauge_buffer);

void fclaw_gauge_buffer_reset(struct fclaw2d_global *glob, 
                              struct fclaw_gauge *g);


void fclaw_gauge_set_user_data(struct fclaw2d_global *glob,
                               struct fclaw_gauge* g,
                               void* user);

void* fclaw_gauge_get_user_data(struct fclaw2d_global *glob,
                                struct fclaw_gauge* g);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
