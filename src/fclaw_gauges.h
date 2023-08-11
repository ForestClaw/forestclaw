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

#ifndef FCLAW_GAUGES_H
#define FCLAW_GAUGES_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/**
 * @file
 * @brief Gauges
 */


/**
 * @brief Gauge structure
 */
typedef struct fclaw_gauge
{
    /** @brief The block that the gauge is in */
    int blockno;
    /** @brief Location of the gauge in the results array */
    int location_in_results;

    /* Some data needed to get around fact that in parallel, we don't communicate
       gauge information */
    /** @brief true if this gauge is on the local processor */
    int is_local;
    /** @brief the patch that the gauge is in */
    int patchno;

    /** @{ @brief Relative to [ax,ay]x[bx,by] set in fclaw2d_options */
    double xc;   
    double yc;
    /** @} */

    /** @brief Tstart */
    double t1;
    /** @brief Tend */
    double t2;
    /** @brief Gauge number */
    int num;

    /* Control output times for gauges */
    /** @brief Output gauges this often */
    double min_time_increment;
    /** @brief Last time we output gauge */
    double last_time;

    /* Store data in buffer before outputting gauges */
    /** @brief Where are we in the gauge output buffer */
    int next_buffer_location;
    /** @brief Buffer for storing gauge data */
    void **buffer;

    /** @brief User data */
    void* user_data;

} fclaw_gauge_t;

struct fclaw2d_global;
struct fclaw_patch;
struct fclaw2d_block;
struct fclaw_gauge;

/**
 * @brief Sets the data for each gauge
 * 
 * @param[in] glob the global context
 * @param[in,out] gauges the array of gauges
 * @param[in] num_gauges the number of gauges
 */
typedef void (*fclaw_gauge_set_data_t)(struct fclaw2d_global *glob, 
                                       struct fclaw_gauge **gauges, 
                                       int *num);

/**
 * @brief Creates files for each gauge
 * 
 * @param[in] glob the global context
 * @param[in] gauges the array of gauges
 * @param[in] num_gauges the number of gauges
 */
typedef void (*fclaw_gauge_create_files_t)(struct fclaw2d_global *glob, 
                                           struct fclaw_gauge *gauges, 
                                           int num_gauges);

/**
 * @brief Maps gauge to normalized coordinates in a global [0,1]x[0,1]  domain.
 * 
 * @param[in] glob the global context
 * @param[in] block the block that the gauge is in
 * @param[in] blockno the block number that the gauge is in
 * @param[in] g the gauge
 * @param[out] xc,yc the normalized coordinates
 */
typedef void (*fclaw_gauge_normalize_t)(struct fclaw2d_global *glob, 
                                       struct fclaw2d_block *block,
                                       int blockno, 
                                       struct fclaw_gauge *g,
                                       double *xc, double *yc);

/**
 * @brief Updates the current buffer entry for the gauge
 * 
 * @param[in] glob the global context
 * @param[in] block the block that the gauge is in
 * @param[in] patch the patch that the gauge is in
 * @param[in] blockno the block number that the gauge is in
 * @param[in] patchno the patch number that the gauge is in
 * @param[in] tcurr the current time
 * @param[in,out] g the gauge
 */
typedef void (*fclaw_gauge_update_t)(struct fclaw2d_global* glob, 
                                     struct fclaw2d_block* block,
                                     struct fclaw_patch* patch, 
                                     int blockno, int patchno,
                                     double tcurr, struct fclaw_gauge *g);

/**
 * @brief Prints the buffer to a file
 * 
 * @param glob the global context
 * @param g the gauge
 */
typedef void (*fclaw_gauge_print_t)(struct fclaw2d_global *glob, 
                                    struct fclaw_gauge *gauge);

/**
 * @brief vtable for gauges
 */
typedef struct fclaw_gauges_vtable
{
    /** @brief Sets the data for each gauge */
    fclaw_gauge_set_data_t      set_gauge_data;
    /** @brief Creates files for each gauge */
    fclaw_gauge_create_files_t  create_gauge_files;
    /** @brief Updates the current buffer entry for the gauge */
    fclaw_gauge_update_t        update_gauge;
    /** @brief Prints the buffer to a file */
    fclaw_gauge_print_t         print_gauge_buffer;

    /** @brief Maps gauge to normalized coordinates in a global [0,1]x[0,1]  domain. */
    fclaw_gauge_normalize_t     normalize_coordinates;

    /** @brief true if vtable has been set */
    int is_set;
} fclaw_gauges_vtable_t;

/**
 * @brief Locate the gauges in the mesh
 * 
 * @param glob the global context
 */
void fclaw_locate_gauges(struct fclaw2d_global *glob);

/**
 * @brief Initialize the gauges vtable
 * 
 * @param glob the global context
 */
void fclaw_gauges_vtable_initialize(struct fclaw2d_global *glob);

/**
 * @brief Get the gauges vtable
 * 
 * @param glob the global context
 * @return fclaw_gauges_vtable_t* the vtable
 */
fclaw_gauges_vtable_t* fclaw_gauges_vt(struct fclaw2d_global *glob);



/* ------------------------ Virtualized gauge functions ------------------------------- */

/**
 * @brief Set the data for each gauge
 * 
 * @param glob the global context
 * @param gauges the array of gauges
 * @param num_gauges the number of gauges
 */
void fclaw_set_gauge_data(struct fclaw2d_global* glob, 
                          struct fclaw_gauge **gauges, 
                          int *num_gauges);

/**
 * @brief Create files for each gauge
 * 
 * @param[in] glob the global context
 * @param[in] gauges the array of gauges
 * @param[in] num_gauges the number of gauges
 */
void fclaw_create_gauge_files(struct fclaw2d_global* glob, 
                              struct fclaw_gauge *gauges, 
                              int num_gauges);

/**
 * @brief Map gauge to normalized coordinates in a global [0,1]x[0,1]  domain.
 * 
 * @param[in] glob the global context
 * @param[in] block the block that the gauge is in
 * @param[in] blockno the block number that the gauge is in
 * @param[in] g the gauge
 * @param[out] xc,yc the normalized coordinates
 */
void fclaw_gauge_normalize_coordinates(struct fclaw2d_global *glob, 
                                      struct fclaw2d_block *block,
                                      int blockno, 
                                      struct fclaw_gauge *g,
                                      double *xc, double *yc);

/**
 * @brief Update the current buffer entry for the gauge
 * 
 * @param[in] glob the global context
 * @param[in] block the block that the gauge is in
 * @param[in] patch the patch that the gauge is in
 * @param[in] blockno the block number that the gauge is in
 * @param[in] patchno the patch number that the gauge is in
 * @param[in] tcurr the current time
 * @param[in,out] g the gauge
 */
void  fclaw_update_gauge(struct fclaw2d_global* glob, 
                         struct fclaw2d_block *block,
                         struct fclaw_patch *patch,
                         int blockno, int patchno,
                         double tcurr, fclaw_gauge_t *g);

/**
 * @brief Print the buffer to a file
 * 
 * @param glob the global context
 * @param g the gauge
 */
void fclaw_print_gauge_buffer(struct fclaw2d_global* glob, 
                              struct fclaw_gauge *g);


/* ---------------------------------- Gauges ------------------------------------------ */

/**
 * @brief Allocate gauge structures
 * 
 * @param[in] glob the global context
 * @param[in] num_gauges the number of gauges
 * @param[out] g allocated array of gauges 
 */
void fclaw_gauge_allocate(struct fclaw2d_global *glob, int num_gauges,
                          struct fclaw_gauge **g);

/**
 * @brief Set the data for a single gauge
 * 
 * @param[in] glob the global context
 * @param[in,out] g the gauge
 * @param[in] num the gauge number (index the the guage array)
 * @param[in] xc, yc Relative to [ax,ay]x[bx,by] set in fclaw2d_options 
 * @param[in] t1 Tstart
 * @param[in] t2 Tend
 * @param[in] min_time_increment How often to output the gauge
 */
void fclaw_gauge_set_data(struct fclaw2d_global *glob, 
                          struct fclaw_gauge *g,
                          int num, 
                          double xc, double yc, 
                          double  t1, double t2,
                          double min_time_increment);

/**
 * @brief Get the data for a single gauge
 * 
 * @param[in] glob the global context
 * @param[in] g the gauge
 * @param[out] num the gauge number (index the the guage array)
 * @param[out] xc, yc Relative to [ax,ay]x[bx,by] set in fclaw2d_options 
 * @param[out] t1 Tstart
 * @param[out] t2 Tend
 */
void fclaw_gauge_get_data(struct fclaw2d_global *glob, 
                          struct fclaw_gauge *g,                             
                          int *num, 
                          double *xc, double *yc, 
                          double  *t1, double *t2);

/**
 * @brief Get the gauge number (index the the guage array)
 * 
 * @param glob the global context
 * @param g the gauge
 * @return int the gauge number
 */
int fclaw_gauge_get_id(struct fclaw2d_global *glob, 
                       struct fclaw_gauge *g);
    
/**
 * @brief Set data for current buffer entry
 * 
 * @param glob the global context
 * @param g the gauge
 * @param guser the data 
 */
void fclaw_gauge_set_buffer_entry(struct fclaw2d_global *glob,
                                  struct fclaw_gauge* g,
                                  void* guser);

/**
 * @brief get the buffer data
 * 
 * @param[in] glob the global context
 * @param[in] g the gauge
 * @param[out] kmax size of buffer
 * @param[out] gauge_buffer the buffer
 */
void fclaw_gauge_get_buffer(struct fclaw2d_global *glob,
                            struct fclaw_gauge *g,
                            int *kmax, void*** gauge_buffer);

/**
 * @brief the the user data for a gauge
 * 
 * @param[in] glob the global context
 * @param[in] g the gauge
 * @param[in] user the user data
 */
void fclaw_gauge_set_user_data(struct fclaw2d_global *glob,
                               struct fclaw_gauge* g,
                               void* user);

/**
 * @brief Get the user data for gauge
 * 
 * @param glob the global context
 * @param g the gauge
 * @return void* the user data
 */
void* fclaw_gauge_get_user_data(struct fclaw2d_global *glob,
                                struct fclaw_gauge* g);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
