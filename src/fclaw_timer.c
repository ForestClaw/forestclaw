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

#include <fclaw_timer.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>

#include <fclaw2d_partition.h>
#include <sc_statistics.h>
#include <fclaw2d_domain.h>

#define  PRIORITY_WALL          FCLAW_TIMER_PRIORITY_WALL
#define  PRIORITY_EXCLUSIVE1    FCLAW_TIMER_PRIORITY_SUMMARY
#define  PRIORITY_EXCLUSIVE2    FCLAW_TIMER_PRIORITY_EXCLUSIVE
#define  PRIORITY_COUNTERS1     FCLAW_TIMER_PRIORITY_COUNTERS
#define  PRIORITY_COUNTERS2     FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_REGRID        FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_PARTITION     FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_ADVANCE       FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_GHOST         FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_SEARCH        FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_COMM          FCLAW_TIMER_PRIORITY_DETAILS
#define  PRIORITY_CUDA          FCLAW_TIMER_PRIORITY_CUDA
#define  PRIORITY_EXTRA         FCLAW_TIMER_PRIORITY_EXTRA


#define FCLAW2D_STATS_SET(stats,glob,NAME) do {               \
    SC_CHECK_ABORT (!(glob)->timers[FCLAW2D_TIMER_ ## NAME].running,              \
                    "Timer " #NAME " still running in fclaw2d_domain_finalize");  \
    sc_stats_set1 ((stats) + FCLAW2D_TIMER_ ## NAME,                              \
                   (glob)->timers[FCLAW2D_TIMER_ ## NAME].cumulative, #NAME);     \
} while (0)

#define FCLAW2D_STATS_SET_GROUP(stats,NAME,GROUP) do {                  \
    sc_stats_set_group_prio ((stats) + FCLAW2D_TIMER_ ## NAME,          \
                   GROUP_ ## GROUP, PRIORITY_ ## GROUP);    \
} while (0)



/* -----------------------------------------------------------------
   Work with timers
   ----------------------------------------------------------------- */

double
fclaw2d_timer_wtime (void)
{
    return sc_MPI_Wtime ();
}

void
fclaw2d_timer_init (fclaw2d_timer_t *timer)
{
    memset (timer, 0, sizeof (fclaw2d_timer_t));
}

void
fclaw2d_timer_start (fclaw2d_timer_t *timer)
{
    if (!timer->running) 
    {
        timer->started = fclaw2d_timer_wtime ();
        timer->stopped = 0.;
        timer->running = 1;
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void
fclaw2d_timer_stop (fclaw2d_timer_t *timer)
{
    if (timer->running) 
    {
        timer->stopped = fclaw2d_timer_wtime ();
        timer->cumulative += timer->stopped - timer->started;
        timer->running = 0;
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

/* Don't put timers inside of any step update functions when using OPENMP */
void fclaw2d_timer_start_threadsafe(fclaw2d_timer_t *timer)
{
#if !defined(_OPENMP)
    fclaw2d_timer_start(timer);
#endif    
}

void fclaw2d_timer_stop_threadsafe(fclaw2d_timer_t *timer)
{
#if !defined(_OPENMP)
    fclaw2d_timer_stop(timer);
#endif    
}

void
fclaw2d_timer_report(fclaw2d_global_t *glob)
{
    sc_statinfo_t stats[FCLAW2D_TIMER_COUNT];


    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_WALLTIME]);

    FCLAW2D_STATS_SET (stats, glob, INIT);
    FCLAW2D_STATS_SET (stats, glob, OUTPUT);
    FCLAW2D_STATS_SET (stats, glob, DIAGNOSTICS);
    FCLAW2D_STATS_SET (stats, glob, REGRID);
    FCLAW2D_STATS_SET (stats, glob, ADVANCE);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL);
    FCLAW2D_STATS_SET (stats, glob, ADAPT_COMM);
    FCLAW2D_STATS_SET (stats, glob, PARTITION_COMM);
    FCLAW2D_STATS_SET (stats, glob, GHOSTPATCH_COMM);
    FCLAW2D_STATS_SET (stats, glob, DIAGNOSTICS_COMM);
    FCLAW2D_STATS_SET (stats, glob, CFL_COMM);
    FCLAW2D_STATS_SET (stats, glob, WALLTIME);
    FCLAW2D_STATS_SET (stats, glob, REGRID_BUILD);
    FCLAW2D_STATS_SET (stats, glob, REGRID_TAGGING);
    FCLAW2D_STATS_SET (stats, glob, TIMESYNC);
    FCLAW2D_STATS_SET (stats, glob, PARTITION);
    FCLAW2D_STATS_SET (stats, glob, PARTITION_BUILD);
    FCLAW2D_STATS_SET (stats, glob, ADVANCE_STEP2);
    FCLAW2D_STATS_SET (stats, glob, ADVANCE_B4STEP2);
    FCLAW2D_STATS_SET (stats, glob, GHOSTPATCH_BUILD);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_COPY);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_AVERAGE);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_INTERP);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_PHYSBC);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_STEP1);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_STEP2);
    FCLAW2D_STATS_SET (stats, glob, GHOSTFILL_STEP3);
    FCLAW2D_STATS_SET (stats, glob, NEIGHBOR_SEARCH);
    FCLAW2D_STATS_SET (stats, glob, CUDA_ALLOCATE);
    FCLAW2D_STATS_SET (stats, glob, CUDA_MEMCOPY_H2H);
    FCLAW2D_STATS_SET (stats, glob, CUDA_MEMCOPY_H2D);
    FCLAW2D_STATS_SET (stats, glob, CUDA_MEMCOPY_D2H);
    FCLAW2D_STATS_SET (stats, glob, EXTRA1);
    FCLAW2D_STATS_SET (stats, glob, EXTRA2);
    FCLAW2D_STATS_SET (stats, glob, EXTRA3);
    FCLAW2D_STATS_SET (stats, glob, EXTRA4);

    int d = glob->count_grids_per_proc;
    glob->count_grids_per_proc = (d > 0) ? d : 1;   /* To avoid division by zero */
    d = glob->count_amr_advance;
    glob->count_amr_advance = (d > 0) ? d : 1;   /* To avoid division by zero */

    double gpp = glob->count_grids_per_proc/       glob->count_amr_advance;
    double glb = glob->count_grids_local_boundary/ glob->count_amr_advance;
    double grb = glob->count_grids_remote_boundary/glob->count_amr_advance;
    double gint = gpp - glb;

    /* compute arithmetic mean of total advance steps per processor */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_ADVANCE_STEPS_COUNTER],
                   glob->count_single_step,"ADVANCE_STEPS_COUNTER");

    /* Compute the arithmetic mean of grids per processor */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_GRIDS_PER_PROC],gpp,"GRIDS_PER_PROC");

    /* Compute the arithmetic mean of grids in the interior */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_GRIDS_INTERIOR],gint,
                   "GRIDS_INTERIOR");

    /* Compute the arithmetic mean of local grids on the boundary */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_GRIDS_LOCAL_BOUNDARY],glb,
                   "GRIDS_LOCAL_BOUNDARY");

    /* Compute the arithmetic mean of remote grids on the boundary */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_GRIDS_REMOTE_BOUNDARY],grb,
                   "GRIDS_REMOTE_BOUNDARY");

    sc_stats_set1 (&stats[FCLAW2D_TIMER_UNACCOUNTED],
                   glob->timers[FCLAW2D_TIMER_WALLTIME].cumulative -
                   (glob->timers[FCLAW2D_TIMER_INIT].cumulative +
                    glob->timers[FCLAW2D_TIMER_REGRID].cumulative +
                    glob->timers[FCLAW2D_TIMER_OUTPUT].cumulative +
                    glob->timers[FCLAW2D_TIMER_DIAGNOSTICS].cumulative +
                    glob->timers[FCLAW2D_TIMER_ADVANCE].cumulative +
                    glob->timers[FCLAW2D_TIMER_GHOSTFILL].cumulative +
                    glob->timers[FCLAW2D_TIMER_ADAPT_COMM].cumulative +
                    glob->timers[FCLAW2D_TIMER_GHOSTPATCH_COMM].cumulative +
                    glob->timers[FCLAW2D_TIMER_PARTITION_COMM].cumulative +
                    glob->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM].cumulative +
                    glob->timers[FCLAW2D_TIMER_CFL_COMM].cumulative),
                   "UNACCOUNTED");

    sc_stats_set1 (&stats[FCLAW2D_TIMER_GLOBAL_COMM],
                   glob->timers[FCLAW2D_TIMER_ADAPT_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_GHOSTPATCH_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_PARTITION_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_CFL_COMM].cumulative,
                   "FCLAW2D_TIMER_GLOBAL_COMM");

    /* Just subtracting FCLAW2D_TIMER_GLOBAL here doesn't work ... */
    sc_stats_set1 (&stats[FCLAW2D_TIMER_LOCAL_COMM],
                   glob->timers[FCLAW2D_TIMER_WALLTIME].cumulative -
                   (glob->timers[FCLAW2D_TIMER_ADAPT_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_GHOSTPATCH_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_PARTITION_COMM].cumulative +
                   glob->timers[FCLAW2D_TIMER_DIAGNOSTICS_COMM].cumulative +
                    glob->timers[FCLAW2D_TIMER_CFL_COMM].cumulative),
                   "FCLAW2D_TIMER_LOCAL_COMM");


    /* --------------------------------- Set stats groups ------------------------------*/

    /* Names for timer groups */
    enum 
    {
        GROUP_NONE = -1,
        GROUP_EXCLUSIVE1,
        GROUP_EXCLUSIVE2,
        GROUP_WALL,
        GROUP_COUNTERS1,
        GROUP_COUNTERS2,        
        GROUP_COMM,
        GROUP_REGRID,
        GROUP_PARTITION,
        GROUP_ADVANCE,
        GROUP_GHOST,
        GROUP_SEARCH,
        GROUP_CUDA,
        GROUP_EXTRA,
        GROUP_COUNT
    };

    FCLAW2D_STATS_SET_GROUP(stats,WALLTIME,              WALL);

    FCLAW2D_STATS_SET_GROUP(stats,ADVANCE,               EXCLUSIVE1);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL,             EXCLUSIVE1);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTPATCH_COMM,       EXCLUSIVE1);
    FCLAW2D_STATS_SET_GROUP(stats,REGRID,                EXCLUSIVE1);
    FCLAW2D_STATS_SET_GROUP(stats,ADAPT_COMM,            EXCLUSIVE1);

    FCLAW2D_STATS_SET_GROUP(stats,INIT,                  EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,OUTPUT,                EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,DIAGNOSTICS,           EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,PARTITION_COMM,        EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,DIAGNOSTICS_COMM,      EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,CFL_COMM,              EXCLUSIVE2);
    FCLAW2D_STATS_SET_GROUP(stats,UNACCOUNTED,           EXCLUSIVE2);

    FCLAW2D_STATS_SET_GROUP(stats,LOCAL_COMM,            COMM);
    FCLAW2D_STATS_SET_GROUP(stats,GLOBAL_COMM,           COMM);

    FCLAW2D_STATS_SET_GROUP(stats,ADVANCE_STEPS_COUNTER, COUNTERS1);
    FCLAW2D_STATS_SET_GROUP(stats,GRIDS_PER_PROC,        COUNTERS1);

    FCLAW2D_STATS_SET_GROUP(stats,GRIDS_INTERIOR,        COUNTERS2);
    FCLAW2D_STATS_SET_GROUP(stats,GRIDS_LOCAL_BOUNDARY,  COUNTERS2);
    FCLAW2D_STATS_SET_GROUP(stats,GRIDS_REMOTE_BOUNDARY, COUNTERS2);

    FCLAW2D_STATS_SET_GROUP(stats,REGRID_BUILD,          REGRID);
    FCLAW2D_STATS_SET_GROUP(stats,REGRID_TAGGING,        REGRID);

    FCLAW2D_STATS_SET_GROUP(stats,PARTITION,             PARTITION);  
    FCLAW2D_STATS_SET_GROUP(stats,PARTITION_BUILD,       PARTITION);

    FCLAW2D_STATS_SET_GROUP(stats,ADVANCE_STEP2,         ADVANCE);
    FCLAW2D_STATS_SET_GROUP(stats,ADVANCE_B4STEP2,       ADVANCE);

    FCLAW2D_STATS_SET_GROUP(stats,GHOSTPATCH_BUILD,      GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,TIMESYNC,              GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_COPY,        GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_AVERAGE,     GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_INTERP,      GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_PHYSBC,      GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_STEP1,       GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_STEP2,       GHOST);
    FCLAW2D_STATS_SET_GROUP(stats,GHOSTFILL_STEP3,       GHOST);

    FCLAW2D_STATS_SET_GROUP(stats,NEIGHBOR_SEARCH,       SEARCH);

    FCLAW2D_STATS_SET_GROUP(stats,CUDA_ALLOCATE,         CUDA);
    FCLAW2D_STATS_SET_GROUP(stats,CUDA_MEMCOPY_H2H,      CUDA);
    FCLAW2D_STATS_SET_GROUP(stats,CUDA_MEMCOPY_H2D,      CUDA);
    FCLAW2D_STATS_SET_GROUP(stats,CUDA_MEMCOPY_D2H,      CUDA);


    FCLAW2D_STATS_SET_GROUP(stats,EXTRA1,                EXTRA);
    FCLAW2D_STATS_SET_GROUP(stats,EXTRA2,                EXTRA);
    FCLAW2D_STATS_SET_GROUP(stats,EXTRA3,                EXTRA);
    FCLAW2D_STATS_SET_GROUP(stats,EXTRA4,                EXTRA);


    /* Get a partial sum of timers not accounted for in reported summary */
    int priority = fclaw_opt->report_timing_verbosity;  /* Make this an option */

    /* ----------------------------------- Compute timers ------------------------------*/

    /* This does all-reduce, etc to set stats */
    sc_stats_compute (glob->mpicomm, FCLAW2D_TIMER_COUNT, stats);

    /* ------------------------------------ Print timers ------------------------------*/


    sc_stats_print_ext(sc_package_id, SC_LP_ESSENTIAL, FCLAW2D_TIMER_COUNT,
                       stats,sc_stats_group_all,priority,1,0);


    SC_GLOBAL_ESSENTIALF ("Procs %d advance %d %g exchange %d %g "
                          "regrid %d %d %g\n", glob->mpisize,
                          glob->count_amr_advance,
                          stats[FCLAW2D_TIMER_ADVANCE].average,
                          glob->count_ghost_exchange,
                          stats[FCLAW2D_TIMER_GHOSTFILL].average,
                          glob->count_amr_regrid,
                          glob->count_amr_new_domain,
                          stats[FCLAW2D_TIMER_REGRID].average);

    SC_GLOBAL_ESSENTIALF ("Max/P %d advance %d %g exchange %d %g "
                          "regrid %d %d %g\n", glob->mpisize,
                          glob->count_amr_advance,
                          stats[FCLAW2D_TIMER_ADVANCE].max,
                          glob->count_ghost_exchange,
                          stats[FCLAW2D_TIMER_GHOSTFILL].max,
                          glob->count_amr_regrid,
                          glob->count_amr_new_domain,
                          stats[FCLAW2D_TIMER_REGRID].max);

#if 0
    /* Find out process rank */
    /* TODO : Fix this so that it doesn't interfere with output printed above. */
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    /* Write out individual processor timers */
    printf("%12s time on proc %d : %12.4f\n","ADVANCE",
           domain->mpirank,glob->timers[FCLAW2D_TIMER_ADVANCE].cumulative);
    printf("%12s time on proc %d : %12.4f\n","GHOSTCOMM",
           domain->mpirank,glob->timers[FCLAW2D_TIMER_GHOSTCOMM].cumulative);
    printf("%12s time on proc %d : %12.4f\n","EXCHANGE",
           domain->mpirank,glob->timers[FCLAW2D_TIMER_EXCHANGE].cumulative);
    printf("%12s time on proc %d : %12.4f\n","REGRID",
           domain->mpirank,glob->timers[FCLAW2D_TIMER_REGRID].cumulative);
    printf("\n");
#endif


}
