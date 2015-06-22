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

/*
 * IDEA 1: add boolean about on_parallel_boundary vs. interior.
 *         1. Do the parallel boundary work here.
 *            To the least possible amount of work before sending.
 *             * copy_ghost_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *               (which subset of cells shall do it?)
 *               (maybe we can skip doing it before sending altogether?)
 *            More specifically, only copy/average across faces
 *            between two parallel boundary patches, not corners.
 *            Argue in paper why no corners needed for 5-point stencils.
 *         2. Start sending: ghost_exchange_begin.
 *         3. Do all the work on interior patches and whatever has not
 *            been computed between parallel boundary patches.
 *            Do the most possible amount of work while communicating.
 *             * copy_ghost_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *               (how much of this can we possibly to here: maximize this.)
 *             * interpolate_ghost_c2f for interior patches and
 *               faces between parallel boundary and interior patches.
 *         4. Recieve messages: ghost_exchange_end.
 *         5. Work on receive buffers / parallel patches with remote data.
 *             * copy_ghost_samelevel
 *             * average_ghost_f2c
 *             * fill_physical_ghost
 *             * interpolate_ghost_c2f
 *         6. Repeat 5. for whatever parts of patches that are not done yet.
 *            Integrate this with 5. so we don't loop over patches twice.
 *
 * All of this goes into the paper as detailed algorithm with explanations.
 */


/** \file
 * Fill ghost cells.
 *
 *
 **/

#include <fclaw2d_forestclaw.h>
#include <fclaw_timer.h>

#include <fclaw2d_ghost_fill.h>
#include <fclaw2d_partition.h>
#include <fclaw2d_exchange.h>
#include <fclaw2d_domain.h>


#ifdef __cplusplus
extern "C" {
#if 0
}
#endif
#endif

/* -------------------------------------------------
   Basic routines - operate on a single level
   ----------------------------------------------- */
static
void copy2ghost(fclaw2d_domain_t *domain,
                int level,
                fclaw_bool time_interp,
                fclaw_bool read_parallel_patches);

static
void average2ghost(fclaw2d_domain_t *domain,
                   int coarse_level,
                   fclaw_bool time_interp,
                   fclaw_bool read_parallel_patches);
static
void interpolate2ghost(fclaw2d_domain_t *domain,
                       int fine_level,
                       fclaw_bool time_interp,
                       fclaw_bool read_parallel_patches);



/* -------------------------------------------------
   Loop over all levels
   ----------------------------------------------- */
static
void copy_ghost_samelevel(fclaw2d_domain_t* domain,
                          int minlevel,
                          int maxlevel,
                          fclaw_bool time_interp,
                          fclaw_bool read_parallel_patches)
{
    /* Copy between grids that are at the same level. */
    for(int level = maxlevel; level >= minlevel; level--)
    {
        int time_interp = 0;
        copy2ghost(domain,level,time_interp,read_parallel_patches);
    }
    if (time_interp)
    {
        int time_interp_level = minlevel-1;
        copy2ghost(domain,time_interp_level,time_interp,read_parallel_patches);
    }

}


/**
 * \ingroup Averaging
 * Fill in coarse grid ghost cells by averaging or copying  from neighboring fine grids.
 **/
static
void average_fine2coarse_ghost(fclaw2d_domain_t *domain,
                               int mincoarse,
                               int maxcoarse,
                               fclaw_bool time_interp,
                               fclaw_bool read_parallel_patches)
{
    /* Average fine grids to coarse grid ghost cells */
    for(int level = maxcoarse; level >= mincoarse; level--)
    {
        /* Don't do any time interpolation yet */
        int time_interp = 0;
        average2ghost(domain,level,time_interp,read_parallel_patches);
    }
    if (time_interp)
    {
        /* Average fine grid to coarse time interpolated level.  Time interpolated
         faces need correct values.  */
        int time_interp_level = mincoarse-1;
        average2ghost(domain,time_interp_level,time_interp,read_parallel_patches);
    }
}

static
void interpolate_coarse2fine_ghost(fclaw2d_domain_t* domain,
                                   int minfine,
                                   int maxfine,
                                   fclaw_bool time_interp,
                                   fclaw_bool read_parallel_patches)
{
    /* Interpolate from coarse grid to fine grid ghost */
    for(int level = maxfine; level >= minfine; level--)
    {
        /* No need to interpolate to coarse time-interpolated grids */
        int time_interp = 0;
        interpolate2ghost(domain,level,time_interp,read_parallel_patches);
    }
    if (time_interp)
    {
        /* interpolate from time interpolated level */
        int mincoarse = minfine-1;
        interpolate2ghost(domain,mincoarse,time_interp,read_parallel_patches);
    }
}

static
void copy2ghost(fclaw2d_domain_t *domain,
                int level,
                fclaw_bool time_interp,
                fclaw_bool read_parallel_patches)
{
    fclaw2d_exchange_info_t e_info;
    e_info.exchange_type = FCLAW2D_COPY;
    e_info.grid_type = FCLAW2D_IS_COARSE;
    e_info.time_interp = time_interp;
    e_info.read_parallel_patches = read_parallel_patches;

    /* face exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_face_fill,
                                 (void *) &e_info);

    /* corner exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_corner_fill,
                                 (void *) &e_info);
}


static
void average2ghost(fclaw2d_domain_t *domain,
                   int coarse_level,
                   fclaw_bool time_interp,
                   fclaw_bool read_parallel_patches)
{
    fclaw2d_exchange_info_t e_info;
    e_info.time_interp = time_interp; /* Does this matter here? */
    e_info.read_parallel_patches = read_parallel_patches;
    e_info.exchange_type = FCLAW2D_AVERAGE;

    /* Only update ghost cells at local boundaries */
    e_info.grid_type = FCLAW2D_IS_COARSE;

    /* Face average */
    fclaw2d_domain_iterate_level(domain, coarse_level,
                                 cb_face_fill, (void *) &e_info);

    /* Corner average */
    fclaw2d_domain_iterate_level(domain, coarse_level, cb_corner_fill,
                                 (void *) &e_info);

    if (read_parallel_patches)
    {
        /* Second pass : average from local fine grids to remote coarse grids. These
           coarse grids might be needed for interpolation later. */
        e_info.grid_type = FCLAW2D_IS_FINE;

        int fine_level = coarse_level + 1;

        /* Face average */
        fclaw2d_domain_iterate_level(domain, fine_level,
                                     cb_face_fill, (void *) &e_info);

        /* Corner average :
           We can skip the corner update, since we don't need the corner ghost cell
           values for doing interpolation (at least not yet) */
#if 0
        fclaw2d_domain_iterate_level(domain, fine_level, cb_corner_average,
                                     (void *) &e_info);
#endif
    }
}

static
void interpolate2ghost(fclaw2d_domain_t *domain,
                       int fine_level,
                       fclaw_bool time_interp,
                       fclaw_bool read_parallal_patches)
{
    fclaw2d_exchange_info_t e_info;
    e_info.time_interp = time_interp;
    e_info.level = fine_level;
    e_info.exchange_type = FCLAW2D_INTERPOLATE;

    int coarse_level = fine_level - 1;

    /* ----------------------------------------------------------
       First pass - look for fine grids to interpolate to. This
       should include include the time interpolated level.
       ---------------------------------------------------------- */

    e_info.grid_type = FCLAW2D_IS_COARSE;
    e_info.read_parallel_patches = read_parallal_patches;

    /* Face interpolate */
    fclaw2d_domain_iterate_level(domain,coarse_level, cb_face_fill,
                                 (void *) &e_info);

    /* Corner interpolate */
    fclaw2d_domain_iterate_level(domain,coarse_level, cb_corner_fill,
                                 (void *) &e_info);

    /* -----------------------------------------------------
       Second pass - Iterate over local fine grids, looking
       for remote coarse grids we can use to fill in BCs at
       fine grid ghost cells along the parallel boundary
       ----------------------------------------------------- */

    e_info.grid_type = FCLAW2D_IS_FINE;

    /* Interpolate to faces at parallel boundaries from coarse grid ghost
       patches */
    fclaw2d_domain_iterate_level(domain, fine_level, cb_face_fill,(void *) &e_info);

    /* Interpolate to corners at parallel boundaries from coarse grid
       ghost patches */
    fclaw2d_domain_iterate_level(domain, fine_level, cb_corner_fill,(void *) &e_info);
}


/* -----------------------------------------------------------------------
   Physical boundaries conditions and misc.
   ---------------------------------------------------------------------*/
static
void fill_physical_ghost(fclaw2d_domain_t* domain,
                         int minlevel,
                         int maxlevel,
                         double t, fclaw_bool time_interp)
{
    for(int level = maxlevel; level >= minlevel; level--)
    {
        /* No time interpolation yet */
        fclaw2d_set_physical_bc(domain,level,t,0);
    }
    if (time_interp)
    {
        /* Fill boundary ghost cells on time interpolated level */
        int time_interp_level = minlevel-1;
        fclaw2d_set_physical_bc(domain,time_interp_level,t,1);
    }
}

/* -----------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------------*/

/**
 * <summary>Complete exchange of all ghost patches at all levels.</summary>
 * <remarks>All parallel ghost patches are also exchanged at all
 * levels.</remarks>
 * <list>
 *    <item>Every level exchanges ghost cells with other patches
 *       at that level</item>
 *    <item>Every finer level exchanges with a coarser level</item>
 *    <item> All levels will be updated in next update step, regardless of
 *       whether we are in the subcycled or non-subcycled case.</item>
 *       </list>
 *   The reason for two separate ghost cell exchange routines is that
 *   the logic here is considerably simpler than for the partial
 *   update used in intermediate steps in the subcycled case.
 **/
void fclaw2d_ghost_update(fclaw2d_domain_t* domain,
                          int minlevel,
                          int maxlevel,
                          int time_interp,
                          fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    fclaw_global_infof("Exchanging ghost patches across all levels\n");

    double t = get_domain_time(domain);

#if 0
    /* uncomment this if debugging ghost cell interpolation */
    fclaw_global_essentialf("WARNING : compute_slopes set to average\n");
#endif

    /* ---------------------------------------------------------
       Get coarse grid ghost cells ready to use for interpolation.
       Coarse ghost cells on ghost patches are not updated in
       this step.  Ghost patches do not yet have valid data, and
       so boundary patches will have to be updated after the exchange.
       ---------------------------------------------------------- */
    /* Copy ghost cells at coarse levels.  Include finest level, although it isn't
       needed immediately */
    fclaw_bool read_parallel_patches = fclaw_false;
    copy_ghost_samelevel(domain,minlevel,maxlevel,time_interp,
                         read_parallel_patches);

    if (time_interp)
    {
        int time_interp_level = minlevel - 1;
        fclaw_global_infof("Time interpolated level is %d\n",   \
                           time_interp_level);
    }

    int mincoarse = minlevel;
    int maxcoarse = maxlevel-1;   /* maxlevel >= minlevel */

    /* Average fine grid data to the coarse grid. If we only have a single
     grid, we have maxcoarse < mincoarse, and this loop does nothing */
    average_fine2coarse_ghost(domain,mincoarse,maxcoarse,
                              time_interp,
                              read_parallel_patches);

    /* Supply physical boundary conditions on coarse grid;  We don't do the
       physical BCs for the finest level grids */
    fill_physical_ghost(domain,mincoarse,maxcoarse,t,time_interp);

    if (domain->mpisize > 1)
    {
        /* -------------------------------------------------------------
           Parallel ghost patch exchange
           ------------------------------------------------------------- */
        fclaw2d_exchange_ghost_patches(domain,minlevel,maxlevel,time_interp);

        /* -------------------------------------------------------------
           Repeat above, but now with parallel ghost cells.
           ------------------------------------------------------------- */

        /* Fill in ghost cells on parallel patch boundaries */
        read_parallel_patches = fclaw_true;

        copy_ghost_samelevel(domain,minlevel,maxlevel,time_interp,
                             read_parallel_patches);

        /* Average fine grid data to the coarse grid. */
        average_fine2coarse_ghost(domain,mincoarse,maxcoarse, time_interp,
                                  read_parallel_patches);

        /* Supply physical boundary conditions on coarse grid */
        fill_physical_ghost(domain,mincoarse,maxcoarse,t,time_interp);
    }

    /* ---------------------------------------------------------
       Fill in fine grid data via interpolation.
       Fine grids that live at parallel boundaries may need to get
       ghost cell data from the coarse ghost patch.
       ---------------------------------------------------------- */
    int minfine = minlevel+1;
    int maxfine = maxlevel;

    read_parallel_patches = fclaw_true;
    interpolate_coarse2fine_ghost(domain,minfine, maxfine,
                                  time_interp,read_parallel_patches);

    /* --------------------------------------------------------- */
    /* Do a final fill in of boundary conditions of all physical
       values */
    fill_physical_ghost(domain,minlevel,maxlevel,t,time_interp);

    /* --------------------------------------------------------- */

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE)
    {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}



#ifdef __cplusplus
#if 0
{
#endif
}
#endif
