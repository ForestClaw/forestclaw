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

/** \file
 * Fill ghost cells.
 *
 *
 **/

#include "amr_includes.H"

#ifdef __cplusplus
extern "C" {
#if 0
}
#endif
#endif

static
void interpolate2ghost(fclaw2d_domain_t *domain,int fine_level,
                       fclaw_bool time_interp);

static
void copy2ghost(fclaw2d_domain_t *domain, int level,
                fclaw_bool read_parallel_patches);

static
void average2ghost(fclaw2d_domain_t *domain, int coarse_level,
                   fclaw_bool time_interp,
                   fclaw_bool read_parallel_patches);


void intralevel_ghost_copy(fclaw2d_domain_t* domain, int minlevel,
                           int maxlevel,
                           fclaw_bool read_parallel_patches)
{
    /* Copy between grids that are at the same level. */
    for(int level = maxlevel; level >= minlevel; level--)
    {
        copy2ghost(domain,level,read_parallel_patches);
    }
}


/**
 * \ingroup Averaging
 * Fill in coarse grid ghost cells by averaging or copying  from neighboring fine grids.
 **/
void fill_coarse_ghost(fclaw2d_domain_t *domain,
                       int mincoarse,
                       int maxcoarse,
                       fclaw_bool time_interp,
                       fclaw_bool read_parallel_patches)
{
    /* Average fine grids to coarse grid ghost cells */
    for(int level = maxcoarse; level >= mincoarse; level--)
    {
        average2ghost(domain,level,time_interp,read_parallel_patches);
    }
}

void fill_fine_ghost(fclaw2d_domain_t* domain, int minfine, int maxfine,
                     fclaw_bool time_interp)
{
    /* Interpolate from coarse grid to fine grid ghost */
    for(int level = maxfine; level >= minfine; level--)
    {
        interpolate2ghost(domain,level,time_interp);
    }
}

void fill_physical_ghost(fclaw2d_domain_t* domain, int minlevel, int maxlevel,
                         double t, fclaw_bool time_interp)
{
    for(int level = maxlevel; level >= minlevel; level--)
    {
        fclaw2d_set_physical_bc(domain,level,t,time_interp);
    }
}

fclaw_bool fclaw2d_patch_is_remote(fclaw2d_patch_t *patch)
{
    return fclaw2d_patch_is_ghost(patch);
}

fclaw_bool fclaw2d_patch_is_local(fclaw2d_patch_t *patch)
{
    return !fclaw2d_patch_is_ghost(patch);
}

static
void copy2ghost(fclaw2d_domain_t *domain, int level,
                fclaw_bool read_parallel_patches)
{
    fclaw2d_exchange_info_t e_info;
    e_info.exchange_type = FCLAW2D_COPY;
    e_info.grid_type = FCLAW2D_IS_COARSE;
    e_info.time_interp = fclaw_false;
    e_info.read_parallel_patches = read_parallel_patches;

    /* face exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_face_fill,
                                 (void *) &e_info);

    /* corner exchanges */
    fclaw2d_domain_iterate_level(domain, level, cb_corner_fill,
                                 (void *) &e_info);

}


static
void average2ghost(fclaw2d_domain_t *domain, int coarse_level,
                   fclaw_bool time_interp,
                   fclaw_bool read_parallel_patches)
{
    fclaw2d_exchange_info_t e_info;
    e_info.time_interp = time_interp;
    e_info.read_parallel_patches = read_parallel_patches;
    e_info.exchange_type = FCLAW2D_AVERAGE;
    int fine_level = coarse_level + 1;

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

        /* Face average */
        fclaw2d_domain_iterate_level(domain, fine_level,
                                     cb_face_fill, (void *) &e_info);

        /* Corner average */
        /* We can skip the corner update, since we don't need the corner ghost cell
           values for doing interpolation */
        /*
        fclaw2d_domain_iterate_level(domain, fine_level, cb_corner_average,
                                     (void *) &e_info);
        */
    }
}

static
void interpolate2ghost(fclaw2d_domain_t *domain,int fine_level,
                       fclaw_bool time_interp)
{
    fclaw2d_exchange_info_t e_info;
    e_info.time_interp = time_interp;
    e_info.level = fine_level;
    e_info.exchange_type = FCLAW2D_INTERPOLATE;

    int coarse_level = fine_level - 1;

    /* ----------------------------------------------------------
       First pass - iterate over coarse grids and update ghost
       cells on local fine grids.
       ---------------------------------------------------------- */

    /* Interpolation done over all patches */
    e_info.read_parallel_patches = fclaw_true;
    e_info.grid_type = FCLAW2D_IS_COARSE;

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
    fclaw2d_domain_iterate_level(domain, fine_level, cb_face_fill,
                                 (void *) &e_info);

    /* Interpolate to corners at parallel boundaries from coarse grid
       ghost patches */
    fclaw2d_domain_iterate_level(domain, fine_level, cb_corner_fill,
                                 (void *) &e_info);
}





/**
 * <summary>Complete exchange of all ghost patches at all levels.</summary>
 * <remarks>All parallel ghost patches are also exchanged at all
 * levels.</remarks>
 * <list>
 *    <item>Every level exchanges ghost cells with other patches
 *       at that level</item>
 *    <item>Every finer level exchanges with a coarser level</item>
 *    <item>No time interpolation is assumed, as all levels are time
 *       synchronized at this point.</item>
 *    <item>This the only routine that is called for the non-subcycled
 *       case.</item>
 *    <item> All levels will be updated in next update step, regardless of
 *       whether we are in the subcycled or non-subcycled case.</item>
 *       </list>
 *   The reason for two separate ghost cell exchange routines is that
 *   the logic here is considerably simpler than for the partial
 *   update used in intermediate steps in the subcycled case.
 **/
void fclaw2d_ghost_update_all_levels(fclaw2d_domain_t* domain,
                                     fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    fclaw_global_infof("Exchanging ghost patches across all levels\n");

    int minlevel = domain->global_minlevel;
    int maxlevel = domain->global_maxlevel;

    /* No time interpolation, since levels are time synchronized */
    fclaw_bool time_interp = fclaw_false;
    double t = get_domain_time(domain);

    /* ---------------------------------------------------------
       Get coarse grid ghost cells ready to use for interpolation.
       Coarse ghost cells on ghost patches are not updated in
       this step.  Ghost patches do not yet have valid data, and
       so boundary patches will have to be updated after the exchange.
       ---------------------------------------------------------- */
    int mincoarse = minlevel;
    int maxcoarse = maxlevel-1;   /* Be careful with minlevel-maxlevel */

    /* Copy ghost cells at coarse levels.  Include finest level, although it isn't
       needed immediately */
    fclaw_bool read_parallel_patches = fclaw_false;
    intralevel_ghost_copy(domain,minlevel,maxlevel,read_parallel_patches);

    /* Average fine grid data to the coarse grid. */
    fill_coarse_ghost(domain,mincoarse,maxcoarse, time_interp,
                      read_parallel_patches);

    /* Supply physical boundary conditions on coarse grid;  We don't do the
       physical BCs for the finest level grids */
    fill_physical_ghost(domain,mincoarse,maxcoarse,t,time_interp);

    if (domain->mpisize > 1)
    {
        /* -------------------------------------------------------------
           Parallel ghost patch exchange
           -------------------------------------------------------------- */
        exchange_ghost_patch_data_all(domain);

        /* -------------------------------------------------------------
           Repeat above, but now with parallel ghost cells.
           This may involve lots of repeated worked.
           -------------------------------------------------------------- */

        /* Fill in ghost cells on parallel patch boundaries */
        read_parallel_patches = fclaw_true;
        intralevel_ghost_copy(domain,minlevel,maxlevel,
                              read_parallel_patches);

        /* Average fine grid data to the coarse grid. */
        fill_coarse_ghost(domain,mincoarse,maxcoarse, time_interp,
                          read_parallel_patches);

        /* Supply physical boundary conditions on coarse grid */
        fill_physical_ghost(domain,mincoarse,maxcoarse,t,time_interp);
    }

    /* ---------------------------------------------------------
       Fill in fine grid data via interpolation.  Fine grid ghost
       cells on ghost patches are updated.  They must then be
       passed back to their home processor.
       ---------------------------------------------------------- */
    int minfine = minlevel+1;
    int maxfine = maxlevel;

    fill_fine_ghost(domain,minfine, maxfine,time_interp);

    /* --------------------------------------------------------- */
    /* Do a final fill in of boundary conditions of all physical
       values */
    fill_physical_ghost(domain,minlevel,maxlevel,t,time_interp);

    /* --------------------------------------------------------- */

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}

#if 0

/* ----------------------------------------------------------
   This does the intermediate ghost patch and ghost cell
   exchange needed by either the subcycled case or the non
   subcycled case.

   Assumptions going into this routine :
   -- coarse_level > global_minlevel
   -- Data at coarse_level-1 is time interpolated data
   -- This routine is never called in the non-subcycled case
   -- This should be more 'lightweight' than doing a complete
      exchange (across all levels) at each fine grid time
      step.  But it should be equivalent.

   ---------------------------------------------------------- */
void update_ghost_partial(fclaw2d_domain_t* domain, int coarse_level,
                          int fine_level,
                          subcycle_manager *a_time_stepper,
                          fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = gparms->verbosity;

    int time_interp_level = coarse_level - 1;

    fclaw_global_infof("Exchanging ghost patches from levels %d to %d\n",\
                       coarse_level, fine_level);
    if (!a_time_stepper->nosubcycle())
    {
        fclaw_global_infof("Time interpolated level is %d\n",   \
                           time_interp_level);
    }

    /* Make available patches from levels coarse to fine */
    set_boundary_patch_ptrs(domain,time_interp_level, fine_level);

    /* Do parallel ghost exchange */
    exchange_ghost_patch_data_levels(domain,time_interp_level,
                                     fine_level);

    /* -------------------------------------------------------
       Do ghost cell exchange.
       ------------------------------------------------------- */
    double t = get_domain_time(domain); /* needed for phys. bc */
    /* Do all level exchanges first */
    for(int level = fine_level; level >= coarse_level; level--)
    {
        level_exchange(domain,level);
    }

    /* Then do exchanges with coarser level */
    fclaw_bool time_interp_vec[fine_level-coarse_level+1];
    for(int level = fine_level; level > coarse_level; level--)
    {
        time_interp_vec[level] = fclaw_false;
    }

    /* coarse level exchanges with a time interpolated level
       time_interp_level = coarse_level-1
    */
    time_interp_vec[coarse_level] = fclaw_true;
    for(int level = fine_level; level >= coarse_level; level--)
    {
        exchange_with_coarse(domain,level,t,time_interp_vec[level]);
    }


    fclaw_bool time_interp = fclaw_false;
    for (int level = fine_level; level >= coarse_level; level--)
    {
        fclaw2d_set_physical_bc(domain,level,t,time_interp);
    }

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}
#endif

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
