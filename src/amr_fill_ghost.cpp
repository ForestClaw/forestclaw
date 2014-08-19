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

#include "amr_includes.H"

void intralevel_ghost_copy(fclaw2d_domain_t* domain, int minlevel,
                           int maxlevel,
                           fclaw2d_ghost_cell_types_t ghost_cells)
{
    /* Copy between grids that are at the same level. */
    for(int level = maxlevel; level >= minlevel; level--)
    {
        copy2ghost(domain,level,ghost_cells);
    }
}


void fill_coarse_ghost(fclaw2d_domain_t *domain,
                       int mincoarse,
                       int maxcoarse,
                       fclaw_bool time_interp,
                       fclaw2d_ghost_cell_types_t ghost_cells)
{
    /* Average fine grids to coarse grid ghost cells */
    for(int level = maxcoarse; level >= mincoarse; level--)
    {
        average2ghost(domain,level,time_interp,ghost_cells);
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
        set_phys_bc(domain,level,t,time_interp);
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



/* ------------------------------------------------------------
   This does a complete exchange over all levels.

   -- Parallel ghost patches are exchanged at all levels
   -- Every level exchanges ghost cells with other patches
      at that level
   -- Every finer level exchanges with a coarser level
   -- No time interpolation is assumed, as all levels are time
      synchronized at this point.
   -- This the only routine that is called for the non-subcycled
      case.
   -- All levels will be updated in next update step, regardless of
      whether we are in the subcycled or non-subcycled case.

  The reason for two separate ghost cell exchange routines is that
  the logic here is considerably simpler than for the partial
  update used in intermediate steps in the subcycled case.
  ------------------------------------------------------------- */
void update_ghost_all_levels(fclaw2d_domain_t* domain,
                             fclaw2d_timer_names_t running)
{
    fclaw2d_domain_data_t *ddata = get_domain_data(domain);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_stop (&ddata->timers[running]);
    }
    fclaw2d_timer_start (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);

    const amr_options_t *gparms = get_domain_parms(domain);
    fclaw_bool verbose = gparms->verbosity;

    if (verbose)
    {
        cout << "Exchanging ghost patches across all levels " << endl;
    }

    int minlevel = domain->global_minlevel;
    int maxlevel = domain->global_maxlevel;

    /* No time interpolation, since levels are time synchronized */
    fclaw_bool no_time_interp = fclaw_false;
    double t = get_domain_time(domain);

    /* ---------------------------------------------------------
       Get coarse grid ghost cells ready to use for interpolation.
       Coarse ghost cells on ghost patches are not updated in
       this step.  Ghost patches do not yet have valid data, and
       so boundary patches will have to be updated after the exchange.
       ---------------------------------------------------------- */
    int mincoarse = minlevel;
    int maxcoarse = maxlevel-1;

    /* Copy ghost cells at coarse levels.  Include finest level, although it isn't
       needed immediately */
    intralevel_ghost_copy(domain,minlevel,maxlevel,FCLAW2D_LOCAL_BOUNDARY);

    /* Average fine grid data to the coarse grid. */
    fill_coarse_ghost(domain,mincoarse,maxcoarse, no_time_interp,
                      FCLAW2D_LOCAL_BOUNDARY);

    /* Supply physical boundary conditions on coarse grid */
    fill_physical_ghost(domain,mincoarse,maxcoarse,t,no_time_interp);

    /* -------------------------------------------------------------
      Parallel ghost patch exchange
      -------------------------------------------------------------- */
    exchange_ghost_patch_data_all(domain);

    /* -------------------------------------------------------------
      Repeat above, but now only with ghost cells.
      -------------------------------------------------------------- */

    /* Fill in ghost cells on parallel patch boundaries */
    intralevel_ghost_copy(domain,minlevel,maxlevel,
                          FCLAW2D_PARALLEL_BOUNDARY);

    /* Average fine grid data to the coarse grid. */
    fill_coarse_ghost(domain,mincoarse,maxcoarse, no_time_interp,
                      FCLAW2D_PARALLEL_BOUNDARY);

    /* Supply physical boundary conditions on coarse grid */
    fill_physical_ghost(domain,mincoarse,maxcoarse,t,no_time_interp);


    /* ---------------------------------------------------------
       Fill in fine grid data via interpolation.  Fine grid ghost
       cells on ghost patches are updated.  They must then be
       passed back to their home processor.
       ---------------------------------------------------------- */
    int minfine = minlevel+1;
    int maxfine = maxlevel;

    fill_fine_ghost(domain,minfine, maxfine, no_time_interp);

    /* Fill in fine grid ghost cells at physical boundaries. This is
     needed to make sure we have valid values in the fine grid corners at
     the physical boundaries. */
    fill_physical_ghost(domain,minfine,maxfine,t, no_time_interp);

    /* --------------------------------------------------------- */

    // Stop timing
    fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_EXCHANGE]);
    if (running != FCLAW2D_TIMER_NONE) {
        fclaw2d_timer_start (&ddata->timers[running]);
    }
}
