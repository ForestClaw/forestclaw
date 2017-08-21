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

#include <fclaw2d_time_sync.h>

void fclaw2d_time_sync_reset(fclaw2d_global_t *glob, 
                             int minlevel,int maxlevel)
{
    /* This is used for updating conservation arrays, for example */
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt();

    if (fclaw_vt->time_sync_reset != NULL)
    {
        fclaw_vt->time_sync_reset(glob,minlevel,maxlevel);
    }
}


static
void copy_at_blockbdry(fclaw2d_global_t *glob,
                       int level,
                       int read_parallel_patches,
                       fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
    fclaw2d_exchange_info_t e_info;
    e_info.exchange_type = FCLAW2D_TIME_SYNC_COPY;
    e_info.grid_type = FCLAW2D_IS_COARSE;
    e_info.time_interp = 0;
    e_info.read_parallel_patches = read_parallel_patches;

    fclaw2d_global_patch_iterator(glob, level, cb_face_fill,
                         &e_info);    
}


static
void fine2coarse(fclaw2d_global_t *glob,
                int level,
                int read_parallel_patches,
                fclaw2d_ghost_fill_parallel_mode_t ghost_mode)
{
    fclaw2d_exchange_info_t e_info;
    e_info.exchange_type = FCLAW2D_TIME_SYNC_FINE_TO_COARSE;
    e_info.grid_type = FCLAW2D_IS_COARSE;
    e_info.time_interp = 0;
    e_info.read_parallel_patches = read_parallel_patches;

    fclaw2d_global_patch_iterator(glob, level, cb_face_fill,
                         &e_info);
}

static
void coarse_correct(fclaw2d_global_t *glob, 
                    int minlevel, int maxlevel,
                    int read_parallel_patches,
                    fclaw2d_ghost_fill_parallel_mode_t ghost_mode)

{
    int level;

    int read_parallel_patches;
    int time_interp = 0;

    int mincoarse = minlevel;
    int maxcoarse = maxlevel-1;
    
    for(level = maxcoarse; level >= mincoarse; level--)
    {
        fine2coarse(glob,level,
                    read_parallel_patches,ghost_mode);
    }

    for(level = maxlevel; level >= minlevel; level--)
    {
        copy_at_blockbdry(glob,level,
                          read_parallel_patches,ghost_mode);
    }

}


void fclaw2d_time_sync(fclaw2d_global_t *glob, int minlevel, int maxlevel)
{

    int time_interp = 0;

    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_TIMESYNC]);


    /* --------------------------------------------------------------
        Start send ...
    ------------------------------------------------------------*/
    fclaw2d_exchange_ghost_patches_begin(glob,minlevel,maxlevel,time_interp,
                                         FCLAW2D_TIMER_TIMESYNC);

    /* -------------------------------------------------------------
        Receive ghost patches ...
    ------------------------------------------------------------- */

    fclaw2d_exchange_ghost_patches_end(glob,minlevel,maxlevel,time_interp,
                                       FCLAW2D_TIMER_TIMESYNC);

    /* -------------------------------------------------------------
        Loop over ghost patches to find indirect neighbors and do
        any necessary face exchanges.
    ------------------------------------------------------------- */
    fclaw2d_face_neighbor_ghost(glob,minlevel,maxlevel,time_interp);


    /* Add time_sync coarse grid corrections (include corrections form 
    the fine grid and from copies) */
    coarse_correct(glob,minlevel,maxlevel,read_parallel_patches,parallel_mode);

    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_TIMESYNC]);

}


