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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_partition.h>

#include <sc_statistics.h>

#define FCLAW2D_STATS_SET(stats,ddata,NAME) do {                               \
    SC_CHECK_ABORT (!(ddata)->timers[FCLAW2D_TIMER_ ## NAME].running,          \
                    "Timer " #NAME " still running in amrreset");              \
    sc_stats_set1 ((stats) + FCLAW2D_TIMER_ ## NAME,                           \
                   (ddata)->timers[FCLAW2D_TIMER_ ## NAME].cumulative, #NAME); \
} while (0)


static
void delete_ghost_patches(fclaw2d_domain_t* domain)
{
    for(int i = 0; i < domain->num_ghost_patches; i++)
    {
        fclaw2d_patch_t* ghost_patch = &domain->ghost_patches[i];

        fclaw2d_patch_delete_cp(ghost_patch);
        fclaw2d_patch_delete_data(ghost_patch);
    }
}


void amrreset(fclaw2d_domain_t **domain)
{
    fclaw2d_domain_data_t *ddata = get_domain_data (*domain);

    for(int i = 0; i < (*domain)->num_blocks; i++)
    {
        fclaw2d_block_t *block = (*domain)->blocks + i;
        fclaw2d_block_data_t *bd = (fclaw2d_block_data_t *) block->user;

        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_delete_cp(patch);
            fclaw2d_patch_delete_data(patch);
#if 0
            fclaw2d_patch_data_t *pdata = (fclaw2d_patch_data_t *) patch->user;
            delete pdata->cp;
            pdata->cp = NULL;
            FCLAW2D_FREE (pdata);
            patch->user = NULL;
#endif
            ++ddata->count_delete_clawpatch;
        }

        FCLAW2D_FREE (bd);
        block->user = NULL;
    }

    // Free old parallel ghost patch data structure, must exist by construction.
    delete_ghost_patches(*domain);
    fclaw2d_domain_exchange_t *e_old = fclaw2d_partition_get_exchange_data(*domain);
    fclaw2d_domain_free_after_exchange (*domain, e_old);

    // Output memory discrepancy for the ClawPatch
    if (ddata->count_set_clawpatch != ddata->count_delete_clawpatch) {
        printf ("[%d] This domain had Clawpatch set %d and deleted %d times\n",
                (*domain)->mpirank,
                ddata->count_set_clawpatch, ddata->count_delete_clawpatch);
    }

    // Evaluate timers if this domain has not been superseded yet.
    if (ddata->is_latest_domain) {
        sc_statinfo_t stats[FCLAW2D_TIMER_COUNT];

        fclaw2d_timer_stop (&ddata->timers[FCLAW2D_TIMER_WALLTIME]);

        FCLAW2D_STATS_SET (stats, ddata, INIT);
        FCLAW2D_STATS_SET (stats, ddata, REGRID);
        FCLAW2D_STATS_SET (stats, ddata, OUTPUT);
        FCLAW2D_STATS_SET (stats, ddata, CHECK);
        FCLAW2D_STATS_SET (stats, ddata, ADVANCE);
        FCLAW2D_STATS_SET (stats, ddata, EXCHANGE);
        FCLAW2D_STATS_SET (stats, ddata, BUILDPATCHES);
        FCLAW2D_STATS_SET (stats, ddata, WALLTIME);
        sc_stats_set1 (&stats[FCLAW2D_TIMER_UNACCOUNTED],
                       ddata->timers[FCLAW2D_TIMER_WALLTIME].cumulative -
                       (ddata->timers[FCLAW2D_TIMER_INIT].cumulative +
                        ddata->timers[FCLAW2D_TIMER_REGRID].cumulative +
                        ddata->timers[FCLAW2D_TIMER_OUTPUT].cumulative +
                        ddata->timers[FCLAW2D_TIMER_CHECK].cumulative +
                        ddata->timers[FCLAW2D_TIMER_ADVANCE].cumulative +
                        ddata->timers[FCLAW2D_TIMER_EXCHANGE].cumulative),
                       "UNACCOUNTED");
        sc_stats_compute ((*domain)->mpicomm, FCLAW2D_TIMER_COUNT, stats);
        sc_stats_print (sc_package_id, SC_LP_ESSENTIAL, FCLAW2D_TIMER_COUNT,
                        stats, 1, 0);
        SC_GLOBAL_ESSENTIALF ("Procs %d advance %d %g exchange %d %g "
                               "regrid %d %g\n", (*domain)->mpisize,
                               ddata->count_amr_advance,
                               stats[FCLAW2D_TIMER_ADVANCE].average,
                               ddata->count_ghost_exchange,
                               stats[FCLAW2D_TIMER_EXCHANGE].average,
                               ddata->count_amr_regrid,
                               stats[FCLAW2D_TIMER_REGRID].average);
        SC_GLOBAL_ESSENTIALF ("Max/P %d advance %d %g exchange %d %g "
                               "regrid %d %g\n", (*domain)->mpisize,
                               ddata->count_amr_advance,
                               stats[FCLAW2D_TIMER_ADVANCE].max,
                               ddata->count_ghost_exchange,
                               stats[FCLAW2D_TIMER_EXCHANGE].max,
                               ddata->count_amr_regrid,
                               stats[FCLAW2D_TIMER_REGRID].max);
    }


    delete_domain_data(*domain);  // Delete allocated pointers to set of functions.

    fclaw2d_domain_destroy(*domain);
    *domain = NULL;
}
