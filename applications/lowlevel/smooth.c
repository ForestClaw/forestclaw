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

#include <fclaw_base.h>
#include <fclaw2d_convenience.h>

typedef struct fclaw_smooth
{
    sc_MPI_Comm mpicomm;
    fclaw_app_t *a;

    int minlevel;
    int maxlevel;
    fclaw2d_domain_t *domain;
}
fclaw_smooth_t;

static void
run_refine (fclaw_smooth_t * s)
{
    int lev;
    int ib, ip;
    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    for (lev = s->minlevel; lev < s->maxlevel; ++lev)
    {
        for (ib = 0; ib < s->domain->num_blocks; ++ib)
        {
            block = s->domain->blocks + ib;
            for (ip = 0; ip < block->num_patches; ++ip)
            {
                patch = block->patches + ip;
                if ((patch->xlower > .3 && patch->xupper < .6) &&
                    (patch->ylower > .2 && patch->yupper < .55))
                {

                    /* prompt refinement of this patch */
                    fclaw2d_patch_mark_refine (s->domain, ib, ip);
                }
            }
        }
    }
}

int
main (int argc, char **argv)
{
    fclaw_smooth_t smoo, *s = &smoo;

    s->a = fclaw_app_new (&argc, &argv, NULL);
    s->mpicomm = fclaw_app_get_mpi_size_rank (s->a, NULL, NULL);

    /* set parameters */
    s->minlevel = 2;
    s->maxlevel = 5;

    /* create a new domain */
    s->domain = fclaw2d_domain_new_unitsquare (s->mpicomm, s->minlevel);

    /* run refinement cycles */
    run_refine (s);

    /* cleanup */
    fclaw2d_domain_destroy (s->domain);
    fclaw_app_destroy (s->a);
    return 0;
}
