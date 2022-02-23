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

#include <fclaw_global.h>
#include <fclaw2d_global.h>
#include <fclaw3d_global.h>

void
fclaw_global_iterate_patches (fclaw_global_t * glob,
                              fclaw_patch_callback_t gpcb, void *user)
{
    fclaw_global_iterate_t gi;
    gi.glob = glob;
    gi.gpcb = gpcb;
    gi.gfcb = NULL;
    gi.user = user;

    if (glob->domain->dim == 2)
    {
        fclaw2d_domain_iterate_patches (glob->domain->d.d2.domain2,
                                        fclaw2d_iterate_patch_cb, &gi);
    }
    else
    {
        FCLAW_ASSERT (glob->domain->dim == 3);
        fclaw3d_domain_iterate_patches (glob->domain->d.d3.domain3,
                                        fclaw3d_iterate_patch_cb, &gi);
    }
}

void
fclaw_global_iterate_families (fclaw_global_t * glob,
                               fclaw_family_callback_t gfcb, void *user)
{
    fclaw_global_iterate_t gi;
    gi.glob = glob;
    gi.gpcb = NULL;
    gi.gfcb = gfcb;
    gi.user = user;

    if (glob->domain->dim == 2)
    {
        fclaw2d_domain_iterate_families (glob->domain->d.d2.domain2,
                                         fclaw2d_iterate_family_cb, &gi);
    }
    else
    {
        FCLAW_ASSERT (glob->domain->dim == 3);
        fclaw3d_domain_iterate_families (glob->domain->d.d3.domain3,
                                         fclaw3d_iterate_family_cb, &gi);
    }
}
