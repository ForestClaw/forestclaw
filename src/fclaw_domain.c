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

#include <fclaw_domain.h>
#include <fclaw2d_domain.h>
#include <fclaw3d_domain.h>

int
fclaw_domain_is_valid (fclaw_domain_t * domain)
{
    if (domain == NULL)
        return 0;
    if (domain->dim == 2)
    {
        if (domain->d.d2.dmagic2 != FCLAW2D_DOMAIN_MAGIC)
            return 0;
        if (domain->d.d2.domain2 == NULL)
            return 0;
        if (domain->d.d2.domain2->mpisize <= 0)
            return 0;
    }
    else
    {
        if (domain->dim != 3)
            return 0;
        if (domain->d.d3.dmagic3 != FCLAW3D_DOMAIN_MAGIC)
            return 0;
        if (domain->d.d3.domain3 == NULL)
            return 0;
        if (domain->d.d3.domain3->mpisize <= 0)
            return 0;
    }
    return 1;
}

void
fclaw_domain_destroy (fclaw_domain_t * domain,
                      fclaw_domain_callback_t dele, void *user)
{
    FCLAW_ASSERT (fclaw_domain_is_valid (domain));
    if (domain->dim == 2)
    {
        fclaw_domain_destroy2d (domain, dele, user);
    }
    else
    {
        FCLAW_ASSERT (domain->dim == 3);
        fclaw_domain_destroy3d (domain, dele, user);
    }
}

void
fclaw_domain_iterate_patches (fclaw_domain_t * domain,
                              fclaw_domain_callback_t iter, void *user)
{
    fclaw_domain_iterate_t di;

    FCLAW_ASSERT (fclaw_domain_is_valid (domain));

    di.d = domain;
    di.iter = iter;
    di.user = user;

    if (domain->dim == 2)
    {
        fclaw2d_domain_iterate_patches (domain->d.d2.domain2,
                                        fclaw2d_domain_iterate_cb, &di);
    }
    else
    {
        FCLAW_ASSERT (domain->dim == 3);
        fclaw3d_domain_iterate_patches (domain->d.d3.domain3,
                                        fclaw3d_domain_iterate_cb, &di);
    }
}

void
fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                            fclaw_domain_callback_t iter, void *user)
{
    fclaw_domain_iterate_t di;

    FCLAW_ASSERT (fclaw_domain_is_valid (domain));

    di.d = domain;
    di.iter = iter;
    di.user = user;

    if (domain->dim == 2)
    {
        fclaw2d_domain_iterate_level (domain->d.d2.domain2, level,
                                      fclaw2d_domain_iterate_cb, &di);
    }
    else
    {
        FCLAW_ASSERT (domain->dim == 3);
        fclaw3d_domain_iterate_level (domain->d.d3.domain3, level,
                                      fclaw3d_domain_iterate_cb, &di);
    }
}
