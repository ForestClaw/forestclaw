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
/** 
 * @file
 * @brief Functions needed for an underlying wrapped fclaw3d_domain_t in fclaw_domain_t
 *        Most users will not need to include the file, unless doing something more advanced.
 */


#ifndef FCLAW3D_DOMAIN_H
#define FCLAW3D_DOMAIN_H

#include <forestclaw3d.h>  /* Contains definition of patch-iterator callback */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_domain;

typedef struct fclaw3d_domain_wrap
{
    fclaw3d_domain_t*          domain;
    fclaw3d_domain_exchange_t* exchange;
    fclaw3d_domain_indirect_t* indirect;
} fclaw3d_domain_wrap_t;

struct fclaw_domain* fclaw_domain_wrap_3d(fclaw3d_domain_t *domain2d);

fclaw3d_domain_t* fclaw_domain_get_3d_domain(const struct fclaw_domain* patch);

fclaw3d_domain_wrap_t* fclaw_domain_get_3d_domain_wrap(struct fclaw_domain* domain);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
