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

#ifndef FCLAW_GLOBAL_H
#define FCLAW_GLOBAL_H

#include <forestclaw.h>

typedef struct fclaw_global
{
  fclaw_domain_t *domain;

  /* lots of dimension-independent stuff */
}
fclaw_global_t;

/** Dimension-independent callback prototype for patch iterators.
 * We iterate over local patches only.
 * \param [in] domain	Dimension-independent global context.
 * \param [in] patch	The local patch currently processed by the iterator.
 * \param [in] blockno  Block number of processed patch.
 * \param [in] patchno  Patch number within block of processed patch.
 * \param [in,out] user	Data that was passed into the iterator functions.
 */
typedef void (*fclaw_global_callback_t)
    (fclaw_global_t * global, fclaw_patch_t * patch,
     int blockno, int patchno, void *user);

typedef struct fclaw_global_iterate
{
  fclaw_global_t *glob;
  fclaw_global_callback_t gpcb;
  void *user;
}
fclaw_global_iterate_t;

void fclaw_global_iterate_patches (fclaw_global_t * glob,
                                   fclaw_global_callback_t gpcb, void *user);

void fclaw_global_iterate_families (fclaw_global_t * glob,
                                    fclaw_global_callback_t gpcb, void *user);

#endif /* FCLAW_GLOBAL_H */
