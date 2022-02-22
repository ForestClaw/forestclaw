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
/** 
 * @file
 * Dimension-independent version of domain, block, patch
 */

#ifndef FORESTCLAW_H
#define FORESTCLAW_H

#include <fclaw_base.h>

typedef struct fclaw_patch
{
    union
    {
        /* avoid including dimension-specific files */
        struct fclaw2d_patch_data *pd2;
        struct fclaw3d_patch_data *pd3;
    }
    pd;

    /* Those below are currently also in pd2, pd3.
       Choose exactly one good place for them to be in */

  /** User defined patch structure */
    void *user_patch;

  /** Additional user data */
    void *user_data;
}
fclaw_patch_t;

typedef struct fclaw_domain
{
    /* store dimension-specific domain structures */
    int dim;
    union
    {
        struct
        {
            /* avoid including dimension-specific files */
            int dmagic2;
            struct fclaw2d_domain *domain2;
        }
        d2;
        struct
        {
            int dmagic3;
            struct fclaw3d_domain *domain3;
        }
        d3;
    }
    d;

    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;
}
fclaw_domain_t;

/** Dimension-independent callback prototype for the patch iterators.
 * We iterate over local patches only.
 * \param [in] domain	Dimension-independent domain structure.
 * \param [in] patch	The local patch currently processed by the iterator.
 * \param [in] blockno  Block number of processed patch.
 * \param [in] patchno  Patch number within block of processed patch.
 * \param [in,out] user	Data that was passed into the iterator functions.
 */
typedef void (*fclaw_patch_callback_t)
    (fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user);

/* Verify that dimension and domain object are set correctly.
 * \param [in] domain           Dimension-independent domain.
 * \return                      True if dimension and magic valid and a
 *                              non-NULL dimension-specific domain assigned.
 */
int fclaw_domain_is_valid (fclaw_domain_t * domain);

/** Destruct a dimension-specific domain and its patch data.
 * \param [in,out] domain       Initialized, valid domain structure.
 *                              On output, the pointer is invalidated.
 * \param [in] dele             This callback is pointed to an existing
 *                              fclaw_patch whose data it is supposed to
 *                              delete.  This includes deleting
 *                              the patch user pointer, which the
 *                              callback can do using the third
 *                              parameter to this function, \a user.
 *                              Just do not delete the patch itself!
 * \param [in,out] user         Pointer passed through to \a init.
 */
void fclaw_domain_destroy (fclaw_domain_t * domain,
                           fclaw_patch_callback_t dele, void *user);

/** Iterate over all local patches on a given level.
 * \param [in] domain	Dimension-independent domain structure.
 * \param [in] level	Level to iterate.  Ignore patches of other levels.
 * \param [in] iter     Function called for each patch of matching level.
 * \param [in,out] user	Data is passed to the \a iter callback.
 */
void fclaw_domain_iterate_level (fclaw_domain_t * domain, int level,
                                 fclaw_patch_callback_t iter, void *user);

#endif /* !FORESTCLAW_H */
