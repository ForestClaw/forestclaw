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
 * Dimension-independent wrapper of a forestclaw patch.
 */

#ifndef FCLAW_PATCH_H
#define FCLAW_PATCH_H

/*
 * Domain-independent header file should not include domain-specific headers.
 * The corresponding source file include the 2d and 3d domain-specific headers.
 */
#include <fclaw_base.h>

typedef struct fclaw_patch fclaw_patch_t;

struct fclaw_patch_bounds_2d
{
    double xlower, xupper;
    double ylower, yupper;
};

struct fclaw_patch_bounds_3d
{
    double xlower, xupper;
    double ylower, yupper;
    double zlower, zupper;
};

/** 
 * @brief The metadata structure for a forest leaf, which is a forestclaw patch.
 * The patch may be either a process-local patch or a ghost patch.
 */
struct fclaw_patch
{
    int dim;
    struct fclaw_patch_bounds_2d* d2;
    struct fclaw_patch_bounds_3d* d3;
    int level;                  /**< 0 is root, increases if refined */
    int target_level;           /**< level desired after adaptation */
    int flags;                  /**< flags that encode tree information */
    /** Union, If this is a local patch, it points to the next patch, otherwise it gives
     * the bock number of this patch */
    union
    {
        fclaw_patch_t *next;  /**< local: next patch same level same block */
        int blockno;            /**< off-proc: this patch's block number */
    }
    u;
    void *user;                 /**< User Pointer */
};
#endif /* !FCLAW_PATCH_H */
