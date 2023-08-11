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

#ifndef FCLAW3D_MAP_H
#define FCLAW3D_MAP_H

#include <fclaw_base.h>
#if 0
/* put whatever is needed into that file */
#include <fclaw3d_map_query_defs.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

struct fclaw3d_map_context;
typedef struct fclaw3d_map_context fclaw_map_context_t;

/** Destructor for a fclaw3d_map_context.
*/
typedef void (*fclaw3d_map_destroy_t) (fclaw_map_context_t * cont);

/** Mapping context that is interpreted by its query and c2m members.
* The callbacks are free to define the meaning of the user_* fields.
*/
struct fclaw3d_map_context
{
  fclaw3d_map_destroy_t destroy;
    void *user_data;
};

/** Deallocate a mapping context.
 * If the \a destroy member is not NULL, it is called on the context.
 * Otherwise, this function calls FCLAW_FREE (cont).
 * \param [in] cont     Mapping context where the \a destroy member is either
 *                      NULL or a valid function that is then called.
 */
void fclaw3d_map_destroy (fclaw_map_context_t * cont);

fclaw_map_context_t* fclaw3d_map_new_nomap (void);

#ifdef __cplusplus
}
#endif

#endif /* !FCLAW3D_MAP_H */
