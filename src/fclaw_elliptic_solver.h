/*
  Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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


#ifndef FCLAW2D_ELLIPTIC_SOLVER_H
#define FCLAW2D_ELLIPTIC_SOLVER_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;

typedef void (*fclaw_elliptic_setup_t)(struct fclaw_global *glob);

typedef void (*fclaw_elliptic_rhs_t)(struct fclaw_global *glob);

typedef void (*fclaw_elliptic_solve_t)(struct fclaw_global *glob);

typedef void (*fclaw_elliptic_physical_bc_t)(struct fclaw_global *glob);



typedef struct fclaw_elliptic_vtable
{
    fclaw_elliptic_setup_t       setup;
    fclaw_elliptic_rhs_t         rhs;
    fclaw_elliptic_solve_t       solve;
    fclaw_elliptic_physical_bc_t apply_bc;

    int is_set;

} fclaw_elliptic_vtable_t;

void fclaw_elliptic_vtable_initialize(struct fclaw_global* glob);

void fclaw_elliptic_solve(struct fclaw_global *glob);

fclaw_elliptic_vtable_t* fclaw_elliptic_vt(struct fclaw_global* glob);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_CLAWPACH46_H */

