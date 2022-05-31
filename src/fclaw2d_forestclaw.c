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

#include <fclaw2d_forestclaw.h>

#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_diagnostics.h>
#include <fclaw_gauges.h>

#include <fclaw2d_elliptic_solver.h>

void fclaw2d_vtables_initialize(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_initialize(glob);
    fclaw2d_patch_vtable_initialize(glob);
    fclaw2d_diagnostics_vtable_initialize();
    fclaw2d_elliptic_vtable_initialize();

    fclaw_gauges_vtable_initialize();    
}

void fclaw2d_problem_setup(fclaw2d_global_t *glob)
{
	fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt(glob);
	
    /* User defined problem setup */
    if (fclaw_vt->problem_setup != NULL)
    {
        fclaw_vt->problem_setup(glob);
    }
}
