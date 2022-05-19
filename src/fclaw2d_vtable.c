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

#include <fclaw_pointer_map.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_output.h>
#include <fclaw2d_global.h>

#include <forestclaw2d.h>

static
fclaw2d_vtable_t* vt_new()
{
    return (fclaw2d_vtable_t*) FCLAW_ALLOC (fclaw2d_vtable_t, 1);
}

static
void vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fclaw2d_vtable_t* fclaw2d_vt(fclaw2d_global_t *glob)
{
	fclaw2d_vtable_t* vt = (fclaw2d_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fclaw2d");
	FCLAW_ASSERT(vt != NULL);
	FCLAW_ASSERT(vt->is_set != 0);
	return vt;
}


void fclaw2d_after_regrid(fclaw2d_global_t *glob)
{
    fclaw2d_vtable_t *fclaw_vt = fclaw2d_vt(glob);
    if (fclaw_vt->after_regrid != NULL)
    {
        fclaw_vt->after_regrid(glob);
    }
}

/* Initialize any settings that can be set here */
void fclaw2d_vtable_initialize(fclaw2d_global_t *glob)
{

    fclaw2d_vtable_t *vt = vt_new();

    /* ------------------------------------------------------------
      Functions below here depend on q and could be solver specific
      ------------------------------------------------------------- */

    /* These may be redefined by the user */
    /* Problem setup */
    vt->problem_setup              = NULL;

    /* Defaults for regridding */
    vt->after_regrid               = NULL;

    /* Defaults for output */
    vt->output_frame               = NULL;

    vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fclaw2d") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fclaw2d", vt, vt_destroy);
}
