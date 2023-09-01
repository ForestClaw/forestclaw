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
#include <fclaw_vtable.h>
#include <fclaw_output.h>
#include <fclaw_global.h>

#include <forestclaw.h>

static
fclaw_vtable_t* vt_new()
{
    return (fclaw_vtable_t*) FCLAW_ALLOC_ZERO (fclaw_vtable_t, 1);
}

static
void vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fclaw_vtable_t* fclaw_vt(fclaw_global_t *glob)
{
	fclaw_vtable_t* vt = (fclaw_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fclaw");
	FCLAW_ASSERT(vt != NULL);
	FCLAW_ASSERT(vt->is_set != 0);
	return vt;
}


void fclaw_after_regrid(fclaw_global_t *glob)
{
    fclaw_vtable_t *vt = fclaw_vt(glob);
    if (vt->after_regrid != NULL)
    {
        vt->after_regrid(glob);
    }
}

/* Initialize any settings that can be set here */
void fclaw_vtable_initialize(fclaw_global_t *glob)
{

    fclaw_vtable_t *vt = vt_new();

    vt->is_set = 1;

	if(fclaw_pointer_map_get(glob->vtables,"fclaw") != NULL)
    {
        fclaw_abortf("fclaw_vtable_initialize : fclaw vtable already initialized\n");
    }
	fclaw_pointer_map_insert(glob->vtables, "fclaw", vt, vt_destroy);
}
