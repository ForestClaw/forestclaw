/*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_clawpatch_vtk_vtable.h>
#include <fclaw_global.h>

static
fclaw_clawpatch_vtk_vtable_entry_t* vtk_vt_entry_new(const char* name,
                                                     fclaw_vtk_entry_type_t type,
                                                     size_t size_per_patch,
                                                     fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		FCLAW_ALLOC(fclaw_clawpatch_vtk_vtable_entry_t, 1);

	entry->name = name;
	entry->type = type;
	entry->size_per_patch = size_per_patch;
	entry->callback = callback;

	return entry;
}

static
void vtk_add_field_entry(fclaw_clawpatch_vtk_vtable_t* vtk_vt,
						 const char* name,
						 fclaw_vtk_entry_type_t type,
						 size_t size_per_patch,
						 fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		vtk_vt_entry_new(name, type, size_per_patch, callback);
	sc_list_append(vtk_vt->field_entries, entry);
}

static
void vtk_add_point_entry(fclaw_clawpatch_vtk_vtable_t* vtk_vt,
						 const char* name,
						 fclaw_vtk_entry_type_t type,
						 size_t size_per_patch,
						 fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		vtk_vt_entry_new(name, type, size_per_patch, callback);
	sc_list_append(vtk_vt->point_entries, entry);
}

static
void vtk_add_cell_entry(fclaw_clawpatch_vtk_vtable_t* vtk_vt,
						const char* name,
						fclaw_vtk_entry_type_t type,
						size_t size_per_patch,
						fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		vtk_vt_entry_new(name, type, size_per_patch, callback);
	sc_list_append(vtk_vt->cell_entries, entry);
}

static
fclaw_clawpatch_vtk_vtable_t* vtk_vt_new(void)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*)
		FCLAW_ALLOC_ZERO(fclaw_clawpatch_vtk_vtable_t, 1);

	vtk_vt->field_entries = sc_list_new(NULL);
	vtk_vt->point_entries = sc_list_new(NULL);
	vtk_vt->cell_entries = sc_list_new(NULL);

	return vtk_vt;
}

static
void vtk_vt_destroy(void* vt)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*) vt;
	fclaw_clawpatch_vtk_vtable_entry_t* entry;

	if (vtk_vt->field_entries != NULL)
	{
		while (vtk_vt->field_entries->elem_count > 0)
		{
			entry = (fclaw_clawpatch_vtk_vtable_entry_t*)
				sc_list_pop(vtk_vt->field_entries);
			FCLAW_FREE(entry);
		}
		sc_list_destroy(vtk_vt->field_entries);
	}
	if (vtk_vt->point_entries != NULL)
	{
		while (vtk_vt->point_entries->elem_count > 0)
		{
			entry = (fclaw_clawpatch_vtk_vtable_entry_t*)
				sc_list_pop(vtk_vt->point_entries);
			FCLAW_FREE(entry);
		}
		sc_list_destroy(vtk_vt->point_entries);
	}
	if (vtk_vt->cell_entries != NULL)
	{
		while (vtk_vt->cell_entries->elem_count > 0)
		{
			entry = (fclaw_clawpatch_vtk_vtable_entry_t*)
				sc_list_pop(vtk_vt->cell_entries);
			FCLAW_FREE(entry);
		}
		sc_list_destroy(vtk_vt->cell_entries);
	}

	FCLAW_FREE(vt);
}

fclaw_clawpatch_vtk_vtable_t*
fclaw_clawpatch_vtk_vtable(struct fclaw_global* glob)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*)
		fclaw_global_get_vtable((fclaw_global_t*) glob,
								"fclaw_clawpatch_vtk_vtable");
	FCLAW_ASSERT(vtk_vt != NULL);
	return vtk_vt;
}

void fclaw_clawpatch_vtk_vtable_initialize(struct fclaw_global* glob)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt = vtk_vt_new();

	fclaw_global_vtable_store((fclaw_global_t*) glob,
							  "fclaw_clawpatch_vtk_vtable",
							  vtk_vt,
							  vtk_vt_destroy);
}


