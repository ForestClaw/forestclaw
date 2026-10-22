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

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_vtk_vtable.h>
#include <fclaw_global.h>
#include <test.hpp>
#include <string.h>

static fclaw_clawpatch_vtk_vtable_entry_t*
find_entry(sc_list_t* entries, const char* name)
{
    for (sc_link_t* link = entries->first; link != NULL; link = link->next)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) link->data;
        if (strcmp(entry->name, name) == 0)
        {
            return entry;
        }
    }
    return NULL;
}

static fclaw_clawpatch_vtk_vtable_t*
initialize_vtk_vtable(fclaw_global_t* glob, fclaw_clawpatch_options_t* options)
{
    fclaw_clawpatch_options_store(glob, options);
    fclaw_clawpatch_vtk_vtable_initialize(glob);
    return fclaw_clawpatch_vtk_vtable(glob);
}

TEST_CASE("VTK vtable contains standard fields and cell data")
{
    fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
    fclaw_clawpatch_options_t* options = fclaw_clawpatch_options_new(2);
    options->meqn = 2;
    options->rhs_fields = 2;
    options->maux = 2;
    options->vtk_aux_out = FCLAW_ALLOC(int, 2);
    options->vtk_aux_out[0] = 1;
    options->vtk_aux_out[1] = 0;

    fclaw_clawpatch_vtk_vtable_t* vtk_vt = initialize_vtk_vtable(glob, options);

    REQUIRE_EQ(vtk_vt->field_entries->elem_count, 3);
    CHECK(find_entry(vtk_vt->field_entries, "levels") != NULL);
    CHECK(find_entry(vtk_vt->field_entries, "patch_starts") != NULL);
    CHECK(find_entry(vtk_vt->field_entries, "patch_spacings") != NULL);

    CHECK_EQ(vtk_vt->position_entry.type, FCLAW_VTK_FLOAT64);
    CHECK_EQ(vtk_vt->position_entry.number_of_components, 3);
    CHECK(vtk_vt->position_entry.elements_in_patch != NULL);
    CHECK(vtk_vt->position_entry.callback != NULL);
    CHECK_EQ(vtk_vt->connectivity_entry.type, FCLAW_VTK_INT64);
    CHECK_EQ(vtk_vt->offsets_entry.type, FCLAW_VTK_INT64);
    CHECK_EQ(vtk_vt->types_entry.type, FCLAW_VTK_UINT8);

    REQUIRE_EQ(vtk_vt->celldata_entries->elem_count, 8);
    CHECK(find_entry(vtk_vt->celldata_entries, "mpirank") != NULL);
    CHECK(find_entry(vtk_vt->celldata_entries, "blockno") != NULL);
    CHECK(find_entry(vtk_vt->celldata_entries, "patchno") != NULL);

    fclaw_clawpatch_vtk_vtable_entry_t* meqn =
        find_entry(vtk_vt->celldata_entries, "meqn");
    REQUIRE(meqn != NULL);
    CHECK_EQ(meqn->number_of_components, 2);
    CHECK(meqn->elements_in_patch != NULL);
    CHECK(meqn->callback != NULL);

    fclaw_clawpatch_vtk_vtable_entry_t* aux =
        find_entry(vtk_vt->celldata_entries, "aux");
    REQUIRE(aux != NULL);
    CHECK_EQ(aux->number_of_components, 1);

    CHECK(find_entry(vtk_vt->celldata_entries, "rhs") != NULL);
    CHECK(find_entry(vtk_vt->celldata_entries, "soln") != NULL);
    CHECK(find_entry(vtk_vt->celldata_entries, "error") != NULL);
    CHECK(find_entry(vtk_vt->celldata_entries, "vtkGhostType") == NULL);

    fclaw_global_destroy(glob);
    fclaw_clawpatch_options_destroy(options);
}

TEST_CASE("VTK vtable includes ghost-cell metadata when requested")
{
    fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
    fclaw_clawpatch_options_t* options = fclaw_clawpatch_options_new(2);
    options->vtk_ghost_out = 1;

    fclaw_clawpatch_vtk_vtable_t* vtk_vt = initialize_vtk_vtable(glob, options);
    fclaw_clawpatch_vtk_vtable_entry_t* ghost =
        find_entry(vtk_vt->celldata_entries, "vtkGhostType");

    REQUIRE(ghost != NULL);
    CHECK_EQ(ghost->type, FCLAW_VTK_UINT8);
    CHECK_EQ(ghost->number_of_components, 1);
    CHECK(ghost->elements_in_patch != NULL);
    CHECK(ghost->callback != NULL);
    CHECK_EQ(vtk_vt->celldata_entries->elem_count, 5);

    fclaw_global_destroy(glob);
    fclaw_clawpatch_options_destroy(options);
}
