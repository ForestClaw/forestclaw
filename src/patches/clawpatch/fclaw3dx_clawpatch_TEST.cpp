/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch46_fort.h>
#include <fclaw3dx_clawpatch_output_ascii.h>
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include <test/catch.hpp>
#include <fstream>

TEST_CASE("fclaw3dx_clawpatch_vtable_initialize","[fclaw3dx][clawpatch]")
{
	fclaw3dx_clawpatch_vtable_initialize(4);

	fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt();

	CHECK(clawpatch_vt->set_user_data             == NULL);

	//ghost filling
	CHECK(clawpatch_vt->fort_copy_face              == &FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE);
	CHECK(clawpatch_vt->fort_average_face           == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE);
	CHECK(clawpatch_vt->fort_interpolate_face       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE);
	CHECK(clawpatch_vt->fort_copy_corner            == &FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER);
	CHECK(clawpatch_vt->fort_average_corner         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER);
	CHECK(clawpatch_vt->fort_interpolate_corner     == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER);

	//regridding
	CHECK(clawpatch_vt->fort_tag4refinement         == &FCLAW3DX_CLAWPATCH46_FORT_TAG4REFINEMENT);
	CHECK(clawpatch_vt->fort_tag4coarsening         == &FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING);
	CHECK(clawpatch_vt->fort_user_exceeds_threshold == NULL);
	CHECK(clawpatch_vt->fort_interpolate2fine       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE);
	CHECK(clawpatch_vt->fort_average2coarse         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE);

	//ascii output
	CHECK(clawpatch_vt->time_header_ascii           == &fclaw3dx_clawpatch_time_header_ascii);
	CHECK(clawpatch_vt->fort_header_ascii           == &FCLAW3DX_CLAWPATCH46_FORT_HEADER_ASCII);
	CHECK(clawpatch_vt->cb_output_ascii             == &fclaw3dx_clawpatch_output_ascii_cb);
	CHECK(clawpatch_vt->fort_output_ascii           == &FCLAW3DX_CLAWPATCH46_FORT_OUTPUT_ASCII);

	//time interpolation
	CHECK(clawpatch_vt->fort_timeinterp             == &FCLAW3DX_CLAWPATCH46_FORT_TIMEINTERP);

	//ghot packing
	CHECK(clawpatch_vt->fort_local_ghost_pack       == &FCLAW3DX_CLAWPATCH46_FORT_LOCAL_GHOST_PACK);
	CHECK(clawpatch_vt->local_ghost_pack_aux        == NULL);

	//diagnostics
	CHECK(clawpatch_vt->conservation_check          == NULL);
	CHECK(clawpatch_vt->compute_error               == NULL);
	CHECK(clawpatch_vt->fort_compute_patch_error    == NULL);
	CHECK(clawpatch_vt->fort_conservation_check     == NULL);
	CHECK(clawpatch_vt->fort_compute_error_norm     == NULL);
	CHECK(clawpatch_vt->fort_compute_patch_area     == NULL);

	CHECK(clawpatch_vt->is_set                      == 1);

	fclaw2d_patch_vtable_t * patch_vt = fclaw2d_patch_vt();
	//create delete build
	//TODO document patch_vt and expose these as part to public api
	CHECK(patch_vt->patch_new                      != NULL);
	CHECK(patch_vt->patch_delete                   != NULL);
	CHECK(patch_vt->build                          != NULL);
	CHECK(patch_vt->build_from_fine                != NULL);
	CHECK(patch_vt->setup                          == NULL);
}