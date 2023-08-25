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

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_packing.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_options.h>
#include <test.hpp>

TEST_CASE("fclaw_clawpatch_options_new 2d")
{
	fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(2);

	CHECK_EQ(opts->dim,     2);
	CHECK_NE(opts->d2,      nullptr);
	CHECK_EQ(opts->d3,      nullptr);
	CHECK_EQ(opts->d2->mx,  0);

	fclaw_clawpatch_options_destroy(opts);
}

TEST_CASE("fclaw_clawpatch_options_new 3d")
{
	fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(3);

	CHECK_EQ(opts->dim,     3);
	CHECK_EQ(opts->d2,      nullptr);
	CHECK_NE(opts->d3,      nullptr);
	CHECK_EQ(opts->d3->mx,  0);

	fclaw_clawpatch_options_destroy(opts);
}

TEST_CASE("fclaw_clawpatch_options can store options in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw_clawpatch_options_t* opts1 = fclaw_clawpatch_options_new(2);
	fclaw_clawpatch_options_t* opts2 = fclaw_clawpatch_options_new(3);

	fclaw_clawpatch_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw_clawpatch_options_store(glob2, opts2);

	CHECK_EQ(fclaw_clawpatch_get_options(glob1), opts1);
	CHECK_EQ(fclaw_clawpatch_get_options(glob2), opts2);

	fclaw_clawpatch_options_destroy(opts1);
	fclaw_clawpatch_options_destroy(opts2);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw_clawpatch_get_options fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	CHECK_SC_ABORTED(fclaw_clawpatch_get_options(glob1));

	CHECK_SC_ABORTED(fclaw_clawpatch_get_options(glob2));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw_clawpatch_options_store fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw_clawpatch_options_t* opts1 = fclaw_clawpatch_options_new(2);
	fclaw_clawpatch_options_t* opts2 = fclaw_clawpatch_options_new(3);

	fclaw_clawpatch_options_store(glob1, opts1);
	CHECK_SC_ABORTED(fclaw_clawpatch_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw_clawpatch_options_t,1)));

	fclaw_clawpatch_options_store(glob2, opts2);
	CHECK_SC_ABORTED(fclaw_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw_clawpatch_options_t,1)));

	fclaw_clawpatch_options_destroy(opts1);
	fclaw_clawpatch_options_destroy(opts2);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

#endif
TEST_CASE("2d fclaw_clawpatch_options packing/unpacking")
{
	fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(2);
	opts->d2->mx = 5;
	opts->d2->my = 6;
	opts->maux = 4;
	opts->mbc = 3;
	opts->meqn = 32;
	opts->rhs_fields = 39;
	opts->refinement_criteria = 1;
	opts->interp_stencil_width = 3;
	opts->ghost_patch_pack_aux = 7;
	opts->save_aux = 1;
	opts->is_registered = 1;

	const fclaw_packing_vtable_t* vt = fclaw_clawpatch_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fclaw_clawpatch_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(output_opts->dim,2);
	CHECK_EQ(output_opts->d2->mx,opts->d2->mx);
	CHECK_EQ(output_opts->d2->my,opts->d2->my);
	CHECK_EQ(output_opts->maux,opts->maux);
	CHECK_EQ(output_opts->mbc,opts->mbc);
	CHECK_EQ(output_opts->meqn,opts->meqn);
	CHECK_EQ(output_opts->rhs_fields,opts->rhs_fields);
	CHECK_EQ(output_opts->refinement_criteria,opts->refinement_criteria);
	CHECK_EQ(output_opts->interp_stencil_width,opts->interp_stencil_width);
	CHECK_EQ(output_opts->ghost_patch_pack_aux,opts->ghost_patch_pack_aux);
	CHECK_EQ(output_opts->save_aux,opts->save_aux);
	CHECK_EQ(output_opts->is_registered,opts->is_registered);

	vt->destroy(output_opts);
	FCLAW_FREE(opts->d2);
	FCLAW_FREE(opts);
}

TEST_CASE("3d fclaw_clawpatch_options packing/unpacking")
{
	fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(3);
	opts->d3->mx = 5;
	opts->d3->my = 6;
	opts->d3->mz = 3;
	opts->maux = 4;
	opts->mbc = 3;
	opts->meqn = 32;
	opts->rhs_fields = 39;
	opts->refinement_criteria = 1;
	opts->interp_stencil_width = 3;
	opts->ghost_patch_pack_aux = 7;
	opts->save_aux = 1;
	opts->is_registered = 1;

	const fclaw_packing_vtable_t* vt = fclaw_clawpatch_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fclaw_clawpatch_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(output_opts->dim,3);
	CHECK_EQ(output_opts->d3->mx,opts->d3->mx);
	CHECK_EQ(output_opts->d3->my,opts->d3->my);
	CHECK_EQ(output_opts->d3->mz,opts->d3->mz);
	CHECK_EQ(output_opts->maux,opts->maux);
	CHECK_EQ(output_opts->mbc,opts->mbc);
	CHECK_EQ(output_opts->meqn,opts->meqn);
	CHECK_EQ(output_opts->rhs_fields,opts->rhs_fields);
	CHECK_EQ(output_opts->refinement_criteria,opts->refinement_criteria);
	CHECK_EQ(output_opts->interp_stencil_width,opts->interp_stencil_width);
	CHECK_EQ(output_opts->ghost_patch_pack_aux,opts->ghost_patch_pack_aux);
	CHECK_EQ(output_opts->save_aux,opts->save_aux);
	CHECK_EQ(output_opts->is_registered,opts->is_registered);

	vt->destroy(output_opts);
	FCLAW_FREE(opts->d3);
	FCLAW_FREE(opts);
}