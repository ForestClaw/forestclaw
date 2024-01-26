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
#include <fclaw2d_clawpatch_options.h>
#include <fclaw3dx_clawpatch_options.h>
#include <fclaw_packing.h>
#include <test.hpp>

TEST_CASE("fclaw2d_clawpatch_options can store options in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw2d_clawpatch_options_t* opts1 = FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1);
	fclaw2d_clawpatch_options_t* opts2 = FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1);

	fclaw2d_clawpatch_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw3dx_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1));
	fclaw2d_clawpatch_options_store(glob2, opts2);

	CHECK_EQ(fclaw2d_clawpatch_get_options(glob1), opts1);
	CHECK_EQ(fclaw2d_clawpatch_get_options(glob2), opts2);

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw3dx_clawpatch_options can store options in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw3dx_clawpatch_options_t* opts1 = FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1);
	fclaw3dx_clawpatch_options_t* opts2 = FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1);

	fclaw3dx_clawpatch_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw2d_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1));
	fclaw3dx_clawpatch_options_store(glob2, opts2);

	CHECK_EQ(fclaw3dx_clawpatch_get_options(glob1), opts1);
	CHECK_EQ(fclaw3dx_clawpatch_get_options(glob2), opts2);

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw2d_clawpatch_get_options fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	CHECK_SC_ABORTED(fclaw2d_clawpatch_get_options(glob1));

	CHECK_SC_ABORTED(fclaw2d_clawpatch_get_options(glob2));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw3dx_clawpatch_get_options fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	CHECK_SC_ABORTED(fclaw3dx_clawpatch_get_options(glob1));

	CHECK_SC_ABORTED(fclaw3dx_clawpatch_get_options(glob2));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw2d_clawpatch_options_store fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw2d_clawpatch_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1));
	CHECK_SC_ABORTED(fclaw2d_clawpatch_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1)));

	fclaw2d_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1));
	CHECK_SC_ABORTED(fclaw2d_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1)));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw3dx_clawpatch_options_store fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw3dx_clawpatch_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1));
	CHECK_SC_ABORTED(fclaw3dx_clawpatch_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1)));

	fclaw3dx_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1));
	CHECK_SC_ABORTED(fclaw3dx_clawpatch_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1)));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fclaw3dx_clawpatch_get_options fails is fclaw2d_clawpatch_options is set")
{
	fclaw_global_t* glob = fclaw_global_new();

	fclaw2d_clawpatch_options_store(glob, FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1));

	CHECK_SC_ABORTED(fclaw3dx_clawpatch_get_options(glob));

	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw2d_clawpatch_options packing/unpacking")
{
	fclaw2d_clawpatch_options_t* opts = FCLAW_ALLOC_ZERO(fclaw2d_clawpatch_options_t,1);
	opts->mx = 5;
	opts->my = 6;
	opts->maux = 4;
	opts->mbc = 3;
	opts->meqn = 32;
	opts->rhs_fields = 39;
	opts->refinement_criteria = 1;
	opts->interp_stencil_width = 3;
	opts->ghost_patch_pack_aux = 7;
	opts->save_aux = 1;
	opts->is_registered = 1;

	const fclaw_packing_vtable_t* vt = fclaw2d_clawpatch_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fclaw2d_clawpatch_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(output_opts->mx,opts->mx);
	CHECK_EQ(output_opts->my,opts->my);
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
	FCLAW_FREE(opts);
}

TEST_CASE("fclaw3dx_clawpatch_options packing/unpacking")
{
	fclaw3dx_clawpatch_options_t* opts = FCLAW_ALLOC_ZERO(fclaw3dx_clawpatch_options_t,1);
	opts->mx = 5;
	opts->my = 6;
	opts->mz = 3;
	opts->maux = 4;
	opts->mbc = 3;
	opts->meqn = 32;
	opts->rhs_fields = 39;
	opts->refinement_criteria = 1;
	opts->interp_stencil_width = 3;
	opts->ghost_patch_pack_aux = 7;
	opts->save_aux = 1;
	opts->is_registered = 1;

	const fclaw_packing_vtable_t* vt = fclaw3dx_clawpatch_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fclaw3dx_clawpatch_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(output_opts->mx,opts->mx);
	CHECK_EQ(output_opts->my,opts->my);
	CHECK_EQ(output_opts->mz,opts->mz);
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
	FCLAW_FREE(opts);
}
#endif