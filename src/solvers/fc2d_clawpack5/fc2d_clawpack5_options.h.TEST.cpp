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
#include <fc2d_clawpack5_options.h>
#include <test.hpp>
#include <fclaw_packing.h>

TEST_CASE("fc2d_clawpack5_options can store options in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fc2d_clawpack5_options_t* opts1 = FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1);
	fc2d_clawpack5_options_t* opts2 = FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1);

	fc2d_clawpack5_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw_options_t,1));
	fc2d_clawpack5_options_store(glob2, opts2);

	CHECK_EQ(fc2d_clawpack5_get_options(glob1), opts1);
	CHECK_EQ(fc2d_clawpack5_get_options(glob2), opts2);

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fc2d_clawpack5_get_options fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	CHECK_SC_ABORTED(fc2d_clawpack5_get_options(glob1));

	CHECK_SC_ABORTED(fc2d_clawpack5_get_options(glob2));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fc2d_clawpack5_options_store fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fc2d_clawpack5_options_store(glob1, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1));
	CHECK_SC_ABORTED(fc2d_clawpack5_options_store(glob1, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1)));

	fc2d_clawpack5_options_store(glob2, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1));
	CHECK_SC_ABORTED(fc2d_clawpack5_options_store(glob2, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1)));

	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fc2d_clawpack5_options packing/unpacking")
{
	fc2d_clawpack5_options_t* opts = FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1);
	opts->mwaves = 5;
	opts->order_string = "[2 2]";
	int order[] = {2,2};
	opts->order = order;
	int mthlim[] = {5,4,3,2,1};
	opts->mthlim = mthlim;
	opts->mthlim_string = "[5 4 3 2 1]";
	int mthbc[] = {0, 1, 1, 1};
	opts->mthbc = mthbc;
	opts->mthbc_string = "[0 1 1 1]";

	int method[] = {0,4,3,2,5,6,99};
	for(size_t i = 0; i < 7; ++i)
		opts->method[i] = method[i];

	opts->mcapa = 3;
	opts->src_term = 2;
	opts->use_fwaves = 1;

	opts->ascii_out = 2;
	opts->vtk_out = 3;

	opts->is_registered = 1;

	opts->is_unpacked = 0;

	const fclaw_packing_vtable_t* vt = fc2d_clawpack5_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fc2d_clawpack5_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(output_opts->mwaves, opts->mwaves);

	CHECK_NE(output_opts->order_string, opts->order_string);
	CHECK_UNARY(!strcmp(output_opts->order_string, opts->order_string));

	CHECK_EQ(output_opts->order[0], opts->order[0]);
	CHECK_EQ(output_opts->order[1], opts->order[1]);

	for(size_t i = 0; i < 5; ++i)
		CHECK_EQ(output_opts->mthlim[i], opts->mthlim[i]);
	
	CHECK_NE(output_opts->mthlim_string, opts->mthlim_string);
	CHECK_UNARY(!strcmp(output_opts->mthlim_string, opts->mthlim_string));

	for(size_t i = 0; i < 4; ++i)
		CHECK_EQ(output_opts->mthbc[i], opts->mthbc[i]);

	CHECK_NE(output_opts->mthbc_string, opts->mthbc_string);
	CHECK_UNARY(!strcmp(output_opts->mthbc_string, opts->mthbc_string));

	for(size_t i = 0; i < 7; ++i)
		CHECK_EQ(output_opts->method[i], opts->method[i]);

	CHECK_EQ(output_opts->mcapa, opts->mcapa);
	CHECK_EQ(output_opts->src_term, opts->src_term);
	CHECK_EQ(output_opts->use_fwaves, opts->use_fwaves);
	CHECK_EQ(output_opts->ascii_out, opts->ascii_out);
	CHECK_EQ(output_opts->vtk_out, opts->vtk_out);
	CHECK_EQ(output_opts->is_registered, opts->is_registered);
	CHECK_UNARY(output_opts->is_unpacked);

	vt->destroy(output_opts);
	FCLAW_FREE(opts);
}
#endif