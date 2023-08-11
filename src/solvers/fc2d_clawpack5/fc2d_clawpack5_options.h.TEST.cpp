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

#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fc2d_clawpack5_options.h>
#include <test.hpp>

TEST_CASE("fc2d_clawpack5_options can store options in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fc2d_clawpack5_options_t* opts1 = FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1);
	fc2d_clawpack5_options_t* opts2 = FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1);

	fc2d_clawpack5_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw2d_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw_options_t,1));
	fc2d_clawpack5_options_store(glob2, opts2);

	CHECK_EQ(fc2d_clawpack5_get_options(glob1), opts1);
	CHECK_EQ(fc2d_clawpack5_get_options(glob2), opts2);

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fc2d_clawpack5_get_options fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	CHECK_SC_ABORTED(fc2d_clawpack5_get_options(glob1));

	CHECK_SC_ABORTED(fc2d_clawpack5_get_options(glob2));

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}

TEST_CASE("fc2d_clawpack5_options_store fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fc2d_clawpack5_options_store(glob1, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1));
	CHECK_SC_ABORTED(fc2d_clawpack5_options_store(glob1, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1)));

	fc2d_clawpack5_options_store(glob2, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1));
	CHECK_SC_ABORTED(fc2d_clawpack5_options_store(glob2, FCLAW_ALLOC_ZERO(fc2d_clawpack5_options_t,1)));

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}

#endif