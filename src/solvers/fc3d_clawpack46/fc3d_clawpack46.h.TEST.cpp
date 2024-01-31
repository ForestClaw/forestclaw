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
#include <fc3d_clawpack46.h>
#include <fc3d_clawpack46_options.h>
#include <fclaw_forestclaw.h>
#include <fclaw_convenience.h>
#include <fclaw_clawpatch_options.h>
#include <test.hpp>

TEST_CASE("fc3d_clawpack46_solver_initialize stores two seperate vtables in two seperate globs")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	fclaw_clawpatch_options_t* opts1 = fclaw_clawpatch_options_new(3);
	fclaw_clawpatch_options_t* opts2 = fclaw_clawpatch_options_new(3);

	/* create some empty options structures */
	fclaw_clawpatch_options_store(glob1, opts1);
	fc3d_clawpack46_options_store(glob1, FCLAW_ALLOC_ZERO(fc3d_clawpack46_options_t,1));

	fclaw_clawpatch_options_store(glob2, opts2);
	fc3d_clawpack46_options_store(glob2, FCLAW_ALLOC_ZERO(fc3d_clawpack46_options_t,1));

	fclaw_vtables_initialize(glob1);
	fc3d_clawpack46_solver_initialize(glob1);

	fclaw_vtables_initialize(glob2);
	fc3d_clawpack46_solver_initialize(glob2);

	CHECK_NE(fc3d_clawpack46_vt(glob1), fc3d_clawpack46_vt(glob2));

	fclaw_clawpatch_options_destroy(opts1);
	fclaw_clawpatch_options_destroy(opts2);
	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fc3d_clawpack46_solver_initialize sets is_set flag")
{
	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	fclaw_global_t* glob = fclaw_global_new();
	fclaw_global_store_domain(glob, domain);

	fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(3);

	/* create some empty options structures */
	fclaw_clawpatch_options_store(glob, opts);
	fc3d_clawpack46_options_store(glob, FCLAW_ALLOC_ZERO(fc3d_clawpack46_options_t,1));

	fclaw_vtables_initialize(glob);
	fc3d_clawpack46_solver_initialize(glob);


	CHECK_UNARY(fc3d_clawpack46_vt(glob)->is_set);

	fclaw_clawpatch_options_destroy(opts);
	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fc3d_clawpack46_vt fails if not intialized")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	CHECK_SC_ABORTED(fc3d_clawpack46_vt(glob1));

	CHECK_SC_ABORTED(fc3d_clawpack46_vt(glob2));

	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

TEST_CASE("fc3d_clawpack46_vtable_initialize fails if called twice on a glob")
{
	fclaw_global_t* glob1 = fclaw_global_new();
	fclaw_global_t* glob2 = fclaw_global_new();

	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	fclaw_global_store_domain(glob1, domain);
	fclaw_global_store_domain(glob2, domain);

	fclaw_clawpatch_options_t* opts1 = fclaw_clawpatch_options_new(3);
	fclaw_clawpatch_options_t* opts2 = fclaw_clawpatch_options_new(3);

	/* create some empty options structures */
	fclaw_clawpatch_options_store(glob1, opts1);
	fc3d_clawpack46_options_store(glob1, FCLAW_ALLOC_ZERO(fc3d_clawpack46_options_t,1));

	fclaw_clawpatch_options_store(glob2, opts2);
	fc3d_clawpack46_options_store(glob2, FCLAW_ALLOC_ZERO(fc3d_clawpack46_options_t,1));

	fclaw_vtables_initialize(glob1);
	fc3d_clawpack46_solver_initialize(glob1);
	CHECK_SC_ABORTED(fc3d_clawpack46_solver_initialize(glob1));

	fclaw_vtables_initialize(glob2);
	fc3d_clawpack46_solver_initialize(glob2);
	CHECK_SC_ABORTED(fc3d_clawpack46_solver_initialize(glob2));

	fclaw_clawpatch_options_destroy(opts1);
	fclaw_clawpatch_options_destroy(opts2);
	fclaw_domain_destroy(domain);
	fclaw_global_destroy(glob1);
	fclaw_global_destroy(glob2);
}

#endif