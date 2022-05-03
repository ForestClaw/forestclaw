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

#include <fclaw3dx_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <test/doctest.h>
#include <fstream>

TEST_CASE("fclaw3dx_clawpatch_set_refinement_criteria, then fclaw3dx_clawpatch_get_refinment_criteria")
{
	for(int r : {1,2,3}){
		fclaw3dx_clawpatch_set_refinement_criteria(r);
		CHECK(r == fclaw3dx_clawpatch_get_refinement_criteria());
	}
}

TEST_CASE("fclaw3dx_clawpatch_options_store, then fclaw3dx_clawpatch_get_options")
{
	fclaw2d_global_t* glob = fclaw2d_global_new();
	fclaw3dx_clawpatch_options_t opts;
	fclaw3dx_clawpatch_options_store(glob, &opts);
	CHECK(fclaw3dx_clawpatch_get_options(glob) == &opts);
	fclaw2d_global_destroy(glob);
}
