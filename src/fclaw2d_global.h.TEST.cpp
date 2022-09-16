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
#include <test.hpp>

TEST_CASE("fclaw2d_global_set_global")
{
	fclaw2d_global_t* glob = (fclaw2d_global_t*)123;
	fclaw2d_global_set_global(glob);
	CHECK_EQ(fclaw2d_global_get_global(), glob);
	fclaw2d_global_unset_global();
}

TEST_CASE("fclaw2d_global_unset_global")
{
	fclaw2d_global_t* glob = (fclaw2d_global_t*)123;
	fclaw2d_global_set_global(glob);
	fclaw2d_global_unset_global();
#ifdef FCLAW_ENABLE_DEBUG
	CHECK_SC_ABORTED(fclaw2d_global_get_global());
#else
	CHECK_EQ(fclaw2d_global_get_global(), nullptr);
#endif
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw2d_global_set_global twice fails")
{
	fclaw2d_global_t* glob = (fclaw2d_global_t*)123;
	fclaw2d_global_set_global(glob);
	CHECK_SC_ABORTED(fclaw2d_global_set_global(glob));
	fclaw2d_global_unset_global();
}

TEST_CASE("fclaw2d_global_unset_global assert fails when NULL")
{
	CHECK_SC_ABORTED(fclaw2d_global_unset_global());
}

TEST_CASE("fclaw2d_global_get_global assert fails when NULL")
{
	CHECK_SC_ABORTED(fclaw2d_global_get_global());
}

#endif
