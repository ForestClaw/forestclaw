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

#include "fclaw_timer.h"
#include <fclaw2d_global.h>
#include <test.hpp>

TEST_CASE("fclaw2d_global_serialize"){
	for(int count_amr_advance : {1, 2})
	for(int count_ghost_exchange : {1, 2})
	for(int count_amr_regrid : {1, 2})
	for(int count_amr_new_domain : {1, 2})
	for(int count_single_step : {1, 2})
	for(int count_elliptic_grids : {1, 2})
	for(int count_multiproc_corner : {1, 2})
	for(int count_grids_per_proc : {1, 2})
	for(int count_grids_remote_boundary : {1, 2})
	for(int count_grids_local_boundary : {1, 2})
	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	for(int timer_to_set : {1, 2})
	{
		fclaw2d_global_t* glob1;
    	glob1 = FCLAW_ALLOC_ZERO (fclaw2d_global_t, 1);
		glob1->count_amr_advance            = count_amr_advance;
		glob1->count_ghost_exchange         = count_ghost_exchange;
		glob1->count_amr_regrid             = count_amr_regrid;
		glob1->count_amr_new_domain         = count_amr_new_domain;
		glob1->count_single_step            = count_single_step;
		glob1->count_elliptic_grids         = count_elliptic_grids;
		glob1->count_multiproc_corner       = count_multiproc_corner;
		glob1->count_grids_per_proc         = count_grids_per_proc;
		glob1->count_grids_remote_boundary  = count_grids_remote_boundary;
		glob1->count_grids_local_boundary   = count_grids_local_boundary;
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;
		glob1->timers[timer_to_set].started = 1;

		size_t packsize = fclaw2d_global_packsize(glob1);
		char buffer[packsize];
		int bytes_written = fclaw2d_global_pack(glob1, buffer);
		fclaw2d_global_t* glob2;
		int bytes_read = fclaw2d_global_unpack(buffer, &glob2);

		CHECK_EQ(bytes_read,                         bytes_written);

		CHECK_EQ(glob1->count_amr_advance,           glob2->count_amr_advance);
		CHECK_EQ(glob1->count_ghost_exchange,        glob2->count_ghost_exchange);
		CHECK_EQ(glob1->count_amr_regrid,            glob2->count_amr_regrid);
		CHECK_EQ(glob1->count_amr_new_domain,        glob2->count_amr_new_domain);
		CHECK_EQ(glob1->count_single_step,           glob2->count_single_step);
		CHECK_EQ(glob1->count_elliptic_grids,        glob2->count_elliptic_grids);
		CHECK_EQ(glob1->count_multiproc_corner,      glob2->count_multiproc_corner);
		CHECK_EQ(glob1->count_grids_per_proc,        glob2->count_grids_per_proc);
		CHECK_EQ(glob1->count_grids_remote_boundary, glob2->count_grids_remote_boundary);
		CHECK_EQ(glob1->count_grids_local_boundary,  glob2->count_grids_local_boundary);
		CHECK_EQ(glob1->curr_time,                   glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,                     glob2->curr_dt);
		for(int i = 0; i < FCLAW2D_TIMER_COUNT; i++){
			CHECK_EQ(glob1->timers[i].cumulative, glob2->timers[i].cumulative);
			CHECK_EQ(glob1->timers[i].running,    glob2->timers[i].running);
			CHECK_EQ(glob1->timers[i].started,    glob2->timers[i].started);
			CHECK_EQ(glob1->timers[i].stopped,    glob2->timers[i].stopped);
		}


		fclaw2d_global_destroy(glob1);
		fclaw2d_global_destroy(glob2);
	}
}
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
