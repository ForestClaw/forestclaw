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

#include <fclaw_base.h>
#include <fclaw_options.h>
#include <fclaw_packing.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <sc_mpi.h>
#include <test.hpp>

TEST_CASE("fclaw2d_options can store options in two seperate globs")
{
	fclaw2d_global_t* glob1 = fclaw2d_global_new();
	fclaw2d_global_t* glob2 = fclaw2d_global_new();

	fclaw_options_t* opts1 = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	fclaw_options_t* opts2 = FCLAW_ALLOC_ZERO(fclaw_options_t,1);

	fclaw2d_options_store(glob1, opts1);
	/* glob1 has one package glob2 has two */
	fclaw2d_options_store(glob2, opts2);

	CHECK_EQ(fclaw2d_get_options(glob1), opts1);
	CHECK_EQ(fclaw2d_get_options(glob2), opts2);

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw2d_get_options fails if not intialized")
{
	fclaw2d_global_t* glob1 = fclaw2d_global_new();
	fclaw2d_global_t* glob2 = fclaw2d_global_new();

	CHECK_SC_ABORTED(fclaw2d_get_options(glob1));

	CHECK_SC_ABORTED(fclaw2d_get_options(glob2));

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}

TEST_CASE("fclaw2d_options_store fails if called twice on a glob")
{
	fclaw2d_global_t* glob1 = fclaw2d_global_new();
	fclaw2d_global_t* glob2 = fclaw2d_global_new();

	fclaw2d_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw_options_t,1));
	CHECK_SC_ABORTED(fclaw2d_options_store(glob1, FCLAW_ALLOC_ZERO(fclaw_options_t,1)));

	fclaw2d_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw_options_t,1));
	CHECK_SC_ABORTED(fclaw2d_options_store(glob2, FCLAW_ALLOC_ZERO(fclaw_options_t,1)));

	fclaw2d_global_destroy(glob1);
	fclaw2d_global_destroy(glob2);
}
TEST_CASE("fclaw2d_options packing/unpacking")
{
	fclaw_options_t* opts = FCLAW_ALLOC_ZERO(fclaw_options_t,1);
	opts->dim = 3;
	opts->run_directory = "test";
	opts->initial_dt = 0.1;
	opts->tfinal = 1.0;
	opts->outstyle = 1;
	opts->nout = 10;
	opts->nstep = 100;
	opts->subcycle = 0;
	opts->use_fixed_dt = 1;
	opts->max_cfl = .2;
	opts->desired_cfl = .3;
	opts->reduce_cfl = 0;
//	double tout[] = {0.5,0.6,0.7,0.8,0.9};
//	opts->tout = tout;
	opts->refratio = 2;
	opts->minlevel = 4;
	opts->maxlevel = 3;
	opts->regrid_interval = 5;
	opts->smooth_refine = 2;
	opts->refine_threshold = 3.0;
	opts->time_sync = 0;
	opts->output_gauges = 1;
	opts->gauge_buffer_length = 300;
	opts->output_rays = 2;
	opts->manifold = 2;
	opts->mi = 3;
	opts->mj = 4;
	opts->periodic_x = 0;
	opts->periodic_y = 1;
	opts->flux_correction = 1;
	opts->fluctuation_correction = 2;
	opts->coarsen_delay = 3;
	opts->init_ghostcell = 0;
	opts->advance_one_step = 1;
	opts->outstyle_uses_maxlevel = 2;
	opts->timeinterp2fillghost = 3;
	opts->scale_string = "blah";
	double scale[3] = {1.0,2.0,3.0};
	opts->scale = scale;
	opts->shift_string = "jbkal";
	double shift[3] = {4.0,5.0,6.0};
	opts->shift = shift;
	opts->phi = 3.5;
	opts->ax = 3.98;
	opts->bx = 32.38;
	opts->ay = 32.90;
	opts->by = 32.90;
	opts->az = 104.9;
	opts->bz = 3290.9023;
	opts->run_user_diagnostics = 1;
	opts->compute_error = 2;
	opts->conservation_check = 3;
	opts->trapfpe = 3;
	opts->report_timing = 3;
	opts->report_timing_verbosity = 3;
	opts->mpi_debug = 2;
	opts->ghost_patch_pack_area = 3;
	opts->ghost_patch_pack_extra = 4;
	opts->ghost_patch_pack_numextrafields=3298;
	opts->verbosity = 3;
	opts->output = 3;
	opts->tikz_out = 3;
	opts->tikz_figsize_string ="jdfajfda";
	double figsize[2] = {3.0,4.0};
	opts->tikz_figsize = figsize;
	opts->tikz_plot_fig = 3;
	opts->tikz_plot_prefix = "plot2jk";
	opts->tikz_plot_suffix = "plot3jk";
	opts->tikz_mesh_only= 0;
	opts->prefix = "jdsjkl";
	opts->vtkspace = 0.328;
	opts->weighted_partition = 0;
	opts->is_registered = 1;
	opts->logging_prefix = "werqreqw";
	opts->unpacked = false;

	const fclaw_userdata_vtable_t* vt = fclaw_options_get_packing_vtable();

	size_t size = vt->size(opts);
	char buffer[size];
	size_t bytes_written = vt->pack(opts,buffer);
	REQUIRE_EQ(bytes_written,size);

	fclaw_options_t* output_opts = nullptr;
	size_t bytes_read = vt->unpack(buffer,(void**)&output_opts);

	REQUIRE_EQ(bytes_read,size);
	REQUIRE_NE(output_opts,nullptr);

	CHECK_EQ(opts->dim                                 , output_opts->dim);

	CHECK_NE(opts->run_directory                       , output_opts->run_directory);
	CHECK_UNARY(!strcmp(opts->run_directory, output_opts->run_directory));

	CHECK_EQ(opts->initial_dt                          , output_opts->initial_dt);
	CHECK_EQ(opts->tfinal                              , output_opts->tfinal);
	CHECK_EQ(opts->outstyle                            , output_opts->outstyle);
	CHECK_EQ(opts->nout                                , output_opts->nout);
	CHECK_EQ(opts->nstep                               , output_opts->nstep);
	CHECK_EQ(opts->subcycle                            , output_opts->subcycle);
	CHECK_EQ(opts->use_fixed_dt                        , output_opts->use_fixed_dt);
	CHECK_EQ(opts->max_cfl                             , output_opts->max_cfl);
	CHECK_EQ(opts->desired_cfl                         , output_opts->desired_cfl );
	CHECK_EQ(opts->reduce_cfl                          , output_opts->reduce_cfl);

//	CHECK_NE(opts->tout                                , output_opts->tout);
//	for(int i = 0; i < 5; i++)
//	{
//		CHECK_EQ(opts->tout[i],output_opts->tout[i]);
//	}

	CHECK_EQ(opts->refratio                            , output_opts->refratio);
	CHECK_EQ(opts->minlevel                            , output_opts->minlevel);
	CHECK_EQ(opts->maxlevel                            , output_opts->maxlevel);
	CHECK_EQ(opts->regrid_interval                     , output_opts->regrid_interval);
	CHECK_EQ(opts->smooth_refine                       , output_opts->smooth_refine);
	CHECK_EQ(opts->refine_threshold                    , output_opts->refine_threshold);
	CHECK_EQ(opts->time_sync                           , output_opts->time_sync);
	CHECK_EQ(opts->output_gauges                       , output_opts->output_gauges);
	CHECK_EQ(opts->gauge_buffer_length                 , output_opts->gauge_buffer_length);
	CHECK_EQ(opts->output_rays                         , output_opts->output_rays);
	CHECK_EQ(opts->manifold                            , output_opts->manifold);
	CHECK_EQ(opts->mi                                  , output_opts->mi);
	CHECK_EQ(opts->mj                                  , output_opts->mj);
	CHECK_EQ(opts->periodic_x                          , output_opts->periodic_x);
	CHECK_EQ(opts->periodic_y                          , output_opts->periodic_y);
	CHECK_EQ(opts->flux_correction                     , output_opts->flux_correction);
	CHECK_EQ(opts->fluctuation_correction              , output_opts->fluctuation_correction);
	CHECK_EQ(opts->coarsen_delay                       , output_opts->coarsen_delay);
	CHECK_EQ(opts->init_ghostcell                      , output_opts->init_ghostcell);
	CHECK_EQ(opts->advance_one_step                    , output_opts->advance_one_step);
	CHECK_EQ(opts->outstyle_uses_maxlevel              , output_opts->outstyle_uses_maxlevel);
	CHECK_EQ(opts->timeinterp2fillghost                , output_opts->timeinterp2fillghost);

	CHECK_NE(opts->scale_string                        , output_opts->scale_string);
	CHECK_UNARY(!strcmp(opts->scale_string, output_opts->scale_string));

	CHECK_NE(opts->scale                               , output_opts->scale);
	for(int i = 0; i < 3; i++)
	{
		CHECK_EQ(opts->scale[i],output_opts->scale[i]);
	}

	CHECK_NE(opts->shift_string                        , output_opts->shift_string);
	CHECK_UNARY(!strcmp(opts->shift_string, output_opts->shift_string));

	CHECK_NE(opts->shift                               , output_opts->shift);
	for(int i = 0; i < 3; i++)
	{
		CHECK_EQ(opts->shift[i],output_opts->shift[i]);
	}

	CHECK_EQ(opts->phi                                 , output_opts->phi);
	CHECK_EQ(opts->ax                                  , output_opts->ax);
	CHECK_EQ(opts->bx                                  , output_opts->bx);
	CHECK_EQ(opts->ay                                  , output_opts->ay);
	CHECK_EQ(opts->by                                  , output_opts->by);
	CHECK_EQ(opts->az                                  , output_opts->az);
	CHECK_EQ(opts->bz                                  , output_opts->bz);
	CHECK_EQ(opts->run_user_diagnostics                , output_opts->run_user_diagnostics);
	CHECK_EQ(opts->compute_error                       , output_opts->compute_error);
	CHECK_EQ(opts->conservation_check                  , output_opts->conservation_check);
	CHECK_EQ(opts->trapfpe                             , output_opts->trapfpe);
	CHECK_EQ(opts->report_timing                       , output_opts->report_timing);
	CHECK_EQ(opts->report_timing_verbosity             , output_opts->report_timing_verbosity);
	CHECK_EQ(opts->mpi_debug                           , output_opts->mpi_debug);
	CHECK_EQ(opts->ghost_patch_pack_area               , output_opts->ghost_patch_pack_area);
	CHECK_EQ(opts->ghost_patch_pack_extra              , output_opts->ghost_patch_pack_extra);
	CHECK_EQ(opts->ghost_patch_pack_numextrafields     , output_opts->ghost_patch_pack_numextrafields);
	CHECK_EQ(opts->verbosity                           , output_opts->verbosity);
	CHECK_EQ(opts->output                              , output_opts->output);
	CHECK_EQ(opts->tikz_out                            , output_opts->tikz_out);

	CHECK_NE(opts->tikz_figsize_string                 , output_opts->tikz_figsize_string);
	CHECK_UNARY(!strcmp(opts->tikz_figsize_string, output_opts->tikz_figsize_string));

	CHECK_NE(opts->tikz_figsize                        , output_opts->tikz_figsize);
	for(int i = 0; i < 2; i++)
	{
		CHECK_EQ(opts->tikz_figsize[i],output_opts->tikz_figsize[i]);
	}

	CHECK_EQ(opts->tikz_plot_fig                       , output_opts->tikz_plot_fig);

	CHECK_NE(opts->tikz_plot_prefix                    , output_opts->tikz_plot_prefix);
	CHECK_UNARY(!strcmp(opts->tikz_plot_prefix, output_opts->tikz_plot_prefix));

	CHECK_NE(opts->tikz_plot_suffix                    , output_opts->tikz_plot_suffix);
	CHECK_UNARY(!strcmp(opts->tikz_plot_suffix, output_opts->tikz_plot_suffix));

	CHECK_EQ(opts->tikz_mesh_only                      , output_opts->tikz_mesh_only);

	CHECK_NE(opts->prefix                              , output_opts->prefix);
	CHECK_UNARY(!strcmp(opts->prefix, output_opts->prefix));

	CHECK_EQ(opts->vtkspace                            , output_opts->vtkspace);
	CHECK_EQ(opts->weighted_partition                  , output_opts->weighted_partition);
	CHECK_EQ(opts->is_registered                       , output_opts->is_registered);

	CHECK_NE(opts->logging_prefix                      , output_opts->logging_prefix);
	CHECK_UNARY(!strcmp(opts->logging_prefix, output_opts->logging_prefix));

	CHECK_UNARY(opts->unpacked);
}

#endif