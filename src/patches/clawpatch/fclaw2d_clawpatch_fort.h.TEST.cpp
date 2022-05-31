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

#include <doctest.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_global.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>

namespace{
    struct exceeds_test_parameters
    {
        int *blockno          = (int*)   1;
        double *qval          = (double*)2;
        double* qmin          = (double*)3;
        double *qmax          = (double*)4;
        double *quad          = (double*)5;
        double *dx            = (double*)6;
        double *dy            = (double*)7;
        double *xc            = (double*)8;
        double *yc            = (double*)9;
        double *tag_threshold = (double*)10;
        int *init_flag        = (int*)   11;
        int *is_ghost         = (int*)   12;
        int return_value      =          13;
    };

    exceeds_test_parameters global_exceeds_test_parameters;
}
TEST_CASE("FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD calls user function")
{
    global_exceeds_test_parameters = exceeds_test_parameters();
    exceeds_test_parameters& params = global_exceeds_test_parameters;

	fclaw2d_global_t* glob = fclaw2d_global_new();

    fclaw2d_vtables_initialize(glob);
    fclaw2d_clawpatch_vtable_initialize(glob, 4);

    fclaw2d_clawpatch_vt(glob)->fort_user_exceeds_threshold = 
        [](const int *blockno,
           const double *qval, 
           const double* qmin, 
           const double *qmax,
           const double quad[], 
           const double *dx, 
           const double *dy, 
           const double *xc, 
           const double *yc, 
           const double *tag_threshold,
           const int *init_flag,
           const int *is_ghost)
        {
            exceeds_test_parameters& params = global_exceeds_test_parameters;
            CHECK_EQ(blockno, params.blockno);
            CHECK_EQ(qval, params.qval);
            CHECK_EQ(qmin, params.qmin);
            CHECK_EQ(qmax, params.qmax);
            CHECK_EQ(quad, params.quad);
            CHECK_EQ(dx, params.dx);
            CHECK_EQ(dy, params.dy);
            CHECK_EQ(xc, params.xc);
            CHECK_EQ(yc, params.yc);
            CHECK_EQ(tag_threshold, params.tag_threshold);
            CHECK_EQ(init_flag, params.init_flag);
            CHECK_EQ(is_ghost, params.is_ghost);
            return params.return_value;
        };

	fclaw2d_clawpatch_options_t opts;
    opts.refinement_criteria = FCLAW_REFINE_CRITERIA_USER;
    fclaw2d_clawpatch_options_store(glob, &opts);

    fclaw2d_global_set_global(glob);
    int ret = FCLAW2D_CLAWPATCH_EXCEEDS_THRESHOLD(params.blockno,
                                                  params.qval,
                                                  params.qmin,
                                                  params.qmax,
                                                  params.quad,
                                                  params.dx,
                                                  params.dy,
                                                  params.xc,
                                                  params.yc,
                                                  params.tag_threshold,
                                                  params.init_flag,
                                                  params.is_ghost);
    fclaw2d_global_unset_global();

    CHECK_EQ(ret, params.return_value);

    fclaw2d_global_destroy(glob);
}