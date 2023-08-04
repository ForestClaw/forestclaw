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

#include <test.hpp>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_global.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>
#include <fclaw3dx_clawpatch_fort.h>

namespace{
    struct exceeds_test_parameters
    {
        int *blockno          = (int*)    1;
        int* meqn             = (int*)    2; 
        double *qval          = (double*) 3;
        double* qmin          = (double*) 4;
        double *qmax          = (double*) 5;
        double *quad          = (double*) 6;
        double *dx            = (double*) 7;
        double *dy            = (double*) 8;
        double *dz            = (double*) 9;
        double *xc            = (double*)10;
        double *yc            = (double*)11;
        double *zc            = (double*)12;
        int* ivar_threshold   = (int*)   13; 
        double *tag_threshold = (double*)14;
        int *init_flag        = (int*)   15;
        int *is_ghost         = (int*)   16;
        int return_value      =          17;
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

    fclaw_clawpatch_vt(glob)->d2->fort_user_exceeds_threshold = 
        [](const int *blockno,
           const int* meqn,
           const double *qval, 
           const double* qmin, 
           const double *qmax,
           const double quad[], 
           const double *dx, 
           const double *dy, 
           const double *xc, 
           const double *yc, 
           const int* ivar_threshold, 
           const double *tag_threshold,
           const int *init_flag,
           const int *is_ghost)
        {
            exceeds_test_parameters& params = global_exceeds_test_parameters;
            CHECK_EQ(blockno, params.blockno);
            //CHECK_EQ(meqn, params.meqn);  // This is changed inside function
            CHECK_EQ(qval, params.qval);
            CHECK_EQ(qmin, params.qmin);
            CHECK_EQ(qmax, params.qmax);
            CHECK_EQ(quad, params.quad);
            CHECK_EQ(dx, params.dx);
            CHECK_EQ(dy, params.dy);
            CHECK_EQ(xc, params.xc);
            CHECK_EQ(yc, params.yc);
            // CHECK_EQ(ivar_threshold, params.ivar_threshold);    //Will change
            CHECK_EQ(tag_threshold, params.tag_threshold);
            CHECK_EQ(init_flag, params.init_flag);
            CHECK_EQ(is_ghost, params.is_ghost);
            return params.return_value;
        };

	fclaw_clawpatch_options_t opts;
    opts.refinement_criteria = FCLAW_REFINE_CRITERIA_USER;
    fclaw_clawpatch_options_store(glob, &opts);

    fclaw2d_global_set_global(glob);
    int ret = FCLAW2D_CLAWPATCH_TAG_CRITERIA(params.blockno,
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
TEST_CASE("FCLAW3DX_CLAWPATCH_EXCEEDS_THRESHOLD calls user function")
{
    global_exceeds_test_parameters = exceeds_test_parameters();
    exceeds_test_parameters& params = global_exceeds_test_parameters;

	fclaw2d_global_t* glob = fclaw2d_global_new();

    fclaw2d_vtables_initialize(glob);
    fclaw3d_clawpatch_vtable_initialize(glob, 4);

    fclaw_clawpatch_vt(glob)->d3->fort_user_exceeds_threshold = 
        [](const int *blockno,
           const int* meqn,
           const double *qval, 
           const double* qmin, 
           const double *qmax,
           const double quad[], 
           const double *dx, 
           const double *dy, 
           const double *dz, 
           const double *xc, 
           const double *yc, 
           const double *zc, 
           const int* ivar_threshold, 
           const double *tag_threshold,
           const int *init_flag,
           const int *is_ghost)
        {
            exceeds_test_parameters& params = global_exceeds_test_parameters;
            CHECK_EQ(blockno, params.blockno);
            //CHECK_EQ(meqn, params.meqn);
            CHECK_EQ(qval, params.qval);
            CHECK_EQ(qmin, params.qmin);
            CHECK_EQ(qmax, params.qmax);
            CHECK_EQ(quad, params.quad);
            CHECK_EQ(dx, params.dx);
            CHECK_EQ(dy, params.dy);
            CHECK_EQ(dz, params.dz);
            CHECK_EQ(xc, params.xc);
            CHECK_EQ(yc, params.yc);
            CHECK_EQ(zc, params.zc);
            //CHECK_EQ(ivar_threshold, params.ivar_threshold); 
            CHECK_EQ(tag_threshold, params.tag_threshold); 
            CHECK_EQ(init_flag, params.init_flag);
            CHECK_EQ(is_ghost, params.is_ghost);
            return params.return_value;
        };

	fclaw_clawpatch_options_t opts;
    opts.refinement_criteria = FCLAW_REFINE_CRITERIA_USER;
    fclaw_clawpatch_options_store(glob, &opts);

    fclaw2d_global_set_global(glob);
    int ret = FCLAW3DX_CLAWPATCH_TAG_CRITERIA(params.blockno,
                                                   params.qval,
                                                   params.qmin,
                                                   params.qmax,
                                                   params.quad,
                                                   params.dx,
                                                   params.dy,
                                                   params.dz,
                                                   params.xc,
                                                   params.yc,
                                                   params.zc,
                                                   params.tag_threshold,
                                                   params.init_flag,
                                                   params.is_ghost);
    fclaw2d_global_unset_global();

    CHECK_EQ(ret, params.return_value);

    fclaw2d_global_destroy(glob);
}