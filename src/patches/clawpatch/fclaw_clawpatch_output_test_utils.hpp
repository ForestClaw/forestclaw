/*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Scott Aiton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FCLAW_CLAWPATCH_OUTPUT_TEST_UTILS_HPP
#define FCLAW_CLAWPATCH_OUTPUT_TEST_UTILS_HPP

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_convenience.h>
#include <fclaw_domain.h>
#include <fclaw_forestclaw.h>
#include <fclaw_global.h>
#include <fclaw_map.h>
#include <fclaw_options.h>
#include <fclaw_patch.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>

namespace fclaw_clawpatch_output_test
{

constexpr float output_value = 2.5f;

class TestData
{
public:
    fclaw_global_t* glob;
    fclaw_domain_t* domain;
    fclaw_map_context_t* map;
    fclaw_options_t fopts;
    fclaw_clawpatch_options_t* opts;
    std::string basename;

    TestData()
    {
        const auto timestamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        basename = (std::filesystem::temp_directory_path() /
                    ("forestclaw-output-test-" + std::to_string(timestamp))).string();

        domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_SELF, 0);
        map = fclaw_map_new_nomap();
        glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
        opts = fclaw_clawpatch_options_new(2);

        std::memset(&fopts, 0, sizeof(fopts));
        fopts.minlevel = 0;
        fopts.maxlevel = 0;
        fopts.refratio = 2;
        fopts.bx = 1;
        fopts.by = 1;
        fopts.bz = 1;
        fopts.regression_check = "";
        fclaw_options_store(glob, &fopts);

        opts->mx = 2;
        opts->my = 2;
        opts->mbc = 1;
        opts->meqn = 1;
        fclaw_clawpatch_options_store(glob, opts);
        fclaw_global_store_domain(glob, domain);
        fclaw_map_store(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);
        fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(glob);
        patch_vt->physical_bc =
            [](fclaw_global_t*, fclaw_patch_t*, int, int, double, double,
               int[], int) {};

        fclaw_initialize(glob);
        fclaw_global_iterate_patches(
            glob,
            [](fclaw_domain_t*, fclaw_patch_t* patch, int, int, void* user)
            {
                auto* test_data = static_cast<TestData*>(user);
                double* q = fclaw_clawpatch_get_q(test_data->glob, patch);
                std::fill(q, q + fclaw_clawpatch_size(test_data->glob),
                          static_cast<double>(output_value));
            },
            this);
    }

    ~TestData()
    {
        std::remove((basename + ".vtu").c_str());
        std::remove((basename + ".vtkhdf").c_str());
        fclaw_finalize(glob);
        fclaw_global_destroy(glob);
        fclaw_clawpatch_options_destroy(opts);
    }
};

} // namespace fclaw_clawpatch_output_test

#endif
