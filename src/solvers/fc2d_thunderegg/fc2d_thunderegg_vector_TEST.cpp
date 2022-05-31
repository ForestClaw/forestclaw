/*
  Copyright (c) 2021-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc2d_thunderegg_vector.hpp"
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch_output_ascii.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_convenience.h>
#include <fclaw2d_metric.hpp>
#include <fclaw2d_metric.h>
#include <fclaw2d_options.h>
#include <test/doctest.h>
#include <test/test.hpp>
using namespace ThunderEgg;

namespace{
struct QuadDomain {
    fclaw2d_global_t* glob;
    fclaw_options_t fopts;
    fclaw2d_domain_t *domain;
    fclaw2d_clawpatch_options_t opts;

    QuadDomain(){
        fclaw2d_clawpatch_vtable_initialize(glob, 4);
        glob = fclaw2d_global_new();
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.minlevel=1;
        fopts.maxlevel=1;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;

        domain = create_test_domain(sc_MPI_COMM_WORLD,&fopts);
        fclaw2d_global_store_domain(glob, domain);
        fclaw2d_options_store(glob, &fopts);

        memset(&opts, 0, sizeof(opts));
        opts.mx   = 5;
        opts.my   = 6;
        opts.mbc  = 2;
        opts.meqn = 1;
        opts.maux = 1;
        opts.rhs_fields = 1;
        fclaw2d_clawpatch_options_store(glob, &opts);

        fclaw2d_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[1], 0, 1, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[2], 0, 2, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[3], 0, 3, &build_mode);
    }
    ~QuadDomain(){
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[1]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[2]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[3]);
        fclaw2d_global_destroy(glob);
    }
};
struct QuadDomainBrick {
    fclaw2d_global_t* glob;
    fclaw_options_t fopts;
    fclaw2d_domain_t *domain;
    fclaw2d_clawpatch_options_t opts;

    QuadDomainBrick(){
        fclaw2d_clawpatch_vtable_initialize(glob, 4);
        glob = fclaw2d_global_new();
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=2;
        fopts.mj=2;
        fopts.minlevel=0;
        fopts.maxlevel=0;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;

        domain = create_test_domain(sc_MPI_COMM_WORLD,&fopts);
        fclaw2d_global_store_domain(glob, domain);
        fclaw2d_options_store(glob, &fopts);

        memset(&opts, 0, sizeof(opts));
        opts.mx   = 5;
        opts.my   = 6;
        opts.mbc  = 2;
        opts.meqn = 1;
        opts.maux = 1;
        opts.rhs_fields = 1;
        fclaw2d_clawpatch_options_store(glob, &opts);

        fclaw2d_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[1].patches[0], 1, 0, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[2].patches[0], 2, 0, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[3].patches[0], 3, 0, &build_mode);
    }
    ~QuadDomainBrick(){
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[1].patches[0]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[2].patches[0]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[3].patches[0]);
        fclaw2d_global_destroy(glob);
    }
};
}
TEST_CASE("fclaw2d_thunderegg_get_vector")
{
    for(fc2d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        QuadDomain test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts.rhs_fields = meqn;
        }else{
            test_data.opts.meqn = meqn;
        }
        test_data.setup();

        //set data
        for(int i=0; i < 4; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw2d_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case SOLN:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*meqn; j++){
                data[j] = i*(mx+2*mbc)*(my+2*mbc)*meqn + j;
            }
        }

        // get vector
        Vector<2> vec = fc2d_thunderegg_get_vector(test_data.glob,data_choice);

        //CHECK
        CHECK(vec.getNumLocalPatches() == 4);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 2> view = vec.getPatchView(patch_idx);
            CHECK(view.getGhostStart()[0] == -mbc);
            CHECK(view.getGhostStart()[1] == -mbc);
            CHECK(view.getGhostStart()[2] == 0);
            CHECK(view.getStart()[0] == 0);
            CHECK(view.getStart()[1] == 0);
            CHECK(view.getStart()[2] == 0);
            CHECK(view.getEnd()[0] == mx-1);
            CHECK(view.getEnd()[1] == my-1);
            CHECK(view.getEnd()[2] == meqn-1);
            CHECK(view.getGhostEnd()[0] == mx-1+mbc);
            CHECK(view.getGhostEnd()[1] == my-1+mbc);
            CHECK(view.getGhostEnd()[2] == meqn-1);
            CHECK(view.getStrides()[0] == 1);
            CHECK(view.getStrides()[1] == (mx+2*mbc));
            CHECK(view.getStrides()[2] == (mx+2*mbc)*(my+2*mbc));
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int j=-mbc; j<my+mbc; j++){
                    for(int i=-mbc; i<mx+mbc; i++){
                        CHECK(view(i,j,eqn) == patch_idx*(mx+2*mbc)*(my+2*mbc)*meqn + idx);
                        idx++;
                    }
                }
            }
        }
    }
}
TEST_CASE("fclaw2d_thunderegg_store_vector")
{
    for(fc2d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        QuadDomain test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts.rhs_fields = meqn;
        }else{
            test_data.opts.meqn = meqn;
        }
        test_data.setup();


        //set data
        Communicator comm(MPI_COMM_WORLD);
        Vector<2> vec(comm,{mx,my},meqn,4,mbc);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 2> view = vec.getPatchView(patch_idx);
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int j=-mbc; j<my+mbc; j++){
                    for(int i=-mbc; i<mx+mbc; i++){
                        view(i,j,eqn) = patch_idx*(mx+2*mbc)*(my+2*mbc)*meqn + idx;
                        idx++;
                    }
                }
            }
        }

        fc2d_thunderegg_store_vector(test_data.glob,data_choice,vec);
        //check
        for(int i=0; i < 4; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw2d_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case SOLN:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*meqn; j++){
                CHECK(data[j] == i*(mx+2*mbc)*(my+2*mbc)*meqn + j);
            }
        }
    }

}
TEST_CASE("fclaw2d_thunderegg_get_vector multiblock")
{
    for(fc2d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        QuadDomainBrick test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts.rhs_fields = meqn;
        }else{
            test_data.opts.meqn = meqn;
        }
        test_data.setup();

        int mx = test_data.opts.mx;
        int my = test_data.opts.my;
        int mbc = test_data.opts.mbc;

        //set data
        for(int i=0; i < 4; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw2d_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case SOLN:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*meqn; j++){
                data[j] = i*(mx+2*mbc)*(my+2*mbc)*meqn + j;
            }
        }

        // get vector
        Vector<2> vec = fc2d_thunderegg_get_vector(test_data.glob,data_choice);

        //CHECK
        CHECK(vec.getNumLocalPatches() == 4);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 2> view = vec.getPatchView(patch_idx);
            CHECK(view.getGhostStart()[0] == -mbc);
            CHECK(view.getGhostStart()[1] == -mbc);
            CHECK(view.getGhostStart()[2] == 0);
            CHECK(view.getStart()[0] == 0);
            CHECK(view.getStart()[1] == 0);
            CHECK(view.getStart()[2] == 0);
            CHECK(view.getEnd()[0] == mx-1);
            CHECK(view.getEnd()[1] == my-1);
            CHECK(view.getEnd()[2] == meqn-1);
            CHECK(view.getGhostEnd()[0] == mx-1+mbc);
            CHECK(view.getGhostEnd()[1] == my-1+mbc);
            CHECK(view.getGhostEnd()[2] == meqn-1);
            CHECK(view.getStrides()[0] == 1);
            CHECK(view.getStrides()[1] == (mx+2*mbc));
            CHECK(view.getStrides()[2] == (mx+2*mbc)*(my+2*mbc));
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int j=-mbc; j<my+mbc; j++){
                    for(int i=-mbc; i<mx+mbc; i++){
                        CHECK(view(i,j,eqn) == patch_idx*(mx+2*mbc)*(my+2*mbc)*meqn + idx);
                        idx++;
                    }
                }
            }
        }
    }
}
TEST_CASE("fclaw2d_thunderegg_store_vector multiblock")
{
    for(fc2d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        QuadDomainBrick test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts.rhs_fields = meqn;
        }else{
            test_data.opts.meqn = meqn;
        }
        test_data.setup();


        //set data
        Communicator comm(MPI_COMM_WORLD);
        Vector<2> vec(comm,{mx,my},meqn,4,mbc);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 2> view = vec.getPatchView(patch_idx);
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int j=-mbc; j<my+mbc; j++){
                    for(int i=-mbc; i<mx+mbc; i++){
                        view(i,j,eqn) = patch_idx*(mx+2*mbc)*(my+2*mbc)*meqn + idx;
                        idx++;
                    }
                }
            }
        }

        fc2d_thunderegg_store_vector(test_data.glob,data_choice,vec);
        //check
        for(int i=0; i < 4; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw2d_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case SOLN:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw2d_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*meqn; j++){
                CHECK(data[j] == i*(mx+2*mbc)*(my+2*mbc)*meqn + j);
            }
        }
    }
}


