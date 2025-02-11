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

#include "fc3d_thunderegg_vector.hpp"
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch_options.h>
#include <fclaw3d_clawpatch46_fort.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_diagnostics.h>
#include <fclaw_map_brick.h>
#include <fclaw_forestclaw.h>
#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw_convenience.h>
#include <fclaw3d_metric.hpp>
#include <fclaw3d_metric.h>
#include <fclaw_map.h>
#include <fclaw_options.h>
#include <test.hpp>
#include <test/test.hpp>
using namespace ThunderEgg;

namespace{
struct OctDomain {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_map_context_t    *map;
    fclaw_domain_t  *domain;
    fclaw_clawpatch_options_t* opts;

    OctDomain(){
        int size, rank;
        sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &size);
        sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
        glob = fclaw_global_new_comm(sc_MPI_COMM_WORLD, size, rank);
        opts = fclaw_clawpatch_options_new(2);
    
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.mk=1;
        fopts.minlevel=1;
        fopts.maxlevel=1;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fclaw_options_store(glob, &fopts);

        opts->patch_dim = 3;
        opts->mx     = 5;
        opts->my     = 6;
        opts->mbc        = 2;
        opts->meqn       = 1;
        opts->maux       = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, 1);
        fclaw_global_store_domain(glob, domain);

        map = fclaw_map_new_nomap();
        fclaw_map_store(glob,map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

        
    }
    void setup(){
        fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[1], 0, 1, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[2], 0, 2, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[3], 0, 3, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[4], 0, 4, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[5], 0, 5, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[6], 0, 6, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[7], 0, 7, &build_mode);
    }
    ~OctDomain(){
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[1]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[2]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[3]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[4]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[5]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[6]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[7]);
        fclaw_clawpatch_options_destroy(opts);
        fclaw_map_destroy(map);
        fclaw_domain_destroy(domain);
        fclaw_global_destroy(glob);
    }
};
struct OctDomainBrick {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_map_context_t* map;
    fclaw_domain_t *domain;
    fclaw_clawpatch_options_t* opts;

    OctDomainBrick(){
        int size, rank;
        sc_MPI_Comm_size(sc_MPI_COMM_WORLD, &size);
        sc_MPI_Comm_rank(sc_MPI_COMM_WORLD, &rank);
        glob = fclaw_global_new_comm(sc_MPI_COMM_WORLD, size, rank);
        opts = fclaw_clawpatch_options_new(2);

        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=2;
        fopts.mj=2;
        fopts.mk=2;
        fopts.minlevel=0;
        fopts.maxlevel=0;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fclaw_options_store(glob, &fopts);


        opts->patch_dim = 3;
        opts->mx     = 5;
        opts->my     = 6;
        opts->mbc        = 2;
        opts->meqn       = 1;
        opts->maux       = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        domain = fclaw_domain_new_3d_brick(sc_MPI_COMM_WORLD, fopts.mi, fopts.mj, fopts.mk, 0,0,0,0);
        fclaw_global_store_domain(glob, domain);

        map = fclaw_map_new_nomap();
        fclaw_map_store(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);
    
    }
    void setup(){
        fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;
        fclaw_patch_build(glob, domain, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[1].patches[0], 1, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[2].patches[0], 2, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[3].patches[0], 3, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[4].patches[0], 4, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[5].patches[0], 5, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[6].patches[0], 6, 0, &build_mode);
        fclaw_patch_build(glob, domain, &domain->blocks[7].patches[0], 7, 0, &build_mode);
    }
    ~OctDomainBrick(){
        fclaw_patch_data_delete(glob, domain, &domain->blocks[0].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[1].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[2].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[3].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[4].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[5].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[6].patches[0]);
        fclaw_patch_data_delete(glob, domain, &domain->blocks[7].patches[0]);
        fclaw_clawpatch_options_destroy(opts);
        fclaw_map_destroy(map);
        fclaw_domain_destroy(domain);
        fclaw_global_destroy(glob);
    }
};
}
TEST_CASE("fclaw3d_thunderegg_get_vector")
{
    for(fc3d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        OctDomain test_data;
        test_data.opts->mx = mx;
        test_data.opts->my = my;
        test_data.opts->mz = mz;
        test_data.opts->mbc    = mbc;
        if(data_choice == RHS){
            test_data.opts->rhs_fields = meqn;
        }else{
            test_data.opts->meqn = meqn;
        }
        test_data.setup();

        //set data
        for(int i=0; i < 8; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case SOLN:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn; j++){
                data[j] = i*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + j;
            }
        }

        // get vector
        Vector<3> vec = fc3d_thunderegg_get_vector(test_data.glob,data_choice);

        //CHECK
        CHECK(vec.getNumLocalPatches() == 8);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 3> view = vec.getPatchView(patch_idx);
            CHECK(view.getGhostStart()[0] == -mbc);
            CHECK(view.getGhostStart()[1] == -mbc);
            CHECK(view.getGhostStart()[2] == -mbc);
            CHECK(view.getGhostStart()[3] == 0);
            CHECK(view.getStart()[0] == 0);
            CHECK(view.getStart()[1] == 0);
            CHECK(view.getStart()[2] == 0);
            CHECK(view.getStart()[3] == 0);
            CHECK(view.getEnd()[0] == mx-1);
            CHECK(view.getEnd()[1] == my-1);
            CHECK(view.getEnd()[2] == mz-1);
            CHECK(view.getEnd()[3] == meqn-1);
            CHECK(view.getGhostEnd()[0] == mx-1+mbc);
            CHECK(view.getGhostEnd()[1] == my-1+mbc);
            CHECK(view.getGhostEnd()[2] == mz-1+mbc);
            CHECK(view.getGhostEnd()[3] == meqn-1);
            CHECK(view.getStrides()[0] == 1);
            CHECK(view.getStrides()[1] == (mx+2*mbc));
            CHECK(view.getStrides()[2] == (mx+2*mbc)*(my+2*mbc));
            CHECK(view.getStrides()[3] == (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc));
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int k=-mbc; k<mz+mbc; k++){
                    for(int j=-mbc; j<my+mbc; j++){
                        for(int i=-mbc; i<mx+mbc; i++){
                            CHECK(view(i,j,k,eqn) == patch_idx*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + idx);
                            idx++;
                        }
                    }
                }
            }
        }
    }
}
TEST_CASE("fclaw3d_thunderegg_store_vector")
{
    for(fc3d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        OctDomain test_data;
        test_data.opts->mx   = mx;
        test_data.opts->my   = my;
        test_data.opts->mz   = mz;
        test_data.opts->mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts->rhs_fields = meqn;
        }else{
            test_data.opts->meqn = meqn;
        }
        test_data.setup();


        //set data
        Communicator comm(MPI_COMM_WORLD);
        Vector<3> vec(comm,{mx,my,mz},meqn,8,mbc);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 3> view = vec.getPatchView(patch_idx);
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int k=-mbc; k<mz+mbc; k++){
                    for(int j=-mbc; j<my+mbc; j++){
                        for(int i=-mbc; i<mx+mbc; i++){
                            view(i,j,k,eqn) = patch_idx*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + idx;
                            idx++;
                        }
                    }
                }
            }
        }

        fc3d_thunderegg_store_vector(test_data.glob,data_choice,vec);
        //check
        for(int i=0; i < 8; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case SOLN:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[i], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn; j++){
                CHECK(data[j] == i*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + j);
            }
        }
    }

}
TEST_CASE("fclaw3d_thunderegg_get_vector multiblock")
{
    for(fc3d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        OctDomainBrick test_data;
        test_data.opts->mx   = mx;
        test_data.opts->my   = my;
        test_data.opts->mz   = mz;
        test_data.opts->mbc  = mbc;
        if(data_choice == RHS){
            test_data.opts->rhs_fields = meqn;
        }else{
            test_data.opts->meqn = meqn;
        }
        test_data.setup();

        //set data
        for(int i=0; i < 8; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case SOLN:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn; j++){
                data[j] = i*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + j;
            }
        }

        // get vector
        Vector<3> vec = fc3d_thunderegg_get_vector(test_data.glob,data_choice);

        //CHECK
        CHECK(vec.getNumLocalPatches() == 8);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 3> view = vec.getPatchView(patch_idx);
            CHECK(view.getGhostStart()[0] == -mbc);
            CHECK(view.getGhostStart()[1] == -mbc);
            CHECK(view.getGhostStart()[2] == -mbc);
            CHECK(view.getGhostStart()[3] == 0);
            CHECK(view.getStart()[0] == 0);
            CHECK(view.getStart()[1] == 0);
            CHECK(view.getStart()[2] == 0);
            CHECK(view.getStart()[3] == 0);
            CHECK(view.getEnd()[0] == mx-1);
            CHECK(view.getEnd()[1] == my-1);
            CHECK(view.getEnd()[2] == mz-1);
            CHECK(view.getEnd()[3] == meqn-1);
            CHECK(view.getGhostEnd()[0] == mx-1+mbc);
            CHECK(view.getGhostEnd()[1] == my-1+mbc);
            CHECK(view.getGhostEnd()[2] == mz-1+mbc);
            CHECK(view.getGhostEnd()[3] == meqn-1);
            CHECK(view.getStrides()[0] == 1);
            CHECK(view.getStrides()[1] == (mx+2*mbc));
            CHECK(view.getStrides()[2] == (mx+2*mbc)*(my+2*mbc));
            CHECK(view.getStrides()[3] == (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc));
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int k=-mbc; k<mz+mbc; k++){
                    for(int j=-mbc; j<my+mbc; j++){
                        for(int i=-mbc; i<mx+mbc; i++){
                            CHECK(view(i,j,k,eqn) == patch_idx*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + idx);
                            idx++;
                        }
                    }
                }
            }
        }
    }
}
TEST_CASE("fclaw3d_thunderegg_store_vector multiblock")
{
    for(fc3d_thunderegg_data_choice_t data_choice : {RHS,SOLN,STORE_STATE})
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        OctDomainBrick test_data;
        test_data.opts->mx = mx;
        test_data.opts->my = my;
        test_data.opts->mz = mz;
        test_data.opts->mbc    = mbc;
        if(data_choice == RHS){
            test_data.opts->rhs_fields = meqn;
        }else{
            test_data.opts->meqn = meqn;
        }
        test_data.setup();


        //set data
        Communicator comm(MPI_COMM_WORLD);
        Vector<3> vec(comm,{mx,my,mz},meqn,8,mbc);
        for(int patch_idx=0; patch_idx < vec.getNumLocalPatches(); patch_idx++){
            PatchView<double, 3> view = vec.getPatchView(patch_idx);
            int idx=0;
            for(int eqn=0; eqn<meqn; eqn++){
                for(int k=-mbc; k<mz+mbc; k++){
                    for(int j=-mbc; j<my+mbc; j++){
                        for(int i=-mbc; i<mx+mbc; i++){
                            view(i,j,k,eqn) = patch_idx*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + idx;
                            idx++;
                        }
                    }
                }
            }
        }

        fc3d_thunderegg_store_vector(test_data.glob,data_choice,vec);
        //check
        for(int i=0; i < 8; i++){
            double * data = nullptr;
            switch(data_choice){
                case RHS:
                  fclaw_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case SOLN:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
                case STORE_STATE:
                  fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[i].patches[0], &data, &meqn);
                break;
            }
            for(int j=0; j<(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn; j++){
                CHECK(data[j] == i*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn + j);
            }
        }
    }
}


