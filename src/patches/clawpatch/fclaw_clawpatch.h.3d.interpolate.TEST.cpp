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

#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch.h>
#include <fclaw_ghost_fill.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_vtk.h>
#include <fclaw3d_clawpatch46_fort.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw3d_map.h>
#include <fclaw_convenience.h>
#include <fclaw3d_metric.hpp>
#include <fclaw3d_metric.h>
#include <fclaw_options.h>
#include <fclaw_regrid.h>
#include <test.hpp>
#include <test/test.hpp>
#include <fstream>
#include <bitset>

#include <fclaw_forestclaw.h>

namespace {
struct CubeDomain {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_domain_t *domain;
    fclaw3d_map_context_t* map;
    fclaw_clawpatch_options_t* opts;

    CubeDomain(int minlevel, int maxlevel){
        glob = fclaw_global_new();
        opts = fclaw_clawpatch_options_new(3);
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.minlevel=minlevel;
        fopts.maxlevel=maxlevel;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 1;
        fopts.bz = 1;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fopts.init_ghostcell = false;
        fclaw_options_store(glob, &fopts);

        opts->d3->mx   = 5;
        opts->d3->my   = 6;
        opts->d3->mz   = 7;
        opts->mbc  = 2;
        opts->meqn = 1;
        opts->maux = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        fclaw_global_store_domain(glob, domain);

        //map = fclaw3d_map_new_nomap();
        //fclaw_global_store_map_3d(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

            }
    void setup(){
        fclaw_initialize(glob);
    }
    ~CubeDomain(){
        fclaw_clawpatch_options_destroy(opts);
        //fclaw_domain_destroy(domain);
        //fclaw3d_map_destroy(map);
        fclaw_global_destroy(glob);
    }
};

int clawpatch_idx(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    int mbc = clawpatch->mbc;
    int idx = (i+mbc) + (j+mbc)*(mx+2*mbc) + (k+mbc)*(mx+2*mbc)*(my+2*mbc) + m*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc);
    return idx;
}
double fill_function(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    double x = clawpatch->d3->xlower + (i+0.5)*clawpatch->d3->dx;
    double y = clawpatch->d3->ylower + (j+0.5)*clawpatch->d3->dy;
    double z = clawpatch->d3->zlower + (k+0.5)*clawpatch->d3->dz;
    switch(m)
    {
        case 0:
            return 1;
        case 1:
            return 8*x + 2*y + 3*z+1;
        case 2:
            return 9*x + 10*y + 2*z+1;
    }
    return 0;
}

}// end anonymous namespace


TEST_CASE("3d clawpatch interpolate")
{
    int test_no= 0;
    for(int mx   : {8})
    for(int my   : {8})
    for(int mz   : {8})
    for(int mbc  : {2})
    {
        CubeDomain cube(1,2);

        cube.opts->d3->mx   = mx;
        cube.opts->d3->my   = my;
        cube.opts->d3->mz   = mz;
        cube.opts->mbc      = mbc;
        cube.opts->meqn     = 3;

        //don't refine on initialize
        fclaw_patch_vtable_t* cube_patch_vt = fclaw_patch_vt(cube.glob);
        cube_patch_vt->tag4refinement = 
            [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            return (int) !initflag;
        };
        cube_patch_vt->tag4coarsening = 
            [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            return 0;
        };

        cube.setup();
        CHECK_EQ(cube.glob->domain->global_num_patches, 8);

        //create output domain with bigger size, so that we can see ghost cells
        //in the vtk output
        CubeDomain cube_output(2,2);

        cube_output.opts->d3->mx   = mx+2*mbc;
        cube_output.opts->d3->my   = my+2*mbc;
        cube_output.opts->d3->mz   = mz+2*mbc;
        cube_output.opts->mbc      = mbc;
        cube_output.opts->meqn     = 3;

        cube_output.setup();

        //initialize patches
        fclaw_global_iterate_patches(
            cube.glob, 
            [](fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
            {
                fclaw_global_iterate_t* g = (fclaw_global_iterate_t*)user;
                fclaw_clawpatch_options_t* opts = fclaw_clawpatch_get_options(g->glob);
                fclaw_clawpatch_t* clawpatch = fclaw_clawpatch_get_clawpatch(patch);
                // clear q
                double* q = fclaw_clawpatch_get_q(g->glob, patch);
                int size = fclaw_clawpatch_size(g->glob);
                memset(q, 0, sizeof(double)*size);
                //loop over interior
                for(int m = 0; m < opts->meqn; m++)
                for(int k = -opts->mbc; k < opts->d3->mz+opts->mbc; k++)
                for(int j = -opts->mbc; j < opts->d3->my+opts->mbc; j++)
                for(int i = -opts->mbc; i < opts->d3->mx+opts->mbc; i++)
                {
                    int idx = clawpatch_idx(clawpatch, i,j,k,m);
                    //fill with some different linear functions
                    q[idx] = fill_function(clawpatch, i, j, k, m);
                }

            }, 
            nullptr
        );

        fclaw_regrid(cube.glob);
        CHECK_EQ(cube.glob->domain->global_num_patches, 64);

        //check ghost cells
        //fill output domain with error
        struct iterate_t {
            fclaw_global_t* error_glob;
            int num_incorrect_interior_cells;
        };
        iterate_t iterate = {cube_output.glob, 0};
        fclaw_global_iterate_patches(
            cube.glob, 
            [](fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
            {
                fclaw_global_iterate_t* g = (fclaw_global_iterate_t*)user;
                iterate_t* iterate = (iterate_t*)g->user;
                //clawpatch to store error in
                fclaw_patch_t* error_patch = &iterate->error_glob->domain->blocks[blockno].patches[patchno];  
                fclaw_clawpatch_t* error_clawpatch = fclaw_clawpatch_get_clawpatch(error_patch);
                int mbc = error_clawpatch->mbc;


                fclaw_clawpatch_options_t* opts = fclaw_clawpatch_get_options(g->glob);
                fclaw_clawpatch_t* clawpatch = fclaw_clawpatch_get_clawpatch(patch);

                // get q
                double* q = fclaw_clawpatch_get_q(g->glob, patch);
                double* error_q = fclaw_clawpatch_get_q(iterate->error_glob, error_patch);

                //get physical boundaries
                int intersects_bc[fclaw_domain_num_faces(domain)];
                for(int i=0;i<fclaw_domain_num_faces(domain);i++)
                {
                    fclaw_patch_relation_t type = fclaw_patch_get_face_type(patch, i);
                    intersects_bc[i] = type == FCLAW_PATCH_BOUNDARY;
                }

                //clear error q
                int size = fclaw_clawpatch_size(iterate->error_glob);
                memset(error_q, 0, sizeof(double)*size);

                for(int m = 0; m < opts->meqn; m++)
                for(int k = 0; k < clawpatch->d3->mz; k++)
                for(int j = 0; j < clawpatch->d3->my; j++)
                for(int i = 0; i < clawpatch->d3->mx; i++)
                {
                    int idx = clawpatch_idx(clawpatch, i,j,k,m);
                    int error_idx = clawpatch_idx(error_clawpatch, i+mbc,j+mbc,k+mbc,m);
                    //get expected value
                    double expected = fill_function(clawpatch, i, j, k, m);
                    error_q[error_idx] = expected - q[idx];
                    if(q[idx] != doctest::Approx(expected))
                    {
                        iterate->num_incorrect_interior_cells++;
                    }
                }

            }, 
            &iterate
        );

        //check that patch was filled correctly
        CHECK_EQ(iterate.num_incorrect_interior_cells, 0);

        int num_incorrect_cells = iterate.num_incorrect_interior_cells;
        //if not write output
        if(test_output_vtk() && num_incorrect_cells >0)
        {
            char test_no_str[5];
            snprintf(test_no_str, 5, "%04d", test_no);
            std::string filename = "3d_clawpatch_interpolate_"+std::string(test_no_str);
            INFO("Test failed output error to " << filename << ".vtu");
            fclaw_clawpatch_output_vtk_to_file(cube_output.glob,filename.c_str());
        }
        test_no++;

    }

}

