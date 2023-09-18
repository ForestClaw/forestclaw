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
#include <fclaw_clawpatch_output_vtpd.h>
#include <fclaw3d_clawpatch46_fort.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_global.h>
#include <fclaw_regrid.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw3d_map.h>
#include <fclaw_convenience.h>
#include <fclaw3d_metric.hpp>
#include <fclaw3d_metric.h>
#include <fclaw_options.h>
#include <test.hpp>
#include <test/test.hpp>
#include <fstream>
#include <bitset>

#include <fclaw_forestclaw.h>

namespace {
struct TestData {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_domain_t *domain;
    fclaw3d_map_context_t* map;
    fclaw_clawpatch_options_t* opts;
    fclaw_patch_vtable_t* patch_vt;
    fclaw_clawpatch_vtable_t * clawpatch_vt;

    TestData(fclaw_domain_t* domain, int minlevel, int maxlevel){
        glob = fclaw_global_new();
        opts = fclaw_clawpatch_options_new(3);
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.mk=1;
        fopts.minlevel=minlevel;
        fopts.maxlevel=maxlevel;
        fopts.manifold=false;
        fopts.refratio=2;
        fopts.bx = 1;
        fopts.by = 1;
        fopts.bz = 1;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fopts.init_ghostcell = false;
        fclaw_options_store(glob, &fopts);

        opts->mx   = 5;
        opts->my   = 6;
        opts->mz   = 7;
        opts->mbc  = 2;
        opts->meqn = 1;
        opts->maux = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        fclaw_global_store_domain(glob, domain);

        //map = fclaw3d_map_new_nomap();
        //fclaw_global_store_map_3d(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);
        patch_vt = fclaw_patch_vt(glob);
        patch_vt->physical_bc = 
            [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               double t,
               double dt,
               int intersects_phys_bdry[],
               int time_interp)
            {
                //do nothing
            };

        clawpatch_vt = fclaw_clawpatch_vt(glob);

    }
    void setup(){
        fclaw_initialize(glob);
    }
    ~TestData(){
        fclaw_clawpatch_options_destroy(opts);
        //fclaw_domain_destroy(domain);
        //fclaw3d_map_destroy(map);
        fclaw_global_destroy(glob);
    }
};

int clawpatch_idx(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->mx;
    int my = clawpatch->my;
    int mz = clawpatch->mz;
    int mbc = clawpatch->mbc;
    int idx = (i+mbc) + (j+mbc)*(mx+2*mbc) + (k+mbc)*(mx+2*mbc)*(my+2*mbc) + m*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc);
    return idx;
}
double fill_function(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    double x = clawpatch->xlower + (i+0.5)*clawpatch->dx;
    double y = clawpatch->ylower + (j+0.5)*clawpatch->dy;
    double z = clawpatch->zlower + (k+0.5)*clawpatch->dz;
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
void get_bounds_with_ghost(fclaw_clawpatch_t* clawpatch, 
                          int *intersects_bc,
                          int *i_start,
                          int *i_stop,
                          int *j_start,
                          int *j_stop,
                          int *k_start,
                          int *k_stop)
{
    *i_start = intersects_bc[0] ? 0 : -clawpatch->mbc;
    *i_stop  = intersects_bc[1] ? clawpatch->mx : clawpatch->mx + clawpatch->mbc;

    *j_start = intersects_bc[2] ? 0 : -clawpatch->mbc;
    *j_stop  = intersects_bc[3] ? clawpatch->my : clawpatch->my + clawpatch->mbc;

    *k_start = intersects_bc[4] ? 0 : -clawpatch->mbc;
    *k_stop  = intersects_bc[5] ? clawpatch->mz : clawpatch->mz + clawpatch->mbc;
}

bool interior_cell(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->mx;
    int my = clawpatch->my;
    int mz = clawpatch->mz;
    return i >= 0 && i < mx && j >= 0 && j < my && k >= 0 && k < mz;
}
bool face_ghost_cell(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->mx;
    int my = clawpatch->my;
    int mz = clawpatch->mz;
    return ((i < 0 || i >= mx) && j >= 0 && j < my && k >= 0 && k < mz) ||
           (i >= 0 && i < mx && (j < 0 || j >= my) && k >= 0 && k < mz) ||
           (i >= 0 && i < mx && j >= 0 && j < my && (k < 0 || k >= mz));
}
bool edge_ghost_cell(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->mx;
    int my = clawpatch->my;
    int mz = clawpatch->mz;
    return ((i < 0 || i >= mx) && (j < 0 || j >= my) && k >= 0 && k < mz) ||
           ((i < 0 || i >= mx) && j >= 0 && j < my && (k < 0 || k >= mz)) ||
           (i >= 0 && i < mx && (j < 0 || j >= my) && (k < 0 || k >= mz));
}
bool corner_ghost_cell(fclaw_clawpatch_t* clawpatch, int i, int j, int k, int m)
{
    int mx = clawpatch->mx;
    int my = clawpatch->my;
    int mz = clawpatch->mz;
    return (i < 0 || i >= mx) && (j < 0 || j >= my) && (k < 0 || k >= mz);
}

void test_ghost_fill(TestData& cube, TestData& cube_output, std::string output_filename)
{
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
            for(int k = 0; k < opts->mz; k++)
            for(int j = 0; j < opts->my; j++)
            for(int i = 0; i < opts->mx; i++)
            {
                int idx = clawpatch_idx(clawpatch, i,j,k,m);
                //fill with some different linear functions
                q[idx] = fill_function(clawpatch, i, j, k, m);
            }

        }, 
        nullptr
    );

    fclaw_ghost_update(cube.glob, 2, 4, 0, 0, FCLAW_TIMER_NONE);

    //check ghost cells
    //fill output domain with error
    struct iterate_t {
        fclaw_global_t* error_glob;
        int num_incorrect_interior_cells;
        int num_incorrect_face_ghost_cells;
        int num_incorrect_edge_ghost_cells;
        int num_incorrect_corner_ghost_cells;
    };
    iterate_t iterate = {cube_output.glob, 0, 0, 0, 0};
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

            //loop over all cells
            int i_start, i_stop, j_start, j_stop, k_start, k_stop;
            get_bounds_with_ghost(clawpatch, intersects_bc,&i_start,&i_stop,&j_start,&j_stop,&k_start,&k_stop);
            for(int m = 0; m < opts->meqn; m++)
            for(int k = k_start; k < k_stop; k++)
            for(int j = j_start; j < j_stop; j++)
            for(int i = i_start; i < i_stop; i++)
            {
                int idx = clawpatch_idx(clawpatch, i,j,k,m);
                int error_idx = clawpatch_idx(error_clawpatch, i+mbc,j+mbc,k+mbc,m);
                //get expected value
                double expected = fill_function(clawpatch, i, j, k, m);
                error_q[error_idx] = expected - q[idx];
                if(q[idx] != doctest::Approx(expected))
                {
                    if(interior_cell(clawpatch, i,j,k,m))
                    {
                        iterate->num_incorrect_interior_cells++;
                    }
                    else if(face_ghost_cell(clawpatch, i,j,k,m))
                    {
                        iterate->num_incorrect_face_ghost_cells++;
                    }
                    else if(edge_ghost_cell(clawpatch, i,j,k,m))
                    {
                        iterate->num_incorrect_edge_ghost_cells++;
                    }
                    else if(corner_ghost_cell(clawpatch, i,j,k,m))
                    {
                        iterate->num_incorrect_corner_ghost_cells++;
                    }
                    else{
                        SC_ABORT_NOT_REACHED();
                    }
                }
            }

        }, 
        &iterate
    );

    //check that patch was filled correctly
    CHECK_EQ(iterate.num_incorrect_interior_cells, 0);
    CHECK_EQ(iterate.num_incorrect_face_ghost_cells, 0);
    CHECK_EQ(iterate.num_incorrect_edge_ghost_cells, 0);
    CHECK_EQ(iterate.num_incorrect_corner_ghost_cells, 0);

    //write output if --vtk option passed to test runner
    if(test_output_vtk())
    {
        INFO("Test failed output error to " << output_filename << ".vtpd");
        fclaw_clawpatch_output_vtpd_to_file(cube_output.glob,output_filename.c_str());
    }
}

}// end anonymous namespace


TEST_CASE("3d clawpatch ghost filling on uniform cube")
{
    int test_no= 0;
    for(int mx   : {10})
    for(int my   : {10})
    for(int mz   : {10})
    for(int mbc  : {1,2,3})
    {
        int minlevel = 2;
        int maxlevel = 2;
        fclaw_domain_t* domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data(domain,minlevel, maxlevel);

        test_data.opts->mx   = mx;
        test_data.opts->my   = my;
        test_data.opts->mz   = mz;
        test_data.opts->mbc      = mbc;
        test_data.opts->meqn     = 3;

        test_data.setup();

        //create output domain with bigger size, so that we can see ghost cells
        //in the vtk output
        fclaw_domain_t* domain_out = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data_out(domain_out, minlevel, maxlevel);

        test_data_out.opts->mx   = mx+2*mbc;
        test_data_out.opts->my   = my+2*mbc;
        test_data_out.opts->mz   = mz+2*mbc;
        test_data_out.opts->mbc      = mbc;
        test_data_out.opts->meqn     = 3;

        test_data_out.setup();

        
        char test_no_str[5];
        snprintf(test_no_str, 5, "%04d", test_no);
        std::string filename = "3d_ghost_fill_uniform_cube_"+std::string(test_no_str);

        test_ghost_fill(test_data, test_data_out, filename);

        test_no++;

    }

}

TEST_CASE("3d ghost fill on cube with refinement")
{
    int test_no= 0;
    for(int mx   : {8})
    for(int my   : {8})
    for(int mz   : {8})
    for(int mbc  : {2})
    {
        int minlevel = 2;
        int maxlevel = 4;
        fclaw_domain_t* domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data(domain,minlevel, maxlevel);

        test_data.opts->mx   = mx;
        test_data.opts->my   = my;
        test_data.opts->mz   = mz;
        test_data.opts->mbc      = mbc;
        test_data.opts->meqn     = 3;

        //refine on if interior patch
        auto tag4refinement 
          = [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            if(patch->xlower == doctest::Approx(0.5)&&
               patch->ylower == doctest::Approx(0.5)&&
               patch->zlower == doctest::Approx(0.5))
            {
                return 1;
            }
            else 
            {
                return 0;
            }
        };
        auto tag4coarsening
          = [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            return 0;
        };
        auto fort_interpolate2fine
          =[](const int* mx, 
              const int* my, 
              const int* mz,
              const int* mbc, 
              const int* meqn,
              double qcoarse[], 
              double qfine[],
              double areacoarse[], 
              double areafine[],
              const int* igrid, 
              const int* manifold)
              {
                //do nothing
              };

        test_data.patch_vt->tag4refinement = tag4refinement;
        test_data.patch_vt->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        test_data.clawpatch_vt->d3->fort_interpolate2fine = fort_interpolate2fine;

        CHECK_EQ(test_data.glob->domain->global_num_patches, 64);
        test_data.setup();
        CHECK_EQ(test_data.glob->domain->global_num_patches, 127);
        domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);

        //create output domain with bigger size, so that we can see ghost cells
        //in the vtk output
        fclaw_domain_t* domain_out = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data_out(domain_out, minlevel, maxlevel);

        test_data_out.opts->mx   = mx+2*mbc;
        test_data_out.opts->my   = my+2*mbc;
        test_data_out.opts->mz   = mz+2*mbc;
        test_data_out.opts->mbc      = mbc;
        test_data_out.opts->meqn     = 3;

        test_data_out.patch_vt->tag4refinement = tag4refinement;
        test_data_out.patch_vt->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        test_data_out.clawpatch_vt->d3->fort_interpolate2fine = fort_interpolate2fine;

        CHECK_EQ(test_data_out.glob->domain->global_num_patches, 64);
        test_data_out.setup();
        CHECK_EQ(test_data_out.glob->domain->global_num_patches, 127);

        char test_no_str[5];
        snprintf(test_no_str, 5, "%04d", test_no);
        std::string filename = "3d_ghost_fill_uniform_cube_"+std::string(test_no_str);

        test_ghost_fill(test_data, test_data_out, filename);

        test_no++;
    }

}
TEST_CASE("3d ghost fill on cube with refinement coarse interior")
{
    int test_no= 0;
    for(int mx   : {8})
    for(int my   : {8})
    for(int mz   : {8})
    for(int mbc  : {2})
    {
        int minlevel = 2;
        int maxlevel = 3;
        fclaw_domain_t* domain = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data(domain,minlevel, maxlevel);

        test_data.opts->mx   = mx;
        test_data.opts->my   = my;
        test_data.opts->mz   = mz;
        test_data.opts->mbc      = mbc;
        test_data.opts->meqn     = 3;

        //refine on if interior patch
        auto tag4refinement 
          = [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            if(patch->xlower == 0 ||
                patch->ylower == 0 ||
                patch->zlower == 0 ||
                patch->xupper == 1 ||
                patch->yupper == 1 ||
                patch->zupper == 1)
            {
                return 1;
            }
            else 
            {
                return 0;
            }
        };
        auto tag4coarsening
          = [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            return 0;
        };
        auto fort_interpolate2fine
          =[](const int* mx, 
              const int* my, 
              const int* mz,
              const int* mbc, 
              const int* meqn,
              double qcoarse[], 
              double qfine[],
              double areacoarse[], 
              double areafine[],
              const int* igrid, 
              const int* manifold)
              {
                //do nothing
              };

        test_data.patch_vt->tag4refinement = tag4refinement;
        test_data.patch_vt->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        test_data.clawpatch_vt->d3->fort_interpolate2fine = fort_interpolate2fine;

        CHECK_EQ(test_data.glob->domain->global_num_patches, 64);
        test_data.setup();
        CHECK_EQ(test_data.glob->domain->global_num_patches, 456);

        //create output domain with bigger size, so that we can see ghost cells
        //in the vtk output
        fclaw_domain_t* domain_out = fclaw_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        TestData test_data_out(domain_out, minlevel, maxlevel);

        test_data_out.opts->mx   = mx+2*mbc;
        test_data_out.opts->my   = my+2*mbc;
        test_data_out.opts->mz   = mz+2*mbc;
        test_data_out.opts->mbc      = mbc;
        test_data_out.opts->meqn     = 3;

        test_data_out.patch_vt->tag4refinement = tag4refinement;
        test_data_out.patch_vt->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        test_data_out.clawpatch_vt->d3->fort_interpolate2fine = fort_interpolate2fine;

        CHECK_EQ(test_data_out.glob->domain->global_num_patches, 64);
        test_data_out.setup();
        CHECK_EQ(test_data_out.glob->domain->global_num_patches, 456);

        char test_no_str[5];
        snprintf(test_no_str, 5, "%04d", test_no);
        std::string filename = "3d_ghost_fill_uniform_cube_"+std::string(test_no_str);

        test_ghost_fill(test_data, test_data_out, filename);

        test_no++;

    }

}