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
#include <fclaw3dx_clawpatch46_fort.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_physical_bc.h>
#include <fclaw_global.h>
#include <fclaw_regrid.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw3d_map.h>
#include <fclaw3d_convenience.h>
#include <fclaw3d_metric.hpp>
#include <fclaw3d_metric.h>
#include <fclaw_options.h>
#include <test.hpp>
#include <test/test.hpp>
#include <fstream>
#include <bitset>

#include <fclaw2d_forestclaw.h>

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
        fopts.refratio=2;
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

        domain = fclaw3d_domain_new_unitcube(sc_MPI_COMM_WORLD, minlevel);
        fclaw_global_store_domain(glob, domain);

        //map = fclaw3d_map_new_nomap();
        //fclaw_global_store_map_3d(glob, map);

        fclaw2d_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

        fclaw_domain_data_new(glob->domain);
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

int clawpatch_idx(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    int mbc = clawpatch->mbc;
    int idx = (i+mbc) + (j+mbc)*(mx+2*mbc) + (k+mbc)*(mx+2*mbc)*(my+2*mbc) + m*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc);
    return idx;
}
double fill_function(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
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
void get_bounds_with_ghost(fclaw_patch_t* patch, 
                          int *i_start,
                          int *i_stop,
                          int *j_start,
                          int *j_stop,
                          int *k_start,
                          int *k_stop)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    //get physical boundaries
    int intersects_bc[6];
    for(int i=0;i<6;i++)
    {
        fclaw_patch_relation_t type = fclaw_patch_get_face_type(patch, i);
        intersects_bc[i] = type == FCLAW_PATCH_BOUNDARY;
    }
    *i_start = intersects_bc[0] ? 0 : -clawpatch->mbc;
    *i_stop  = intersects_bc[1] ? clawpatch->d3->mx : clawpatch->d3->mx + clawpatch->mbc;

    *j_start = intersects_bc[2] ? 0 : -clawpatch->mbc;
    *j_stop  = intersects_bc[3] ? clawpatch->d3->my : clawpatch->d3->my + clawpatch->mbc;

    *k_start = intersects_bc[4] ? 0 : -clawpatch->mbc;
    *k_stop  = intersects_bc[5] ? clawpatch->d3->mz : clawpatch->d3->mz + clawpatch->mbc;
}
void get_face_ghost_bounds(fclaw_patch_t* patch, 
                          int iface,
                          int *i_start,
                          int *i_stop,
                          int *j_start,
                          int *j_stop,
                          int *k_start,
                          int *k_stop)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    switch(iface)
    {
        case 0:
            *i_start = -clawpatch->mbc;
            *i_stop  = 0;
            *j_start = 0;
            *j_stop  = clawpatch->d3->my;
            *k_start = 0;
            *k_stop  = clawpatch->d3->mz;
            break;
        case 1:
            *i_start = clawpatch->d3->mx;
            *i_stop  = clawpatch->d3->mx + clawpatch->mbc;
            *j_start = 0;
            *j_stop  = clawpatch->d3->my;
            *k_start = 0;
            *k_stop  = clawpatch->d3->mz;
            break;
        case 2:
            *i_start = 0;
            *i_stop  = clawpatch->d3->mx;
            *j_start = -clawpatch->mbc;
            *j_stop  = 0;
            *k_start = 0;
            *k_stop  = clawpatch->d3->mz;
            break;
        case 3:
            *i_start = 0;
            *i_stop  = clawpatch->d3->mx;
            *j_start = clawpatch->d3->my;
            *j_stop  = clawpatch->d3->my + clawpatch->mbc;
            *k_start = 0;
            *k_stop  = clawpatch->d3->mz;
            break;
        case 4:
            *i_start = 0;
            *i_stop  = clawpatch->d3->mx;
            *j_start = 0;
            *j_stop  = clawpatch->d3->my;
            *k_start = -clawpatch->mbc;
            *k_stop  = 0;
            break;
        case 5:
            *i_start = 0;
            *i_stop  = clawpatch->d3->mx;
            *j_start = 0;
            *j_stop  = clawpatch->d3->my;
            *k_start = clawpatch->d3->mz;
            *k_stop  = clawpatch->d3->mz + clawpatch->mbc;
            break;
    }
}
void get_edge_ghost_bounds(fclaw_patch_t* patch, 
                          int iedge,
                          int *i_start,
                          int *i_stop,
                          int *j_start,
                          int *j_stop,
                          int *k_start,
                          int *k_stop)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int axis = iedge/4;
    int upper_1 = iedge&0b1;
    int upper_2 = iedge&0b10;
    if(axis == 0)
    {
        *i_start = 0;
        *i_stop = clawpatch->d3->mx;
        if(upper_1)
        {
            *j_start = clawpatch->d3->my;
            *j_stop  = clawpatch->d3->my + clawpatch->mbc;
        }
        else
        {
            *j_start = -clawpatch->mbc;
            *j_stop  = 0;
        }
        if(upper_2)
        {
            *k_start = clawpatch->d3->mz;
            *k_stop  = clawpatch->d3->mz + clawpatch->mbc;
        }
        else
        {
            *k_start = -clawpatch->mbc;
            *k_stop  = 0;
        }
    }
    else if(axis ==1)
    {
        *j_start = 0;
        *j_stop = clawpatch->d3->my;
        if(upper_1)
        {
            *i_start = clawpatch->d3->mx;
            *i_stop  = clawpatch->d3->mx + clawpatch->mbc;
        }
        else
        {
            *i_start = -clawpatch->mbc;
            *i_stop  = 0;
        }
        if(upper_2)
        {
            *k_start = clawpatch->d3->mz;
            *k_stop  = clawpatch->d3->mz + clawpatch->mbc;
        }
        else
        {
            *k_start = -clawpatch->mbc;
            *k_stop  = 0;
        }
    }
    else if(axis ==2)
    {
        *k_start = 0;
        *k_stop = clawpatch->d3->mz;
        if(upper_1)
        {
            *i_start = clawpatch->d3->mx;
            *i_stop  = clawpatch->d3->mx + clawpatch->mbc;
        }
        else
        {
            *i_start = -clawpatch->mbc;
            *i_stop  = 0;
        }
        if(upper_2)
        {
            *j_start = clawpatch->d3->my;
            *j_stop  = clawpatch->d3->my + clawpatch->mbc;
        }
        else
        {
            *j_start = -clawpatch->mbc;
            *j_stop  = 0;
        }
    }
}
void get_corner_ghost_bounds(fclaw_global_t* glob,
                             fclaw_patch_t* patch, 
                             int blockno,
                             int patchno,
                             int icorner,
                             int *i_start,
                             int *i_stop,
                             int *j_start,
                             int *j_stop,
                             int *k_start,
                             int *k_stop)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    fclaw_patch_relation_t type = fclaw_patch_get_corner_type(patch, icorner);
    if(type == FCLAW_PATCH_BOUNDARY)
    {
        //if corner is a boundary, then we need to determine if it touches an
        //actual boundary
        int faces[3];
        fclaw_domain_corner_faces(glob->domain, icorner, faces);
        int intersects_bdry[6];
        fclaw_physical_get_bc(glob,blockno,patchno,
                            intersects_bdry);
        int touches_boundary;
        for(int i=0; i<3 ;i++)
        {
            touches_boundary = intersects_bdry[faces[i]];
            if(touches_boundary)
            {
                break;
            }
        }
        if(!touches_boundary)
        {
            //ghost_cells will be checked by face_interpolate or edge_interpolate
            *i_start=0;
            *i_stop=0;
            *j_start=0;
            *j_stop=0;
            *k_start=0;
            *k_stop=0;
            return;
        }
        
    }
    int upper_x = icorner&0b1;
    int upper_y = icorner&0b10;
    int upper_z = icorner&0b100;
    if(upper_x)
    {
        *i_start = clawpatch->d3->mx;
        *i_stop  = clawpatch->d3->mx + clawpatch->mbc;
    }
    else
    {
        *i_start = -clawpatch->mbc;
        *i_stop  = 0;
    }
    if(upper_y)
    {
        *j_start = clawpatch->d3->my;
        *j_stop  = clawpatch->d3->my + clawpatch->mbc;
    }
    else
    {
        *j_start = -clawpatch->mbc;
        *j_stop  = 0;
    }
    if(upper_z)
    {
        *k_start = clawpatch->d3->mz;
        *k_stop  = clawpatch->d3->mz + clawpatch->mbc;
    }
    else 
    {
        *k_start = -clawpatch->mbc;
        *k_stop  = 0;
    }
}

bool interior_cell(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    return i >= 0 && i < mx && j >= 0 && j < my && k >= 0 && k < mz;
}
bool face_ghost_cell(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    return ((i < 0 || i >= mx) && j >= 0 && j < my && k >= 0 && k < mz) ||
           (i >= 0 && i < mx && (j < 0 || j >= my) && k >= 0 && k < mz) ||
           (i >= 0 && i < mx && j >= 0 && j < my && (k < 0 || k >= mz));
}
bool edge_ghost_cell(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    return ((i < 0 || i >= mx) && (j < 0 || j >= my) && k >= 0 && k < mz) ||
           ((i < 0 || i >= mx) && j >= 0 && j < my && (k < 0 || k >= mz)) ||
           (i >= 0 && i < mx && (j < 0 || j >= my) && (k < 0 || k >= mz));
}
bool corner_ghost_cell(fclaw_patch_t* patch, int i, int j, int k, int m)
{
    fclaw_clawpatch_t* clawpatch = fclaw2d_clawpatch_get_clawpatch(patch);
    int mx = clawpatch->d3->mx;
    int my = clawpatch->d3->my;
    int mz = clawpatch->d3->mz;
    return (i < 0 || i >= mx) && (j < 0 || j >= my) && (k < 0 || k >= mz);
}


struct iterate_t {
            std::string output_name;
            fclaw_global_t* error_glob;
            int num_incorrect_interior_cells=0;

            int num_incorrect_boundary_face_ghost_cells=0;
            int num_incorrect_copied_face_ghost_cells=0;
            int num_incorrect_averaged_face_ghost_cells=0;
            int num_incorrect_interpolated_face_ghost_cells=0;

            int num_incorrect_boundary_edge_ghost_cells=0;
            int num_incorrect_copied_edge_ghost_cells=0;
            int num_incorrect_averaged_edge_ghost_cells=0;
            int num_incorrect_interpolated_edge_ghost_cells=0;

            int num_incorrect_boundary_corner_ghost_cells=0;
            int num_incorrect_copied_corner_ghost_cells=0;
            int num_incorrect_averaged_corner_ghost_cells=0;
            int num_incorrect_interpolated_corner_ghost_cells=0;

            bool test_face_copy=false;
            bool test_face_average=false;
            bool test_face_interpolate=false;
            bool test_edge_copy=false;
            bool test_edge_average=false;
            bool test_edge_interpolate=false;
            bool test_corner_copy=false;
            bool test_corner_average=false;
            bool test_corner_interpolate=false;
            void reset_counters(){
            num_incorrect_boundary_face_ghost_cells=0;
            num_incorrect_copied_face_ghost_cells=0;
            num_incorrect_averaged_face_ghost_cells=0;
            num_incorrect_interpolated_face_ghost_cells=0;

            num_incorrect_boundary_edge_ghost_cells=0;
            num_incorrect_copied_edge_ghost_cells=0;
            num_incorrect_averaged_edge_ghost_cells=0;
            num_incorrect_interpolated_edge_ghost_cells=0;

            num_incorrect_boundary_corner_ghost_cells=0;
            num_incorrect_copied_corner_ghost_cells=0;
            num_incorrect_averaged_corner_ghost_cells=0;
            num_incorrect_interpolated_corner_ghost_cells=0;
            }
        };
}

void run_tests(iterate_t* iterate)
{
    int test_no= 0;
    for(int mx   : {8})
    for(int my   : {8})
    for(int mz   : {8})
    for(int mbc  : {2})
    {
        iterate->reset_counters();
        CubeDomain cube(2,4);

        cube.opts->d3->mx   = mx;
        cube.opts->d3->my   = my;
        cube.opts->d3->mz   = mz;
        cube.opts->mbc      = mbc;
        cube.opts->meqn     = 3;

        //refine on if interior patch
        auto tag4refinement 
          = [](fclaw_global_t* glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               int initflag)
        {
            if(patch->d3->xlower == doctest::Approx(0.5)&&
               patch->d3->ylower == doctest::Approx(0.5)&&
               patch->d3->zlower == doctest::Approx(0.5))
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
        fclaw_patch_vt(cube.glob)->physical_bc
            = [](struct fclaw_global *glob,
                 struct fclaw_patch *this_patch,
                 int blockno,
                 int patchno,
                 double t,
                 double dt,
                 int *intersects_bc,
                 int time_interp)
                 {
                    //do nothing
                 };




        fclaw_patch_vt(cube.glob)->tag4refinement = tag4refinement;
        fclaw_patch_vt(cube.glob)->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        fclaw_clawpatch_vt(cube.glob)->d3->fort_interpolate2fine = fort_interpolate2fine;

        CHECK_EQ(cube.glob->domain->global_num_patches, 64);
        cube.setup();
        CHECK_EQ(cube.glob->domain->global_num_patches, 127);

        //create output domain with bigger size, so that we can see ghost cells
        //in the vtk output
        CubeDomain cube_output(2,4);

        cube_output.opts->d3->mx   = mx+2*mbc;
        cube_output.opts->d3->my   = my+2*mbc;
        cube_output.opts->d3->mz   = mz+2*mbc;
        cube_output.opts->mbc      = mbc;
        cube_output.opts->meqn     = 3;

        fclaw_patch_vt(cube_output.glob)->tag4refinement = tag4refinement;
        fclaw_patch_vt(cube_output.glob)->tag4coarsening = tag4coarsening;
        // don't want to test interpolation
        fclaw_clawpatch_vt(cube_output.glob)->d3->fort_interpolate2fine = fort_interpolate2fine;
        CHECK_EQ(cube_output.glob->domain->global_num_patches, 64);
        cube_output.setup();
        CHECK_EQ(cube_output.glob->domain->global_num_patches, 127);

        //initialize patches
        fclaw_global_iterate_patches(
            cube.glob, 
            [](fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
            {
                fclaw_global_iterate_t* g = (fclaw_global_iterate_t*)user;
                fclaw_clawpatch_options_t* opts = fclaw_clawpatch_get_options(g->glob);
                // clear q
                double* q = fclaw_clawpatch_get_q(g->glob, patch);
                int size = fclaw_clawpatch_size(g->glob);
                for (int i = 0; i < size; i++)
                {
                    q[i] = 0;
                }
                //memset(q, 0, sizeof(double)*size);
                //loop over interior
                for(int m = 0; m < opts->meqn; m++)
                for(int k = 0; k < opts->d3->mz; k++)
                for(int j = 0; j < opts->d3->my; j++)
                for(int i = 0; i < opts->d3->mx; i++)
                {
                    int idx = clawpatch_idx(patch, i,j,k,m);
                    //fill with some different linear functions
                    q[idx] = fill_function(patch, i, j, k, m);
                }

            }, 
            nullptr
        );

        //check ghost cells
        iterate->error_glob = cube_output.glob;


        fclaw_patch_vtable_t* patch_vt = fclaw_patch_vt(cube.glob);
        if(!iterate->test_face_copy)
        {
            patch_vt->d3->copy_face 
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_face_average)
        {
            patch_vt->d3->average_face
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_face_interpolate)
        {
            patch_vt->d3->interpolate_face
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_edge_copy)
        {
            patch_vt->d3->copy_edge 
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_edge_average)
        {
            patch_vt->d3->average_edge
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_edge_interpolate)
        {
            patch_vt->d3->interpolate_edge
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_corner_copy)
        {
            patch_vt->d3->copy_corner 
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_corner_average)
        {
            patch_vt->d3->average_corner
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }
        if(!iterate->test_corner_interpolate)
        {
            patch_vt->d3->interpolate_corner
                = [](fclaw_global_t*,
                     fclaw_patch_t *,
                     fclaw_patch_t *,
                     int,int,int,int,
                     struct fclaw_patch_transform_data*
                    )
                     {
                        //do nothing
                     };
        }

        CHECK_EQ(cube.glob->domain->global_num_patches, 127);
        fclaw_ghost_update(cube.glob, 2, 4, 0, 0, FCLAW_TIMER_NONE);
        CHECK_EQ(cube.glob->domain->global_num_patches, 127);

        fclaw_global_iterate_patches(
            cube.glob, 
            [](fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
            {
                fclaw_global_iterate_t* g = (fclaw_global_iterate_t*)user;
                iterate_t* iterate = (iterate_t*)g->user;
                fclaw_clawpatch_options_t* opts = fclaw_clawpatch_get_options(g->glob);

                //clawpatch to store error in
                fclaw_patch_t* error_patch = &iterate->error_glob->domain->blocks[blockno].patches[patchno];  
                fclaw_clawpatch_t* error_clawpatch = fclaw2d_clawpatch_get_clawpatch(error_patch);

                int mbc = error_clawpatch->mbc;

                // get q
                double* q = fclaw_clawpatch_get_q(g->glob, patch);
                double* error_q = fclaw_clawpatch_get_q(iterate->error_glob, error_patch);

                //clear error q
                int size = fclaw_clawpatch_size(iterate->error_glob);
                memset(error_q, 0, sizeof(double)*size);

                //loop over interior
                for(int m = 0; m < opts->meqn; m++)
                for(int k = 0; k < opts->d3->mz; k++)
                for(int j = 0; j < opts->d3->my; j++)
                for(int i = 0; i < opts->d3->mx; i++)
                {
                    int idx = clawpatch_idx(patch, i,j,k,m);
                    int error_idx = clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m);
                    //get expected value
                    double expected = fill_function(patch, i, j, k, m);
                    error_q[error_idx] = abs(expected - q[idx]);
                    if(q[idx] != doctest::Approx(expected))
                    {
                        iterate->num_incorrect_interior_cells++;
                    }
                }
                for(int iface = 0; iface < 6; iface++)
                {
                    fclaw_patch_relation_t type = fclaw_patch_get_face_type(patch, iface);
                    int i_start, i_stop, j_start, j_stop, k_start, k_stop;
                    get_face_ghost_bounds(patch,iface,&i_start,&i_stop,&j_start,&j_stop,&k_start,&k_stop);
                    for(int m = 0; m < opts->meqn; m++)
                    for(int k = k_start; k < k_stop; k++)
                    for(int j = j_start; j < j_stop; j++)
                    for(int i = i_start; i < i_stop; i++)
                    {
                        double fill_val = fill_function(patch, i, j, k, m);
                        double val = fclaw_clawpatch_get_q(g->glob, patch)[clawpatch_idx(patch, i,j,k,m)];
                        if(type == FCLAW_PATCH_BOUNDARY)
                        {
                            //should be exactly zero
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = 0-val;
                            if(val != 0)
                            {
                                iterate->num_incorrect_boundary_face_ghost_cells++;
                            }
                        }
                        else if(type == FCLAW_PATCH_SAMESIZE)
                        {
                            double error = 0;
                            if(iterate->test_face_copy && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_copied_face_ghost_cells++;
                            }
                            if(!iterate->test_face_copy && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_copied_face_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_HALFSIZE)
                        {
                            double error = 0;
                            if(iterate->test_face_average && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_averaged_face_ghost_cells++;
                            }
                            if(!iterate->test_face_average && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_averaged_face_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_DOUBLESIZE)
                        {
                            double error = 0;
                            if(iterate->test_face_interpolate && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_interpolated_face_ghost_cells++;
                            }
                            if(!iterate->test_face_interpolate && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_interpolated_face_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                    }
                }
                for(int iedge = 0; iedge < 12; iedge++)
                {
                    fclaw_patch_relation_t type = fclaw_patch_get_edge_type(patch, iedge);
                    int i_start, i_stop, j_start, j_stop, k_start, k_stop;
                    get_edge_ghost_bounds(patch,iedge,&i_start,&i_stop,&j_start,&j_stop,&k_start,&k_stop);
                    for(int m = 0; m < opts->meqn; m++)
                    for(int k = k_start; k < k_stop; k++)
                    for(int j = j_start; j < j_stop; j++)
                    for(int i = i_start; i < i_stop; i++)
                    {
                        double fill_val = fill_function(patch, i, j, k, m);
                        double val = fclaw_clawpatch_get_q(g->glob, patch)[clawpatch_idx(patch, i,j,k,m)];
                        if(type == FCLAW_PATCH_BOUNDARY)
                        {
                            //should be exactly zero
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = 0-val;
                            if(val != 0)
                            {
                                iterate->num_incorrect_boundary_edge_ghost_cells++;
                            }
                        }
                        else if(type == FCLAW_PATCH_SAMESIZE)
                        {
                            double error = 0;
                            if(iterate->test_edge_copy && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_copied_edge_ghost_cells++;
                            }
                            if(!iterate->test_edge_copy && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_copied_edge_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_HALFSIZE)
                        {
                            double error = 0;
                            if(iterate->test_edge_average && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_averaged_edge_ghost_cells++;
                            }
                            if(!iterate->test_edge_average && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_averaged_edge_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_DOUBLESIZE)
                        {
                            double error = 0;
                            if(iterate->test_edge_interpolate && val != doctest::Approx(fill_val))
                            {
                                error = fill_val - val;
                                iterate->num_incorrect_interpolated_edge_ghost_cells++;
                            }
                            if(!iterate->test_edge_interpolate && val != 0)
                            {
                                error = 0 - val;
                                iterate->num_incorrect_interpolated_edge_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                    }
                }
                for(int icorner = 0; icorner < 8; icorner++)
                {
                    fclaw_patch_relation_t type = fclaw_patch_get_corner_type(patch, icorner);
                    int i_start, i_stop, j_start, j_stop, k_start, k_stop;
                    get_corner_ghost_bounds(g->glob,patch,blockno,patchno,icorner,&i_start,&i_stop,&j_start,&j_stop,&k_start,&k_stop);
                    for(int m = 0; m < opts->meqn; m++)
                    for(int k = k_start; k < k_stop; k++)
                    for(int j = j_start; j < j_stop; j++)
                    for(int i = i_start; i < i_stop; i++)
                    {
                        double fill_val = fill_function(patch, i, j, k, m);
                        double val = fclaw_clawpatch_get_q(g->glob, patch)[clawpatch_idx(patch, i,j,k,m)];
                        if(type == FCLAW_PATCH_BOUNDARY)
                        {
                            //should be exactly zero
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = abs(0-val);
                            if(val != 0)
                            {
                                iterate->num_incorrect_boundary_corner_ghost_cells++;
                            }
                        }
                        else if(type == FCLAW_PATCH_SAMESIZE)
                        {
                            double error = 0;
                            if(iterate->test_corner_copy && val != doctest::Approx(fill_val))
                            {
                                error = abs(fill_val - val);
                                iterate->num_incorrect_copied_corner_ghost_cells++;
                            }
                            if(!iterate->test_corner_copy && val != 0)
                            {
                                error = abs(0 - val);
                                iterate->num_incorrect_copied_corner_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_HALFSIZE)
                        {
                            double error = 0;
                            if(iterate->test_corner_average && val != doctest::Approx(fill_val))
                            {
                                error = abs(fill_val - val);
                                iterate->num_incorrect_averaged_corner_ghost_cells++;
                            }
                            if(!iterate->test_corner_average && val != 0)
                            {
                                error = abs(0 - val);
                                iterate->num_incorrect_averaged_corner_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                        else if(type == FCLAW_PATCH_DOUBLESIZE)
                        {
                            double error = 0;
                            if(iterate->test_corner_interpolate && val != doctest::Approx(fill_val))
                            {
                                error = abs(fill_val - val);
                                iterate->num_incorrect_interpolated_corner_ghost_cells++;
                            }
                            if(!iterate->test_corner_interpolate && val != 0)
                            {
                                error = abs(0 - val);
                                iterate->num_incorrect_interpolated_corner_ghost_cells++;
                            }
                            error_q[clawpatch_idx(error_patch, i+mbc,j+mbc,k+mbc,m)] = error;
                        }
                    }
                }

            }, 
            iterate
        );

        //check that patch was filled correctly
        CHECK_EQ(iterate->num_incorrect_interior_cells, 0);

        CHECK_EQ(iterate->num_incorrect_boundary_face_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_copied_face_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_averaged_face_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_interpolated_face_ghost_cells, 0);

        CHECK_EQ(iterate->num_incorrect_boundary_edge_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_copied_edge_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_averaged_edge_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_interpolated_edge_ghost_cells, 0);

        CHECK_EQ(iterate->num_incorrect_boundary_corner_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_copied_corner_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_averaged_corner_ghost_cells, 0);
        CHECK_EQ(iterate->num_incorrect_interpolated_corner_ghost_cells, 0);

        //if not write output
        if(test_output_vtk())
        {
            char test_no_str[5];
            snprintf(test_no_str, 5, "%04d", test_no);
            std::string filename = iterate->output_name+"_"+std::string(test_no_str);
            INFO("Test failed output error to " << filename << ".vtu");
            fclaw_clawpatch_output_vtpd_to_file(cube_output.glob,filename.c_str());
        }
        test_no++;

    }

}

TEST_CASE("3d clawpatch ghost fill copy face")
{
    iterate_t iterate;
    iterate.output_name = "3d_clawpatch_ghost_fill_copy_face";
    iterate.test_face_copy = true;
    run_tests(&iterate);
}

TEST_CASE("3d clawpatch ghost fill average face")
{
    iterate_t iterate;
    iterate.output_name = "3d_clawpatch_ghost_fill_average_face";
    iterate.test_face_average = true;
    run_tests(&iterate);
}

TEST_CASE("3d clawpatch ghost fill interpolate face")
{
    iterate_t iterate;
    iterate.output_name = "3d_clawpatch_ghost_fill_interpolate_face";
    iterate.test_face_interpolate = true;
    run_tests(&iterate);
}