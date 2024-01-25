/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_convenience.h>
#include <fclaw2d_convenience.h>
#include <fclaw3d_convenience.h>
#include <fclaw2d_wrap.h>
#include <fclaw3d_wrap.h>

fclaw_domain_t *fclaw_domain_new_unitsquare (sc_MPI_Comm mpicomm,
                                             int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_unitsquare(mpicomm, initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_2d_torus (sc_MPI_Comm mpicomm,
                                            int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_torus(mpicomm, initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_2d_twosphere (sc_MPI_Comm mpicomm,
                                               int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_twosphere(mpicomm, initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_2d_cubedsphere (sc_MPI_Comm mpicomm,
                                                 int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_cubedsphere(mpicomm, initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_2d_disk (sc_MPI_Comm mpicomm,
                                          int periodic_in_x,
                                          int periodic_in_y,
                                          int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_disk(mpicomm, periodic_in_x, periodic_in_y,
                                       initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_2d_brick (sc_MPI_Comm mpicomm,
                                           int blocks_in_x, int blocks_in_y,
                                           int periodic_in_x,
                                           int periodic_in_y,
                                           int initial_level)
{
    fclaw2d_domain_t* domain2d;
    domain2d = fclaw2d_domain_new_brick(mpicomm, blocks_in_x, blocks_in_y,
                                        periodic_in_x, periodic_in_y,
                                        initial_level);
    return fclaw_domain_wrap_2d(domain2d);
}

fclaw_domain_t *fclaw_domain_new_unitcube (sc_MPI_Comm mpicomm,
                                           int initial_level)
{
    fclaw3d_domain_t* domain3d;
    domain3d = fclaw3d_domain_new_unitcube(mpicomm,initial_level);
    return fclaw_domain_wrap_3d(domain3d);
}

fclaw_domain_t *fclaw_domain_new_3d_brick (sc_MPI_Comm mpicomm,
                                           int blocks_in_x, int blocks_in_y,
                                           int blocks_in_z,
                                           int periodic_in_x,
                                           int periodic_in_y,
                                           int periodic_in_z,
                                           int initial_level)
{
    fclaw3d_domain_t* domain3d;
    domain3d = fclaw3d_domain_new_brick(mpicomm, blocks_in_x, blocks_in_y, blocks_in_z,
                                        periodic_in_x, periodic_in_y, periodic_in_z,
                                        initial_level);
    return fclaw_domain_wrap_3d(domain3d);
}


void fclaw_domain_destroy (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_destroy(domain->d2);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_destroy(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
    for(int blockno = 0; blockno < domain->num_blocks; ++blockno)
    {
        fclaw_block_t *block = &domain->blocks[blockno];
        FCLAW_FREE(block->patches);
    }
    FCLAW_FREE(domain->blocks);
    FCLAW_FREE(domain->ghost_patches);

    FCLAW_FREE(domain);
}

fclaw_domain_t* fclaw_domain_adapt (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_t* new_domain = fclaw2d_domain_adapt(domain->d2);
        return fclaw_domain_wrap_2d(new_domain);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_t* new_domain = fclaw3d_domain_adapt(domain->d3);

        return fclaw_domain_wrap_3d(new_domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

fclaw_domain_t* fclaw_domain_partition (fclaw_domain_t * domain,
                                            int weight_exponent)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_t* new_domain;
        new_domain = fclaw2d_domain_partition(domain->d2, weight_exponent);

        return fclaw_domain_wrap_2d(new_domain);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_t* new_domain;
        new_domain = fclaw3d_domain_partition(domain->d3, weight_exponent);

        return fclaw_domain_wrap_3d(new_domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_partition_unchanged (fclaw_domain_t * domain,
                                         int *unchanged_first,
                                         int *unchanged_length,
                                         int *unchanged_old_first)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_partition_unchanged(domain->d2,
                                           unchanged_first,
                                           unchanged_length,
                                           unchanged_old_first);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_partition_unchanged(domain->d3,
                                           unchanged_first,
                                           unchanged_length,
                                           unchanged_old_first);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_complete (fclaw_domain_t * domain)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_complete(domain->d2);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_complete(domain->d3);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_write_vtk (fclaw_domain_t * domain,
                             const char *basename)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_write_vtk(domain->d2,basename);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_write_vtk(domain->d3,basename);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_levels (fclaw_domain_t * domain, int log_priority)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_list_levels(domain->d2,log_priority);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_list_levels(domain->d3,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_neighbors (fclaw_domain_t * domain, int log_priority)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_list_neighbors(domain->d2,log_priority);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_list_neighbors(domain->d3,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_adapted (fclaw_domain_t * old_domain,
                                fclaw_domain_t * new_domain,
                                int log_priority)
{
    FCLAW_ASSERT(old_domain->refine_dim == new_domain->refine_dim);
    if(old_domain->refine_dim == 2)
    {
        fclaw2d_domain_list_adapted(old_domain->d2,new_domain->d2,log_priority);
    }
    else if(old_domain->refine_dim == 3)
    {
        fclaw3d_domain_list_adapted(old_domain->d3,new_domain->d3,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_search_points (fclaw_domain_t * domain,
                                 sc_array_t * block_offsets,
                                 sc_array_t * coordinates,
                                 sc_array_t * results)
{
    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_search_points(domain->d2,block_offsets,coordinates,results);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_search_points(domain->d3,block_offsets,coordinates,results);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_integrate_rays (fclaw_domain_t * domain,
                                  fclaw_integrate_ray_t intersect,
                                  sc_array_t * rays,
                                  sc_array_t * integrals,
                                  void * user)
{
    fclaw_integrate_ray_wrap_user_t wrap_user;
    wrap_user.intersect = intersect;
    wrap_user.user = user;

    if(domain->refine_dim == 2)
    {
        fclaw2d_domain_integrate_rays(domain->d2,fclaw2d_intersect_wrap_cb,
                                      rays,integrals,&wrap_user);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_domain_integrate_rays(domain->d3,fclaw3d_intersect_wrap_cb,
                                      rays,integrals,&wrap_user);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_overlap_exchange (fclaw_domain_t * domain,
                             sc_array_t * query_points,
                             fclaw_interpolate_point_t interpolate,
                             void *user)
{
    fclaw_interpolate_point_user_wrap_t wrap_user;
    wrap_user.interpolate = interpolate;
    wrap_user.user = user;

    if(domain->refine_dim == 2)
    {
        fclaw2d_overlap_exchange(domain->d2,query_points,fclaw2d_interpolate_point_wrap_cb,&wrap_user);
    }
    else if(domain->refine_dim == 3)
    {
        fclaw3d_overlap_exchange(domain->d3,query_points,fclaw3d_interpolate_point_wrap_cb,&wrap_user);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}