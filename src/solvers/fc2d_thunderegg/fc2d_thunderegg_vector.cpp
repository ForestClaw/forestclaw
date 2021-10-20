/*
  Copyright (c) 2019-2020 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"

using namespace ThunderEgg;

#if 0
/**
 * @brief Get the number of local cells in the forestclaw domain
 * 
 * @param glob the forestclaw glob
 * @return int the number of local cells
 */
static int get_num_local_cells(fclaw2d_global_t* glob){
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_domain_t *domain = glob->domain;

    return domain->local_num_patches * clawpatch_opt->mx * clawpatch_opt->my;
}
#endif

static void get_data(struct fclaw2d_global* glob, int patchno, fc2d_thunderegg_data_choice_t data_choice, double** q, int* meqn)
{
    switch(data_choice){
        case RHS:
          fclaw2d_clawpatch_rhs_data(glob, &glob->domain->blocks[0].patches[patchno], q, meqn);
        break;
        case SOLN:
          fclaw2d_clawpatch_soln_data(glob, &glob->domain->blocks[0].patches[patchno], q, meqn);
        break;
        case STORE_STATE:
          fclaw2d_clawpatch_soln_data(glob, &glob->domain->blocks[0].patches[patchno], q, meqn);
        break;
    }
}

ThunderEgg::Vector<2> fc2d_thunderegg_get_vector(struct fclaw2d_global *glob, fc2d_thunderegg_data_choice_t data_choice)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    std::array<int,3> ns;
    ns[0] = clawpatch_opt->mx;
    ns[1] = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    std::array<int,3> strides;
    strides[0] = 1;
    strides[1] = strides[0]*(ns[0] + 2 * mbc);
    strides[2] = strides[1]*(ns[1] + 2 * mbc);
    std::vector<double *> starts(glob->domain->local_num_patches);
    for(int i=0; i<glob->domain->local_num_patches; i++){
        get_data(glob, i, data_choice, &starts[i], &ns[2]);
    }
    Communicator comm(glob->mpicomm);
    return ThunderEgg::Vector<2>(comm,starts,strides,ns,mbc);
}
void fc2d_thunderegg_store_vector(struct fclaw2d_global *glob, fc2d_thunderegg_data_choice_t data_choice, const ThunderEgg::Vector<2>& vec)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;

    for(int patchno = 0; patchno < glob->domain->local_num_patches; patchno++){
        PatchView<const double, 2> view = vec.getPatchView(patchno);
        int meqn;
        double *data;
        get_data(glob, patchno, data_choice, &data, &meqn);
        for(int eqn = 0; eqn < meqn; eqn++){
            for(int j = -mbc; j < my + mbc; j++){
                for(int i = -mbc; i < mx + mbc; i++){
                    *data = view(i,j,eqn);
                    data++;
                }
            }
        }
    }
}

#if 0
void fc2d_thunderegg_vector(fclaw2d_global_t *glob, int data_choice)
{

    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    fclaw2d_domain_t *domain = glob->domain;

    mfields = clawpatch_opt->rhs_fields;

    ns[0] = clawpatch_opt->mx;
    ns[1] = clawpatch_opt->my;

    mbc = clawpatch_opt->mbc;

    strides[0] = 1;
    strides[1] = ns[0] + 2 * mbc;
    eqn_stride = strides[1] * (ns[1] + 2 * mbc);
    /*
    eqn_stride = 1;
    strides[0] = clawpatch_opt->meqn;
    strides[1] = strides[0] * (ns[0] + 2 * mbc);
    */

    patch_data.resize(domain->local_num_patches * clawpatch_opt->rhs_fields);
    switch (data_choice)
    {
        case RHS:
            fclaw2d_global_iterate_patches(glob, enumeratePatchData, this);
            break;
        case SOLN:
            /* Copy solution to start guess */
            fclaw2d_global_iterate_patches(glob, set_start_value, this);
            break;
#if 1            
        case STORE_STATE:
            /* Store phi at time level n for use in defining operator */
            fclaw2d_global_iterate_patches(glob, store_state, this);
            break;
#endif            
        default:
            fclaw_global_essentialf("fc2d_thunderegg_vector : no valid data_choice specified\n");
            exit(0);
            break;
    }    
}
void fc2d_thunderegg_vector::enumeratePatchData(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *patch,
                                               int blockno, int patchno,
                                               void *user) 
{
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *)user;


    fc2d_thunderegg_vector &vec = *(fc2d_thunderegg_vector *)g->user;

    int global_num, local_num, level;
    fclaw2d_patch_get_info(domain, patch, blockno, patchno, &global_num,
                           &local_num, &level);

    int mfields;
    double* rhs_ptr;
    fclaw2d_clawpatch_rhs_data(g->glob, patch, &rhs_ptr, &mfields);
    // update point to point to non ghost cell
    for(int c = 0; c < mfields; c++){
        vec.patch_data[local_num * mfields + c] = rhs_ptr + vec.eqn_stride * c +
                                                  vec.mbc * vec.strides[0] +
                                                  vec.mbc * vec.strides[1];
    }
}

void fc2d_thunderegg_vector::set_start_value(fclaw2d_domain_t *domain,
                                             fclaw2d_patch_t *patch,
                                             int blockno, int patchno,
                                             void *user) 
{
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *)user;


    fc2d_thunderegg_vector &vec = *(fc2d_thunderegg_vector *)g->user;

    int global_num, local_num, level;
    fclaw2d_patch_get_info(domain, patch, blockno, patchno, &global_num,
                           &local_num, &level);

    int meqn;
    double* q_ptr;
    fclaw2d_clawpatch_soln_data(g->glob, patch, &q_ptr, &meqn);

    int mfields;
    double* rhs_ptr;
    fclaw2d_clawpatch_rhs_data(g->glob, patch, &rhs_ptr, &mfields);
    FCLAW_ASSERT(mfields == meqn);

    // update point to point to non ghost cell
    for(int c = 0; c < meqn; c++)
    {
        vec.patch_data[local_num * meqn + c] = q_ptr + vec.eqn_stride * c +
                                                  vec.mbc * vec.strides[0] +
                                                  vec.mbc * vec.strides[1];
    }
}

#if 1
void fc2d_thunderegg_vector::store_state(fclaw2d_domain_t *domain,
                                         fclaw2d_patch_t *patch,
                                         int blockno, int patchno,
                                         void *user) 
{
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *) user;


    fc2d_thunderegg_vector &vec = *(fc2d_thunderegg_vector *) g->user;

    int global_num, local_num, level;
    fclaw2d_patch_get_info(domain, patch, blockno, patchno, &global_num,
                           &local_num, &level);

    int meqn;
    double* q_ptr;
    fclaw2d_clawpatch_soln_data(g->glob, patch, &q_ptr, &meqn);

    // Store points to state on each patch
    // update point to point to non ghost cell
    for(int c = 0; c < meqn; c++)
    {
        vec.patch_data[local_num * meqn + c] = q_ptr + vec.eqn_stride * c +
                                                  vec.mbc * vec.strides[0] +
                                                  vec.mbc * vec.strides[1];
    }
}
#endif


LocalData<2> fc2d_thunderegg_vector::getLocalDataPriv(int component_index, int local_patch_id) const 
{
    return LocalData<2>(patch_data[local_patch_id * mfields + component_index], strides, ns,mbc);
}
LocalData<2> fc2d_thunderegg_vector::getLocalData(int component_index, int local_patch_id) 
{
    return getLocalDataPriv(component_index, local_patch_id);
}
const LocalData<2> fc2d_thunderegg_vector::getLocalData(int component_index, int local_patch_id) const 
{
    return getLocalDataPriv(component_index, local_patch_id);
}
#endif
