/*
  Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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
#include <fclaw_domain.h>
#include <fclaw_global.h>
#include <fclaw_patch.h>

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"

using namespace ThunderEgg;

static void get_data(struct fclaw_global* glob, fclaw_patch_t* patch, fc2d_thunderegg_data_choice_t data_choice, double** q, int* meqn)
{
    switch(data_choice){
        case RHS:
          fclaw2d_clawpatch_rhs_data(glob, patch, q, meqn);
        break;
        case SOLN:
          fclaw2d_clawpatch_soln_data(glob, patch, q, meqn);
        break;
        case STORE_STATE:
          fclaw2d_clawpatch_soln_data(glob, patch, q, meqn);
        break;
    }
}
ThunderEgg::Vector<2> fc2d_thunderegg_get_vector(struct fclaw_global *glob, fc2d_thunderegg_data_choice_t data_choice)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    std::array<int,3> ns;
    ns[0] = clawpatch_opt->mx;
    ns[1] = clawpatch_opt->my;
    switch(data_choice){
        case RHS:
          ns[2] = clawpatch_opt->rhs_fields;
        break;
        case SOLN:
          ns[2] = clawpatch_opt->meqn;
        break;
        case STORE_STATE:
          ns[2] = clawpatch_opt->meqn;
        break;
    }
    int mbc = clawpatch_opt->mbc;
    std::array<int,3> strides;
    strides[0] = 1;
    strides[1] = strides[0]*(ns[0] + 2 * mbc);
    strides[2] = strides[1]*(ns[1] + 2 * mbc);
    std::vector<double*> starts(glob->domain->local_num_patches);
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++){
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++){
            int meqn;
            get_data(glob, &block->patches[patchno], data_choice, &starts[block->num_patches_before+patchno], &meqn);
        }
    }
    Communicator comm(glob->mpicomm);
    return ThunderEgg::Vector<2>(comm,starts,strides,ns,mbc);
}
void fc2d_thunderegg_store_vector(struct fclaw_global *glob, fc2d_thunderegg_data_choice_t data_choice, const ThunderEgg::Vector<2>& vec)
{
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++){
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++){
            PatchView<const double, 2> view = vec.getPatchView(block->num_patches_before+patchno);
            int meqn;
            double *data;
            get_data(glob, &block->patches[patchno], data_choice, &data, &meqn);
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
}