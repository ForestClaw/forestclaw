#include "fc2d_multigrid_vector.hpp"
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include "fc2d_multigrid.h"
#include "fc2d_multigrid_options.h"
fc2d_multigrid_vector::fc2d_multigrid_vector(fclaw2d_global_t *glob, int eqn) {
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    fclaw2d_domain_t *domain = glob->domain;

    this->eqn = eqn;

    ns[0] = clawpatch_opt->mx;
    ns[1] = clawpatch_opt->my;

    mbc = clawpatch_opt->mbc;

    strides[0] = 1;
    strides[1] = ns[0] + 2 * mbc;
    eqn_stride = strides[1] * (ns[1] + 2 * mbc);

    num_local_patches = domain->local_num_patches;
    patch_data.resize(num_local_patches);
    fclaw2d_global_iterate_patches(glob, enumeratePatchData, this);
}
void fc2d_multigrid_vector::enumeratePatchData(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *patch,
                                               int blockno, int patchno,
                                               void *user) {
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *)user;

    int meqn;

    int global_num, local_num, level;

    fc2d_multigrid_vector &vec = *(fc2d_multigrid_vector *)g->user;

    fclaw2d_patch_get_info(domain, patch, blockno, patchno, &global_num,
                           &local_num, &level);

    fclaw2d_clawpatch_soln_data(g->glob, patch, &vec.patch_data[local_num],
                                &meqn);
    // update point to point to non ghost cell
    vec.patch_data[local_num] += vec.eqn_stride * vec.eqn +
                                 vec.mbc * vec.strides[0] +
                                 vec.mbc * vec.strides[1];
}
LocalData<2> fc2d_multigrid_vector::getLocalDataPriv(int local_patch_id) const {
    return LocalData<2>(patch_data[local_patch_id], strides, ns);
}
LocalData<2> fc2d_multigrid_vector::getLocalData(int local_patch_id) {
    return getLocalDataPriv(local_patch_id);
}
const LocalData<2> fc2d_multigrid_vector::getLocalData(
    int local_patch_id) const {
    return getLocalDataPriv(local_patch_id);
}
