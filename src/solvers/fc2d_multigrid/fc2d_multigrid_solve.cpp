/*
Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#include "fc2d_multigrid.h"
#include "fc2d_multigrid_options.h"
#include "fc2d_multigrid_vector.hpp"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_output_ascii.h>
#include <fclaw2d_clawpatch_output_vtk.h>

#include <fclaw2d_global.h>
#include <fclaw2d_map.h>
#include <fclaw2d_options.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_vtable.h>

#include <p4est_bits.h>
#include <p4est_wrap.h>

#include <Thunderegg/BiCGStab.h>
#include <Thunderegg/BilinearInterpolator.h>
#include <Thunderegg/FivePtPatchOperator.h>
#include <Thunderegg/GMG/CycleFactory2d.h>
#include <Thunderegg/Operators/SchurDomainOp.h>
#include <Thunderegg/P4estDCG.h>
#include <Thunderegg/PatchSolvers/FftwPatchSolver.h>
using namespace std;
extern "C" {

void fc2d_multigrid_solve(fclaw2d_global_t *glob) {
    // this should be moved somewhere else
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

    // create thunderegg vector for eqn 0
    shared_ptr<Vector<2>> f(new fc2d_multigrid_vector(glob, 0));

    // get patch size
    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);
    array<int, 2> ns = {clawpatch_opt->mx, clawpatch_opt->my};

    // get p4est structure
    fclaw2d_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *)domain->pp;

    // create map function
    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    P4estDCG::BlockMapFunc bmf = [&](int block_no, double unit_x, double unit_y,
                                     double &x, double &y) {
        if (fclaw_opt->manifold) {
            // map?
        } else {
            x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * unit_x;
            y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * unit_y;
        }
    };

    // generates levels of patches for GMG
    shared_ptr<P4estDCG> dcg(new P4estDCG(wrap->p4est, ns, false, bmf));

    // get finest level
    shared_ptr<DomainCollection<2>> dc = dcg->getFinestDC();

    // define operators for problems
    // set the patch solver
    shared_ptr<PatchSolver<2>> solver(new FftwPatchSolver<2>(*dc));

    // patch operator
    shared_ptr<PatchOperator<2>> op(new FivePtPatchOperator());

    // interface interpolator
    shared_ptr<IfaceInterp<2>> interp(new BilinearInterpolator());

    // create SchurHelper
    shared_ptr<SchurHelper<2>> sh(new SchurHelper<2>(dc, solver, op, interp));

    // create matrix
    shared_ptr<Operator<2>> A(new SchurDomainOp<2>(sh));

    // create gmg preconditioner
    fc2d_multigrid_options_t *mg_opt = fc2d_multigrid_get_options(glob);
    GMG::CycleOpts copts;
    copts.max_levels = mg_opt->max_levels;
    copts.patches_per_proc = mg_opt->patches_per_proc;
    copts.pre_sweeps = mg_opt->pre_sweeps;
    copts.post_sweeps = mg_opt->post_sweeps;
    copts.mid_sweeps = mg_opt->mid_sweeps;
    copts.coarse_sweeps = mg_opt->coarse_sweeps;
    copts.cycle_type = mg_opt->cycle_type;
    shared_ptr<Operator<2>> M =
        GMG::CycleFactory2d::getCycle(copts, dcg, solver, op, interp);

    // solve
    shared_ptr<VectorGenerator<2>> vg(new DomainCollectionVG<2>(dc));
    shared_ptr<Vector<2>> u = vg->getNewVector();

    int its = BiCGStab<2>::solve(vg, A, u, f, M);

    fclaw_global_productionf("Iterations: %i\n", its);
    fclaw_global_productionf("f-2norm: %f\n", f->twoNorm());
    fclaw_global_productionf("f-infnorm: %f\n", f->infNorm());
    fclaw_global_productionf("u-2norm: %f\n", u->twoNorm());
    fclaw_global_productionf("u-infnorm: %f\n\n", u->infNorm());

    // copy solution into rhs
    f->copy(u);
    fclaw_global_productionf("Checking if copy function works:\n");
    fclaw_global_productionf("fcopy-2norm: %f\n", f->twoNorm());
    fclaw_global_productionf("fcopy-infnorm: %f\n\n", f->infNorm());
}
}
