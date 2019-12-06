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

#include "fc2d_multigrid_solve.h"

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
#include <fclaw2d_map_brick.h>
#include <fclaw2d_options.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_vtable.h>

#include <p4est_bits.h>
#include <p4est_wrap.h>

#include <Thunderegg/BiCGStab.h>
#include <Thunderegg/BiCGStabPatchSolver.h>
#include <Thunderegg/VarPoisson/StarPatchOperator.h>
#include <Thunderegg/GMG/LinearRestrictor.h>
#include <Thunderegg/GMG/DrctIntp.h>
#include <Thunderegg/SchurDomainOp.h>
#include <Thunderegg/P4estDomGen.h>
#include <Thunderegg/GMG/CycleFactory.h>
#include <Thunderegg/BiLinearGhostFiller.h>

using namespace std;
using namespace Thunderegg;
using namespace Thunderegg::Schur;
using namespace Thunderegg::VarPoisson;

void fc2d_multigrid_solve(fclaw2d_global_t *glob) {
  // get needed options
  fclaw2d_clawpatch_options_t *clawpatch_opt =
      fclaw2d_clawpatch_get_options(glob);
  fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
  fc2d_multigrid_options_t *mg_opt = fc2d_multigrid_get_options(glob);
  fc2d_multigrid_vtable_t *mg_vt = fc2d_multigrid_vt();

  // this should be moved somewhere else
  PetscInitialize(nullptr, nullptr, nullptr, nullptr);

  // create thunderegg vector for eqn 0
  shared_ptr<Vector<2>> f(new fc2d_multigrid_vector(glob, 0));

  // get patch size
  array<int, 2> ns = {clawpatch_opt->mx, clawpatch_opt->my};

  // get p4est structure
  fclaw2d_domain_t *domain = glob->domain;
  p4est_wrap_t *wrap = (p4est_wrap_t *)domain->pp;

  // create map function
  P4estDomGen::BlockMapFunc bmf = [&](int block_no, double unit_x,      
                                      double unit_y, double &x, double &y) {
    double x1,y1,z1;
    FCLAW2D_MAP_BRICK2C(&glob->cont,&block_no,&unit_x, &unit_y, &x1, &y1, &z1);
    x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
    y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
  };

  // create neumann function
  IsNeumannFunc<2> inf = [&](Side<2> s, const array<double, 2> &lower,
                             const array<double, 2> &upper) {
    return mg_opt->boundary_conditions[s.toInt()] == 2;
  };

  // generates levels of patches for GMG
  shared_ptr<P4estDomGen> domain_gen(
      new P4estDomGen(wrap->p4est, ns, clawpatch_opt->mbc,inf, bmf));

  // get finest level
  shared_ptr<Domain<2>> te_domain = domain_gen->getFinestDomain();

  // define operators for problems

  // get beta function
  auto beta_func = [&](const std::array<double,2>& coord){
      double beta;
      double grad[2];
      mg_vt->fort_beta(&coord[0],&coord[1],&beta,grad);
      return beta;
  };

  // create vector for beta
  auto beta_vec = PetscVector<2>::GetNewVector(te_domain);
  DomainTools<2>::setValuesWithGhost(te_domain, beta_vec, beta_func);

  // ghost filler
  auto ghost_filler = make_shared<BiLinearGhostFiller>(te_domain);
  // patch operator
  auto op = make_shared<StarPatchOperator<2>>(beta_vec,te_domain,ghost_filler);

  // set the patch solver
  auto  solver = make_shared<BiCGStabPatchSolver<2>>(te_domain, ghost_filler, op,
                                                          mg_opt->patch_bcgs_tol,
                                                          mg_opt->patch_bcgs_max_it);

  // create matrix
  shared_ptr<Operator<2>> A = op;

  // create gmg preconditioner
  shared_ptr<Operator<2>> M;
  if(mg_opt->mg_prec){
    // options
    GMG::CycleOpts copts;
    copts.max_levels = mg_opt->max_levels;
    copts.patches_per_proc = mg_opt->patches_per_proc;
    copts.pre_sweeps = mg_opt->pre_sweeps;
    copts.post_sweeps = mg_opt->post_sweeps;
    copts.mid_sweeps = mg_opt->mid_sweeps;
    copts.coarse_sweeps = mg_opt->coarse_sweeps;
    copts.cycle_type = mg_opt->cycle_type;

    //generators for levels
    BiLinearGhostFiller::Generator filler_gen(ghost_filler);
    StarPatchOperator<2>::Generator   op_gen(op, filler_gen);
		BiCGStabPatchSolver<2>::Generator smooth_gen(solver, filler_gen, op_gen);
		GMG::LinearRestrictor<2>::Generator        restrictor_gen;
		GMG::DrctIntp<2>::Generator       interpolator_gen;
    M = GMG::CycleFactory<2>::getCycle(copts, domain_gen, restrictor_gen, interpolator_gen,
                                       smooth_gen, op_gen);
  }

  // solve
  shared_ptr<VectorGenerator<2>> vg(new DomainVG<2>(te_domain));
  shared_ptr<Vector<2>> u = vg->getNewVector();

  int its = BiCGStab<2>::solve(vg, A, u, f, M, mg_opt->max_it, mg_opt->tol);

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

