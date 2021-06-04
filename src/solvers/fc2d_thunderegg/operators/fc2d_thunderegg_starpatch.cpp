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

#include "fc2d_thunderegg_starpatch.h"

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"
#include "fc2d_thunderegg_vector.hpp"

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

#include <ThunderEgg.h>

using namespace std;
using namespace ThunderEgg;


/**
 * @brief Restrict the beta coefficient vector for the new domain
 * 
 * @param prev_beta_vec the finer beta vector
 * @param prev_domain the previous (finer) domain
 * @param curr_domain the current (coarser) domain
 * @return shared_ptr<Vector<2>> the restricted beta vector
 */
shared_ptr<ValVector<2>> restrict_beta_vec(shared_ptr<Vector<2>> prev_beta_vec, 
                                        shared_ptr<Domain<2>> prev_domain, 
                                        shared_ptr<Domain<2>> curr_domain)
{
    GMG::LinearRestrictor<2> restrictor(prev_domain,curr_domain, prev_beta_vec->getNumComponents(), true);
    auto new_beta_vec = ValVector<2>::GetNewVector(curr_domain, prev_beta_vec->getNumComponents());
    restrictor.restrict(prev_beta_vec, new_beta_vec);
    return new_beta_vec;
}

void fc2d_thunderegg_starpatch_solve(fclaw2d_global_t *glob) 
{
    // get needed options
    fclaw2d_clawpatch_options_t *clawpatch_opt =
    fclaw2d_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);
  
#if 0  
    fc2d_thunderegg_vtable_t *mg_vt = fc2d_thunderegg_vt();
#endif  

    // create thunderegg vector for eqn 0
    shared_ptr<Vector<2>> f = make_shared<fc2d_thunderegg_vector>(glob,RHS);

    // get patch size
    array<int, 2> ns = {clawpatch_opt->mx, clawpatch_opt->my};

    // get p4est structure
    fclaw2d_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *)domain->pp;

    // create map function
    P4estDomGen::BlockMapFunc bmf = [&](int block_no, double unit_x,      
                                        double unit_y, double &x, double &y) 
    {
        double x1,y1,z1;
        FCLAW2D_MAP_BRICK2C(&glob->cont,&block_no,&unit_x, &unit_y, &x1, &y1, &z1);
        x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
        y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
    };

    // create neumann function
    IsNeumannFunc<2> inf = [&](Side<2> s, const array<double, 2> &lower,
                               const array<double, 2> &upper) 
    {
        return mg_opt->boundary_conditions[s.getIndex()] == 2;
    };

    // generates levels of patches for GMG
    P4estDomGen domain_gen(wrap->p4est, ns, 1, inf, bmf);

    // get finest level
    shared_ptr<Domain<2>> te_domain = domain_gen.getFinestDomain();

    // define operators for problems

#if 1
    // get beta function
    auto beta_func = [&](const std::array<double,2>& coord){
        double x = coord[0];
        double y = coord[1];
        double beta = 1 + 0*x*y;
        return beta;
    };
#endif  

    // create vector for beta
    auto beta_vec = ValVector<2>::GetNewVector(te_domain, clawpatch_opt->rhs_fields);
    DomainTools::SetValuesWithGhost<2>(te_domain, beta_vec, beta_func);

    // ghost filler
    auto ghost_filler = make_shared<BiLinearGhostFiller>(te_domain);

    // patch operator
    auto op = make_shared<VarPoisson::StarPatchOperator<2>>(beta_vec,te_domain,ghost_filler);

    // set the patch solver
    auto p_bcgs = make_shared<Iterative::BiCGStab<2>>();
    p_bcgs->setTolerance(mg_opt->patch_bcgs_tol);
    p_bcgs->setMaxIterations(mg_opt->patch_bcgs_max_it);
    shared_ptr<PatchSolver<2>>  solver;
    if(strcmp(mg_opt->patch_solver_type , "BCGS") == 0){
        solver = make_shared<Iterative::PatchSolver<2>>(p_bcgs, op);
    }else if(strcmp(mg_opt->patch_solver_type , "FFT") == 0){
        solver = make_shared<Poisson::FFTWPatchSolver<2>>(op);
    }

    // create matrix
    shared_ptr<Operator<2>> A = op;

    // create gmg preconditioner
    shared_ptr<Operator<2>> M;

    if(mg_opt->mg_prec && domain_gen.hasCoarserDomain())
    {
        // options
        GMG::CycleOpts copts;
        copts.max_levels = mg_opt->max_levels;
        copts.patches_per_proc = mg_opt->patches_per_proc;
        copts.pre_sweeps = mg_opt->pre_sweeps;
        copts.post_sweeps = mg_opt->post_sweeps;
        copts.mid_sweeps = mg_opt->mid_sweeps;
        copts.coarse_sweeps = mg_opt->coarse_sweeps;
        copts.cycle_type = mg_opt->cycle_type;

        //GMG cycle builder
        GMG::CycleBuilder<2> builder(copts);
        
        //add finest level

        //next domain
        auto curr_domain = te_domain;
        auto next_domain = domain_gen.getCoarserDomain();

        //operator
        auto patch_operator = op;

        //smoother
        shared_ptr<GMG::Smoother<2>> smoother = solver;

        //restrictor
        auto restrictor = make_shared<GMG::LinearRestrictor<2>>(curr_domain, 
                                                                next_domain, clawpatch_opt->rhs_fields);

        //vector generator
        auto vg = make_shared<ValVectorGenerator<2>>(curr_domain, clawpatch_opt->rhs_fields);

        builder.addFinestLevel(patch_operator, smoother, restrictor, vg);

        //add intermediate levels
        auto prev_beta_vec = beta_vec;
        auto prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            auto ghost_filler = make_shared<BiLinearGhostFiller>(curr_domain);
            auto restricted_beta_vec = restrict_beta_vec(prev_beta_vec, prev_domain, curr_domain);
            patch_operator = make_shared<VarPoisson::StarPatchOperator<2>>(restricted_beta_vec, curr_domain, ghost_filler);
            prev_beta_vec = restricted_beta_vec;

            //smoother
            shared_ptr<GMG::Smoother<2>> smoother;
            if(strcmp(mg_opt->patch_solver_type , "BCGS") == 0){
                smoother = make_shared<Iterative::PatchSolver<2>>(p_bcgs, patch_operator);
            }else if(strcmp(mg_opt->patch_solver_type , "FFT") == 0){
                smoother = make_shared<Poisson::FFTWPatchSolver<2>>(patch_operator);
            }

            //restrictor
            auto restrictor = make_shared<GMG::LinearRestrictor<2>>(curr_domain, 
                                                                    next_domain, 
                                                                    clawpatch_opt->rhs_fields);

            //interpolator
            auto interpolator = make_shared<GMG::DirectInterpolator<2>>(curr_domain, 
                                                                        prev_domain, 
                                                                        clawpatch_opt->rhs_fields);

            //vector generator
            vg = make_shared<ValVectorGenerator<2>>(curr_domain, clawpatch_opt->rhs_fields);

            builder.addIntermediateLevel(patch_operator, smoother, restrictor, interpolator, vg);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        auto ghost_filler = make_shared<BiLinearGhostFiller>(curr_domain);
        auto restricted_beta_vec = restrict_beta_vec(prev_beta_vec, prev_domain, curr_domain);
        patch_operator = make_shared<VarPoisson::StarPatchOperator<2>>(restricted_beta_vec, curr_domain, ghost_filler);

        //smoother
        if(strcmp(mg_opt->patch_solver_type , "BCGS") == 0){
            smoother = make_shared<Iterative::PatchSolver<2>>(p_bcgs, patch_operator);
        }else if(strcmp(mg_opt->patch_solver_type , "FFT") == 0){
            smoother = make_shared<Poisson::FFTWPatchSolver<2>>(patch_operator);
        }

        //interpolator
        auto interpolator = make_shared<GMG::DirectInterpolator<2>>(curr_domain, 
                                                                    prev_domain, 
                                                                    clawpatch_opt->rhs_fields);

        //vector generator
        vg = make_shared<ValVectorGenerator<2>>(curr_domain, clawpatch_opt->rhs_fields);

        builder.addCoarsestLevel(patch_operator, smoother, interpolator, vg);

        M = builder.getCycle();
    }

    // solve
    auto vg = make_shared<ValVectorGenerator<2>>(te_domain, clawpatch_opt->rhs_fields);
    shared_ptr<Vector<2>> u = vg->getNewVector();

    Iterative::BiCGStab<2> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);
    bool vl = mg_opt->verbosity_level != 0;
    int its = iter_solver.solve(vg, A, u, f, M, vl);

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

