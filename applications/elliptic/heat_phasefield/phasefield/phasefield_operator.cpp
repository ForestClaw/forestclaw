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

#include "phasefield_operator.h"

#include "phasefield_patch_operator.h"

#include "phasefield_options.h"

#include <fc2d_thunderegg.h>
#include <fc2d_thunderegg_vector.hpp>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_clawpatch_output_vtk.h>

#include <fclaw_global.h>
#include <fclaw2d_map.h>
#include <fclaw2d_map_brick.h>
#include <fclaw_options.h>
#include <fclaw_patch.h>
#include <fclaw_vtable.h>

#include <p4est_bits.h>
#include <p4est_wrap.h>

#include <stdlib.h>

using namespace std;
using namespace ThunderEgg;

void phasefield_set_lambda(double lambda)
{
    phasefield::setLambda(lambda);

}

double phasefield_get_lambda()
{
    return phasefield::getLambda();
}
Vector<2> restrict_phi_n_vec(const Vector<2>& prev_beta_vec, 
                             const Domain<2>& prev_domain, 
                             const Domain<2>& curr_domain)
{
    GMG::LinearRestrictor<2> restrictor(prev_domain,curr_domain, true);
    return restrictor.restrict(prev_beta_vec);
}

void phasefield_solve(fclaw_global_t *glob) 
{
    // get needed options
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
  
    GhostFillingType fill_type = GhostFillingType::Faces;
#if 0  
    fc2d_thunderegg_vtable_t *mg_vt = fc2d_thunderegg_vt(glob);
#endif  

    // create thunderegg vector for eqn 0
    Vector<2> f = fc2d_thunderegg_get_vector(glob,RHS);

    // get patch size
    array<int, 2> ns = {clawpatch_opt->mx, clawpatch_opt->my};

    // get p4est structure
    fclaw_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *)domain->d2->pp;

    // create map function
    fclaw2d_map_context_t *cont = fclaw2d_map_get(glob);
    P4estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x,      
                                        double unit_y, double &x, double &y) 
    {
        double x1,y1,z1;
        FCLAW2D_MAP_BRICK2C(&cont,&block_no,&unit_x, &unit_y, &x1, &y1, &z1);
        x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
        y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
    };

    // generates levels of patches for GMG;  2 layers of ghost cells    
    P4estDomainGenerator domain_gen(wrap->p4est, ns, clawpatch_opt->mbc, bmf);

    // get finest level
    Domain<2> te_domain = domain_gen.getFinestDomain();

    /* Store phi at time level n */
    Vector<2> beta_vec = fc2d_thunderegg_get_vector(glob,STORE_STATE);    

    // ghost filler
    BiLinearGhostFiller ghost_filler(te_domain, fill_type);

    // patch operator
    phasefield op(glob,beta_vec,te_domain,ghost_filler);

    // set the patch solver
    Iterative::CG<2> patch_cg;
    patch_cg.setTolerance(mg_opt->patch_iter_tol);
    patch_cg.setMaxIterations(mg_opt->patch_iter_max_it);
    Iterative::BiCGStab<2> patch_bicg;
    patch_bicg.setTolerance(mg_opt->patch_iter_tol);
    patch_bicg.setMaxIterations(mg_opt->patch_iter_max_it);

    Iterative::Solver<2>* patch_iterative_solver = nullptr;
    switch(mg_opt->patch_solver){
        case CG:
            patch_iterative_solver = &patch_cg;
            break; 
        case BICG:
            patch_iterative_solver = &patch_bicg;
            break;
        default:
            fclaw_global_essentialf("phasefield : No valid " \
                                    "patch solver specified\n");
            exit(0);            
    }

    Iterative::PatchSolver<2> solver(*patch_iterative_solver,op,true);

    // create gmg preconditioner
    shared_ptr<Operator<2>> M;

    if (mg_opt->mg_prec && domain_gen.hasCoarserDomain())
    {
        // options
        GMG::CycleOpts copts;
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

        //restrictor
        GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                            next_domain); 

        builder.addFinestLevel(op, solver, restrictor);

        //add intermediate levels
        auto prev_phi_n_vec = beta_vec;
        auto prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
            Vector<2> restricted_phi_n_vec = restrict_phi_n_vec(prev_phi_n_vec, 
                                                                prev_domain, curr_domain);
            phasefield patch_operator(glob,restricted_phi_n_vec,curr_domain, 
                                      ghost_filler);
            prev_phi_n_vec = restricted_phi_n_vec;

            //smoother
            Iterative::PatchSolver<2> smoother(*patch_iterative_solver, patch_operator, true);

            //restrictor
            GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                                next_domain);

            //interpolator
            GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                    prev_domain);

            builder.addIntermediateLevel(patch_operator, smoother, restrictor, 
                                         interpolator);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
        Vector<2> restricted_phi_n_vec = restrict_phi_n_vec(prev_phi_n_vec, 
                                                            prev_domain, curr_domain);
        phasefield coarse_patch_operator(glob,restricted_phi_n_vec, curr_domain, 
                                         ghost_filler);

        //smoother
        Iterative::PatchSolver<2> smoother(*patch_iterative_solver, coarse_patch_operator, true);
        //interpolator
        GMG::DirectInterpolator<2> interpolator(curr_domain, prev_domain);

        builder.addCoarsestLevel(coarse_patch_operator, smoother, interpolator);

        M = builder.getCycle();
    }

#if 0
    // Set starting conditions
    Vector<2> u = fc2d_thunderegg_get_vector(glob,SOLN);
#else
    Vector<2> u = f.getZeroClone();
#endif    


    Iterative::BiCGStab<2> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);    
    bool prt_output = glob->mpirank == 0;
    int its = iter_solver.solve(op, u, f, M.get(), prt_output);

    fclaw_global_productionf("Iterations: %i\n", its);    

    /* Solution is copied to right hand side */
    fc2d_thunderegg_store_vector(glob, RHS, u);
}

