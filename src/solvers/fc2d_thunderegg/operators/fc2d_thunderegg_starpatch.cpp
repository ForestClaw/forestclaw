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

#include "fc2d_thunderegg_starpatch.h"

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"
#include "fc2d_thunderegg_vector.hpp"

#include <fclaw2d_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_clawpatch_output_vtk.h>

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
Vector<2> restrict_beta_vec(const Vector<2>& prev_beta_vec, 
                            const Domain<2>& prev_domain, 
                            const Domain<2>& curr_domain)
{
    GMG::LinearRestrictor<2> restrictor(prev_domain, curr_domain,  true);
    return restrictor.restrict(prev_beta_vec);
}

void fc2d_thunderegg_starpatch_solve(fclaw2d_global_t *glob) 
{
    // get needed options
    fclaw_clawpatch_options_t *clawpatch_opt =
    fclaw_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    // ghost filling type
    GhostFillingType fill_type = GhostFillingType::Faces;
  
#if 0  
    fc2d_thunderegg_vtable_t *mg_vt = fc2d_thunderegg_vt(glob);
#endif  

    // create thunderegg vector for eqn 0
    Vector<2> f = fc2d_thunderegg_get_vector(glob,RHS);

    // get patch size
    array<int, 2> ns = {clawpatch_opt->d2->mx, clawpatch_opt->d2->my};
    int mbc = clawpatch_opt->mbc;

    // get p4est structure
    fclaw_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *)domain->pp;

    // create map function
    fclaw2d_map_context_t* cont = fclaw2d_global_get_map(glob);
    P4estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x,      
                                        double unit_y, double &x, double &y) 
    {
        double x1,y1,z1;
        FCLAW2D_MAP_BRICK2C(&cont,&block_no,&unit_x, &unit_y, &x1, &y1, &z1);
        x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
        y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
    };

    bitset<4> neumann_bitset;
    for(int i=0;i<4;i++){
        neumann_bitset[i]= mg_opt->boundary_conditions[i] == 2;
    }

    // generates levels of patches for GMG
    P4estDomainGenerator domain_gen(wrap->p4est, ns, mbc, bmf);

    // get finest level
    Domain<2> te_domain = domain_gen.getFinestDomain();

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
    Vector<2> beta_vec = f.getZeroClone();
    DomainTools::SetValuesWithGhost<2>(te_domain, beta_vec, beta_func);

    // ghost filler
    BiLinearGhostFiller ghost_filler(te_domain, fill_type);

    // patch operator
    VarPoisson::StarPatchOperator op(beta_vec, te_domain, ghost_filler);

    // set the patch solver
    Iterative::BiCGStab<2> p_cg;
    p_cg.setTolerance(mg_opt->patch_iter_tol);
    p_cg.setMaxIterations(mg_opt->patch_iter_max_it);
    Iterative::BiCGStab<2> p_bicg;
    p_bicg.setTolerance(mg_opt->patch_iter_tol);
    p_bicg.setMaxIterations(mg_opt->patch_iter_max_it);

    unique_ptr<PatchSolver<2>>  solver;
    switch(mg_opt->patch_solver){
        case CG:
            solver.reset(new Iterative::PatchSolver<2>(p_cg, op));
            break;
        case BICG:
            solver.reset(new Iterative::PatchSolver<2>(p_bicg, op));
            break;
#ifdef THUNDEREGG_FFTW_ENABLED
        case FFT:
            solver.reset(new Poisson::FFTWPatchSolver<2>(op, neumann_bitset));
            break;
#endif
        default:
            fclaw_global_essentialf("thunderegg_starpatch : No valid " \
                                    "patch solver specified\n");
            exit(0);            
    }

    // create gmg preconditioner
    shared_ptr<Operator<2>> M;

    if(mg_opt->mg_prec && domain_gen.hasCoarserDomain())
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
        Domain<2> curr_domain = te_domain;
        Domain<2> next_domain = domain_gen.getCoarserDomain();

        //restrictor
        GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                            next_domain);

        builder.addFinestLevel(op, *solver, restrictor);

        //add intermediate levels
        auto prev_beta_vec = beta_vec;
        Domain<2> prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
            Vector<2> restricted_beta_vec = restrict_beta_vec(prev_beta_vec, prev_domain, curr_domain);
            VarPoisson::StarPatchOperator<2> patch_operator(restricted_beta_vec, curr_domain, ghost_filler);
            prev_beta_vec = restricted_beta_vec;

            //smoother
            unique_ptr<GMG::Smoother<2>> smoother;
            switch(mg_opt->patch_solver){
                case CG:
                    smoother.reset(new Iterative::PatchSolver<2>(p_cg, patch_operator));
                    break;
                case BICG:
                    smoother.reset(new Iterative::PatchSolver<2>(p_bicg, patch_operator));
                    break;
#ifdef THUNDEREGG_FFTW_ENABLED
                case FFT:
                    smoother.reset(new Poisson::FFTWPatchSolver<2>(patch_operator, neumann_bitset));
                    break;
#endif
                default:
                    fclaw_global_essentialf("thunderegg_starpatch : No valid " \
                                            "patch solver specified\n");
                    exit(0);            
            }

            //restrictor
            GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                                next_domain);

            //interpolator
            GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                    prev_domain);

            builder.addIntermediateLevel(patch_operator, *smoother, restrictor, interpolator);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
        Vector<2> restricted_beta_vec = restrict_beta_vec(prev_beta_vec, prev_domain, curr_domain);
        VarPoisson::StarPatchOperator<2> patch_operator(restricted_beta_vec, curr_domain, ghost_filler);

        //smoother
        unique_ptr<GMG::Smoother<2>> smoother;
        switch(mg_opt->patch_solver){
            case CG:
                smoother.reset(new Iterative::PatchSolver<2>(p_cg, patch_operator));
                break;
            case BICG:
                smoother.reset(new Iterative::PatchSolver<2>(p_bicg, patch_operator));
                break;
#ifdef THUNDEREGG_FFTW_ENABLED
            case FFT:
                smoother.reset(new Poisson::FFTWPatchSolver<2>(patch_operator, neumann_bitset));
                break;
#endif
            default:
                fclaw_global_essentialf("thunderegg_starpatch : No valid " \
                                        "patch solver specified\n");
                exit(0);            
        }



        //interpolator
        GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                prev_domain);

        //vector generator

        builder.addCoarsestLevel(patch_operator, *smoother, interpolator);

        M = builder.getCycle();
    }

    // solve
    Vector<2> u = f.getZeroClone();

    Iterative::BiCGStab<2> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);
    bool vl = mg_opt->verbosity_level > 0 && glob->mpirank == 0;
    int its = iter_solver.solve(op, u, f, M.get(), vl);

    fclaw_global_productionf("Iterations: %i\n", its);

    // copy solution into rhs
    fc2d_thunderegg_store_vector(glob, RHS, u);
}

