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

#include "operators/fc2d_thunderegg_varpoisson.h"

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"
#include "fc2d_thunderegg_vector.hpp"

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

#include <ThunderEgg.h>

using namespace std;
using namespace ThunderEgg;

#if 0
double varpoisson_beta_coeff(const std::array<double,2>& coord)
{
    //double grad[2];
    return 1.0;
}  
#endif


class varpoisson : public PatchOperator<2>
{
    protected:
    Vector<2> beta;

    public:
    varpoisson(const Vector<2>& beta_vec, 
               const Domain<2>& domain,
               const GhostFiller<2>& ghost_filler);

    varpoisson* clone() const override
    {
        return new varpoisson(*this);
    }
    void applySinglePatch(const PatchInfo<2>& pinfo,
                          const PatchView<const double, 2>& u,
                          const PatchView<double, 2>& f) const override;

    void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<2>& pinfo,
                                                        const PatchView<const double, 2>& u,
                                                        const PatchView<double, 2>& f) const override;

    void modifyRHSForInternalBoundaryConditions(const PatchInfo<2> &pinfo,
	                                            const PatchView<const double, 2> &u,
	                                            const PatchView<double, 2> &f) const override;


};


varpoisson::varpoisson(const Vector<2>& coeffs,
                       const Domain<2>& domain,
                       const GhostFiller<2>& ghost_filler) : 
                            PatchOperator<2>(domain,ghost_filler),
                            beta(coeffs)
{
    ghost_filler.fillGhost(this->beta);
}


void varpoisson::applySinglePatchWithInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                                const PatchView<const double, 2>& u,
                                                                const PatchView<double, 2>& f) const 
{
    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];

    int mfields = u.getEnd()[2]+1;

    if (pinfo.hasNbr(Side<2>::west()))
    {
        auto ghosts = u.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(0,j,m);
            }
        }
    }
    if (pinfo.hasNbr(Side<2>::east())){
        auto ghosts = u.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(mx-1,j,m);
            }
        }
    }

    if (pinfo.hasNbr(Side<2>::south())){
        auto ghosts = u.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,0,m);
            }
        }
    }
    if (pinfo.hasNbr(Side<2>::north())){
        auto ghosts = u.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,my-1,m);
            }
        }
    }

    applySinglePatch(pinfo, u, f);
}
void varpoisson::applySinglePatch(const PatchInfo<2>& pinfo, 
                                  const PatchView<const double, 2>& u,
                                  const PatchView<double, 2>& f) const 
{
    const ComponentView<const double,2>  b  = beta.getComponentView(0, pinfo.local_index);
    //LocalData<2>& b = const_cast<LocalData<2>&>(beta[0]);

    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];

    int mfields = u.getEnd()[2]+1;

#if 0    
    int mbc = pinfo.num_ghost_cells;
    double xlower = pinfo.starts[0];
    double ylower = pinfo.starts[1];
#endif    
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];

    //if physical boundary
    if (!pinfo.hasNbr(Side<2>::west()))
    {
        auto ghosts = u.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(0,j,m);
            }
        }
    }
    if (!pinfo.hasNbr(Side<2>::east())){
        auto ghosts = u.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(mx-1,j,m);
            }
        }
    }

    if (!pinfo.hasNbr(Side<2>::south())){
        auto ghosts = u.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,0,m);
            }
        }
    }
    if (!pinfo.hasNbr(Side<2>::north())){
        auto ghosts = u.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++)
        {
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,my-1,m);
            }
        }
    }

#if 1
    double dx2 = 2*dx*dx;
    double dy2 = 2*dy*dy;
    for(int m = 0; m < mfields; m++)
        for(int j = 0; j < my; j++)
            for(int i = 0; i < mx; i++)
            {
                double b0 = (b(i,j)   + b(i-1,j));
                double b1 = (b(i+1,j) + b(i,j));
                double b2 = (b(i,j)   + b(i,j-1));
                double b3 = (b(i,j+1) + b(i,j));

                double flux[4];
                flux[0] = b0*(u(i,j,m) - u(i-1,j,m));
                flux[1] = b1*(u(i+1,j,m) - u(i,j,m));
                flux[2] = b2*(u(i,j,m) - u(i,j-1,m));
                flux[3] = b3*(u(i,j+1,m) - u(i,j,m));

                f(i,j,m) = (flux[1]-flux[0])/dx2 + (flux[3] - flux[2])/dy2;

            }
#endif

#if 0
        double dx2 = 2*dx*dx;
        double dy2 = 2*dy*dy;
        /* This is considerably slower, despite the fact that it
            computes half as many fluxes.  5.95s vs. 4.85s */
        for(int j = 0; j < my+1; j++)
            for(int i = 0; i < mx+1; i++)
            {
                double b0 = (b(i,j)   + b(i-1,j));
                double flux0 = b0*(u(i,j) - u(i-1,j));

                double b2 = (b(i,j) + b(i,j-1));
                double flux2 = b2*(u(i,j) - u(i,j-1));

                if (i < mx)
                    f(i,j) = -flux0/dx2;

                if (i > 0)
                    f(i-1,j) += flux0/dx2;

                if (j < my)
                    f(i,j) += -flux2/dy2;

                if (j > 0)
                    f(i,j-1) += flux2/dy2;

            }
#endif            

}


void varpoisson::modifyRHSForInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                        const PatchView<const double,2>& u, 
                                                        const PatchView<double,2>& f) const 
{

    const ComponentView<const double,2> b = beta.getComponentView(0, pinfo.local_index);

    int mfields = u.getEnd()[2]+1;

    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];
#if 0    
    int mbc = pinfo.num_ghost_cells;
    double xlower = pinfo.starts[0];
    double ylower = pinfo.starts[1];
#endif    
    double dx = pinfo.spacings[0];
    double dx2 = dx*dx;

    double dy = pinfo.spacings[1];
    double dy2 = dy*dy;

    for(int m = 0; m < mfields; m++)
    {
        for(int j = 0; j < my; j++)
        {
            double b0 = (b(0,j) + b(-1,j))/2;
            double b1 = (b(mx,j) + b(mx-1,j))/2;
            /* bool hasNbr(Side<2> s) */
            if (pinfo.hasNbr(Side<2>::west()))
                f(0,j,m) += -b0*(u(-1,j,m)+u(0,j,m))/dx2;

            if (pinfo.hasNbr(Side<2>::east()))
                f(mx-1,j,m) += -b1*(u(mx-1,j,m)+u(mx,j,m))/dx2;
        }

        for(int i = 0; i < mx; i++)
        {
            double b2 = (b(i,0) + b(i,-1))/2;
            double b3 = (b(i,my) + b(i,my-1))/2;
            if (pinfo.hasNbr(Side<2>::south()))
                f(i,0,m) += -b2*(u(i,-1,m)+u(i,0,m))/dy2;

            if (pinfo.hasNbr(Side<2>::north()))
                f(i,my-1,m) += -b3*(u(i,my-1,m)+u(i,my,m))/dy2;
        }
    }
}
 

Vector<2> varpoisson_restrict_beta_vec(const Vector<2>& prev_beta_vec, 
                                       const Domain<2>& prev_domain, 
                                       const Domain<2>& curr_domain)
{
    GMG::LinearRestrictor<2> restrictor(prev_domain,curr_domain, true);
    return restrictor.restrict(prev_beta_vec);
}



/* Public interface - this function is virtualized */

void fc2d_thunderegg_varpoisson_solve(fclaw_global_t *glob) 
{
    // get needed options
    fclaw_clawpatch_options_t *clawpatch_opt =
                        fclaw_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    GhostFillingType fill_type = GhostFillingType::Faces;
  
    // create thunderegg vector for eqn 0
    Vector<2> f = fc2d_thunderegg_get_vector(glob,RHS);

    // get patch size
    array<int, 2> ns = {clawpatch_opt->d2->mx, clawpatch_opt->d2->my};
    int mbc = clawpatch_opt->mbc;

    // get p4est structure
    fclaw_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *)domain->pp;

    // create map function
    fclaw2d_map_context_t* cont = fclaw2d_map_get(glob);
    P4estDomainGenerator::BlockMapFunc bmf = [&](int block_no, double unit_x,      
                                        double unit_y, double &x, double &y) 
    {
        double x1,y1,z1;
        FCLAW2D_MAP_BRICK2C(&cont,&block_no,&unit_x, &unit_y, &x1, &y1, &z1);
        x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
        y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
    };

    // generates levels of patches for GMG
    P4estDomainGenerator domain_gen(wrap->p4est, ns, mbc, bmf);

    // get finest level
    Domain<2> te_domain = domain_gen.getFinestDomain();

    // define operators for problems


    fc2d_thunderegg_vtable_t*  mg_vt = fc2d_thunderegg_vt(glob);

#if 1
    // get beta function
    auto beta_func = [&](const std::array<double,2>& coord){
        double beta = 0.0;
        double grad[2];
        mg_vt->fort_beta(&coord[0],&coord[1],&beta,grad);
        //beta = 1;
        return beta;
    };

    // create vector for beta

    //auto beta_vec = ValVector<2>::GetNewVector(te_domain, clawpatch_opt->rhs_fields);
    Vector<2> beta_vec = f.getZeroClone();
    DomainTools::SetValuesWithGhost<2>(te_domain, beta_vec, beta_func);
#endif

    // ghost filler
    BiLinearGhostFiller ghost_filler(te_domain, fill_type);

    // patch operator
    varpoisson op(beta_vec,te_domain,ghost_filler);

    // set the patch solver
    Iterative::CG<2> patch_cg;
    patch_cg.setTolerance(mg_opt->patch_iter_tol);
    patch_cg.setTolerance(mg_opt->patch_iter_max_it);
    Iterative::BiCGStab<2> patch_bicg;
    patch_bicg.setTolerance(mg_opt->patch_iter_tol);
    patch_bicg.setTolerance(mg_opt->patch_iter_max_it);

    Iterative::Solver<2>* patch_iterative_solver = nullptr;
    switch(mg_opt->patch_solver){
        case CG:
            patch_iterative_solver = &patch_cg;
            break; 
        case BICG:
            patch_iterative_solver = &patch_bicg;
            break;
        default:
            fclaw_global_essentialf("thunderegg_varpoisson : No valid " \
                                    "patch solver specified\n");
            exit(0);            
    }
    Iterative::PatchSolver<2> solver(*patch_iterative_solver, op);

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

        builder.addFinestLevel(op, solver, restrictor);

        //add intermediate levels
        Vector<2> prev_beta_vec = beta_vec;
        Domain<2> prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
            Vector<2> restricted_beta_vec = varpoisson_restrict_beta_vec(prev_beta_vec, 
                                                                         prev_domain, curr_domain);
            varpoisson patch_operator(restricted_beta_vec, 
                                      curr_domain, ghost_filler);
            prev_beta_vec = restricted_beta_vec;

            //smoother
           Iterative::PatchSolver<2> smoother(*patch_iterative_solver, patch_operator);

            //restrictor
            GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                                next_domain);

            //interpolator
            GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                    prev_domain);

            builder.addIntermediateLevel(patch_operator, smoother, restrictor, interpolator);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
        Vector<2> restricted_beta_vec = varpoisson_restrict_beta_vec(prev_beta_vec, prev_domain, curr_domain);
        varpoisson patch_operator(restricted_beta_vec, curr_domain, 
                                  ghost_filler);

        //smoother
        Iterative::PatchSolver<2> smoother(*patch_iterative_solver, patch_operator);

        //interpolator
        GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                prev_domain);

        builder.addCoarsestLevel(patch_operator, smoother, interpolator);

        M = builder.getCycle();
    }

    // solve
    Vector<2> u = f.getZeroClone();

    Iterative::BiCGStab<2> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);

    bool vl = mg_opt->verbosity_level > 0 && glob->mpirank == 0;
    int its = iter_solver.solve(op, u, f, M.get(),vl);

    // copy solution into rhs
    fc2d_thunderegg_store_vector(glob, RHS, u);
    fclaw_global_productionf("Iterations: %i\n", its);    
}

