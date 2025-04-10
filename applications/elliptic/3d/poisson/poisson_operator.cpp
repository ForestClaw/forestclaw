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

#include "poisson_operator.h"

#include <fc3d_thunderegg.h>
#include <fc3d_thunderegg_options.h>
#include <fc3d_thunderegg_vector.hpp>

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_clawpatch_output_vtk.h>

#include <fclaw_global.h>
#include <fclaw_map.h>
#include <fclaw_map_brick.h>
#include <fclaw_options.h>
#include <fclaw_patch.h>
#include <fclaw_vtable.h>

#include <p8est_bits.h>
#include <p8est_wrap.h>

#include <forestclaw3d.h>

#include <ThunderEgg.h>

using namespace std;
using namespace ThunderEgg;

class fivePoint : public PatchOperator<3>
{
public:
    fivePoint(const Domain<3>&      domain,
              const GhostFiller<3>& ghost_filler);
    
    fivePoint* clone() const override{
        return new fivePoint(*this);
    }

    void applySinglePatch(const PatchInfo<3>& pinfo,
                          const PatchView<const double, 3>& u,
                          const PatchView<double, 3>& f) const override;

    void applySinglePatchWithInternalBoundaryConditions(const PatchInfo<3>& pinfo,
                                                        const PatchView<const double, 3>& u,
                                                        const PatchView<double, 3>& f) const override;

    void modifyRHSForInternalBoundaryConditions(const PatchInfo<3> &pinfo,
	                                            const PatchView<const double, 3> &u,
	                                            const PatchView<double, 3> &f) const override;



};


fivePoint::fivePoint(const Domain<3>&      domain,
                     const GhostFiller<3>& ghost_filler) : PatchOperator<3>(domain,ghost_filler)
{
    /* Nothing to construct yet */
}


void fivePoint::applySinglePatchWithInternalBoundaryConditions(const PatchInfo<3>& pinfo, 
                                                               const PatchView<const double, 3>& u,
                                                               const PatchView<double, 3>& f) const 
{
    applySinglePatch(pinfo,u,f);
}
void fivePoint::applySinglePatch(const PatchInfo<3>& pinfo, 
                                 const PatchView<const double, 3>& u,
                                 const PatchView<double, 3>& f) const 
{
    //const cast since u ghost values have to be modified
    //ThunderEgg doesn't care if ghost values are modified, just don't modify the interior values.

    //fc3d_thunderegg_options_t *mg_opt = fc3d_thunderegg_get_options(glob);

    int mfields = u.getEnd()[3] + 1;
    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];
    int mz = pinfo.ns[2];

#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];
    double dz = pinfo.spacings[2];


    //if physical boundary
    for(auto side : Side<3>::getValues())
    {
        if (!pinfo.hasNbr(side)){
            auto ghosts = u.getGhostSliceOn(side,{0});
            auto interior = u.getSliceOn(side,{0});
            Loop::OverInteriorIndexes<3>(ghosts, [&](auto &coord)
            {
                ghosts[coord] = -interior[coord];
            });
        }
    }

    double dx2 = dx*dx;
    double dy2 = dy*dy;
    double dz2 = dz*dz;

#if 1
    /* Five-point Laplacian */
    for(int m = 0; m < mfields; m++)
        for(int k = 0; k < mz; k++)
            for(int j = 0; j < my; j++)
                for(int i = 0; i < mx; i++)
                {
                    double uijk = u(i,j,k,m);
                    f(i,j,k,m) = (u(i+1,j,k,m) - 2*uijk + u(i-1,j,k,m))/dx2 + 
                                 (u(i,j+1,k,m) - 2*uijk + u(i,j-1,k,m))/dy2 +
                                 (u(i,j,k+1,m) - 2*uijk + u(i,j,k-1,m))/dz2;
                }
    
#else

    /* Five-point Laplacian : Slightly slower than then above.*/
    for(int j = 0; j < my; j++)
        for(int i = 0; i < mx; i++)
        {
            double uij = u(i,j);
            double flux[4];
            flux[0] = (uij - u(i-1,j));
            flux[1] = (u(i+1,j) - uij);
            flux[2] = (uij - u(i,j-1));
            flux[3] = (u(i,j+1) - uij);;
            f(i,j) = (flux[1]-flux[0])/dx2 + (flux[3] - flux[2])/dy2;
        }
#endif
    
}


void fivePoint::modifyRHSForInternalBoundaryConditions(const PatchInfo<3>& pinfo, 
                                                       const PatchView<const double,3>& u, 
                                                       const PatchView<double,3>& f) const 
{
}
 

void poisson_solve(fclaw_global_t *glob) 
{
    // get needed options
    fclaw_clawpatch_options_t *clawpatch_opt =
    fclaw_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fc3d_thunderegg_options_t *mg_opt = fc3d_thunderegg_get_options(glob);
  
    GhostFillingType fill_type = GhostFillingType::Faces;
#if 0  
    fc3d_thunderegg_vtable_t *mg_vt = fc3d_thunderegg_vt(glob);
#endif  

    // create thunderegg vector for eqn 0
    Vector<3> f = fc3d_thunderegg_get_vector(glob,RHS);

    // get patch size
    array<int, 3> ns = {clawpatch_opt->mx, clawpatch_opt->my, clawpatch_opt->mz};
    int mbc = clawpatch_opt->mbc;

    // get p4est structure
    fclaw_domain_t *domain = glob->domain;
    p8est_wrap_t *wrap = (p8est_wrap_t *)domain->d3->pp;

    // create map function
    P8estDomainGenerator::BlockMapFunc bmf = [&](int block_no, 
                                                 double unit_x, double unit_y, double unit_z, 
                                                 double &x, double &y, double& z) 
    {
        double x1,y1,z1;
        FCLAW_MAP_3D_BRICK2C(&glob->cont,&block_no,&unit_x, &unit_y, &unit_z, &x1, &y1, &z1);
        x = fclaw_opt->ax + (fclaw_opt->bx - fclaw_opt->ax) * x1;
        y = fclaw_opt->ay + (fclaw_opt->by - fclaw_opt->ay) * y1;
        z = fclaw_opt->az + (fclaw_opt->bz - fclaw_opt->az) * z1;
    };

    // generates levels of patches for GMG
    P8estDomainGenerator domain_gen(wrap->p4est, ns, mbc, bmf);

    // get finest level
    Domain<3> te_domain = domain_gen.getCoarserDomain();
    std::shared_ptr<Timer> timer = make_shared<Timer>(te_domain.getCommunicator());
    te_domain.setTimer(timer);

    // ghost filler
    TriLinearGhostFiller ghost_filler(te_domain, fill_type);

    // patch operator
    fivePoint op(te_domain,ghost_filler);

    // set the patch solver
    Iterative::BiCGStab<3> p_bicg;
    p_bicg.setTolerance(mg_opt->patch_iter_tol);
    p_bicg.setMaxIterations(mg_opt->patch_iter_max_it);
    Iterative::BiCGStab<3> p_cg;
    p_cg.setTolerance(mg_opt->patch_iter_tol);
    p_cg.setMaxIterations(mg_opt->patch_iter_max_it);
    shared_ptr<PatchSolver<3>>  solver;

    bitset<6> neumann_bitset;
    for(int i=0;i<6;i++){
        neumann_bitset[i]= mg_opt->boundary_conditions[i] == 2;
    }

    switch (mg_opt->patch_solver)
    {
        case BICG:
            solver = make_shared<Iterative::PatchSolver<3>>(p_bicg, op);
            break;
        case CG:
            solver = make_shared<Iterative::PatchSolver<3>>(p_cg, op);
            break;
#ifdef THUNDEREGG_FFTW_ENABLED
        case FFT:
            /* This ignores the five point operator defined above and just uses the 
               ThunderEgg operator 'Poisson'. */
            solver = make_shared<Poisson::FFTWPatchSolver<3>>(op,neumann_bitset);
            break;
#endif
        default:
            fclaw_global_essentialf("thunderegg_fivepoint : No valid patch solver specified\n");
            exit(0);            
    }

    // create gmg preconditioner
    shared_ptr<Operator<3>> M;

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
        GMG::CycleBuilder<3> builder(copts);
        
        //add finest level

        //next domain
        Domain<3> curr_domain = te_domain;
        Domain<3> next_domain = domain_gen.getCoarserDomain();
        next_domain.setTimer(timer);

        //restrictor
        GMG::LinearRestrictor<3> restrictor(curr_domain, 
                                            next_domain);

        builder.addFinestLevel(op, *solver, restrictor);

        //add intermediate levels
        Domain<3> prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();
            next_domain.setTimer(timer);

            //operator
            TriLinearGhostFiller ghost_filler(curr_domain, fill_type);
            fivePoint patch_operator(curr_domain, ghost_filler);

            //smoother
            unique_ptr<GMG::Smoother<3>> smoother;
            switch (mg_opt->patch_solver)
            {
                case BICG:
                    smoother.reset(new Iterative::PatchSolver<3>(p_bicg, patch_operator));
                    break;
                case CG:
                    smoother.reset(new Iterative::PatchSolver<3>(p_cg, patch_operator));
                    break;
#ifdef THUNDEREGG_FFTW_ENABLED
                case FFT:
                    smoother.reset(new Poisson::FFTWPatchSolver<3>(patch_operator, neumann_bitset));
                    break;
#endif
                default:
                    fclaw_global_essentialf("thunderegg_fivepoint : No valid " \
                                            "patch solver specified\n");
                    exit(0);            
            }


            //restrictor
            GMG::LinearRestrictor<3> restrictor(curr_domain, 
                                                next_domain);

            //interpolator
            GMG::DirectInterpolator<3> interpolator(curr_domain, 
                                                    prev_domain);

            builder.addIntermediateLevel(patch_operator, *smoother, restrictor, 
                                         interpolator);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        TriLinearGhostFiller ghost_filler(curr_domain, fill_type);
        fivePoint patch_operator(curr_domain, ghost_filler);

        //smoot
        unique_ptr<GMG::Smoother<3>> smoother;
        switch (mg_opt->patch_solver)
        {
            case BICG:
                smoother.reset(new Iterative::PatchSolver<3>(p_bicg, patch_operator));
                break;
            case CG:
                smoother.reset(new Iterative::PatchSolver<3>(p_cg, patch_operator));
                break;
#ifdef THUNDEREGG_FFTW_ENABLED
            case FFT:
                smoother.reset(new Poisson::FFTWPatchSolver<3>(patch_operator, neumann_bitset));
                break;
#endif
            default:
                fclaw_global_essentialf("thunderegg_fivepoint : No valid " \
                                        "patch solver specified\n");
                exit(0);            
        }


        //interpolator
        GMG::DirectInterpolator<3> interpolator(curr_domain, prev_domain);

        builder.addCoarsestLevel(patch_operator, *smoother, interpolator);

        M = builder.getCycle();
    }

    // solve
    Vector<3> u = f.getZeroClone();

    Iterative::BiCGStab<3> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);
    bool prt_output = mg_opt->verbosity_level > 0 && glob->mpirank == 0;
    int its = iter_solver.solve(op, u, f, M.get(),prt_output);

    fclaw_global_productionf("Iterations: %i\n", its);    
    //std::cout << *timer;

    /* Solution is copied to right hand side */
    fc3d_thunderegg_store_vector(glob, RHS, u);

}

