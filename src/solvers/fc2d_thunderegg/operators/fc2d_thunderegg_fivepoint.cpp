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

#include "operators/fc2d_thunderegg_fivepoint.h"

#include "fc2d_thunderegg.h"
#include "fc2d_thunderegg_options.h"
#include "fc2d_thunderegg_vector.hpp"

#include <fclaw_elliptic_solver.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_clawpatch_output_vtk.h>

#include <fclaw_global.h>
#include <fclaw2d_domain.h>
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

class fivePoint : public PatchOperator<2>
{
public:
    fivePoint(const Domain<2>&      domain,
              const GhostFiller<2>& ghost_filler);
    
    fivePoint* clone() const override{
        return new fivePoint(*this);
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


fivePoint::fivePoint(const Domain<2>&      domain,
                     const GhostFiller<2>& ghost_filler) : PatchOperator<2>(domain,ghost_filler)
{
    /* Nothing to construct yet */
}


void fivePoint::applySinglePatchWithInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                               const PatchView<const double, 2>& u,
                                                               const PatchView<double, 2>& f) const 
{
    int mfields = u.getEnd()[2] + 1;
    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];
    //if physical boundary
    if (pinfo.hasNbr(Side<2>::west())){
        auto ghosts = u.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++){
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(0,j,m);
            }
        }
    }
    if (pinfo.hasNbr(Side<2>::east())){
        auto ghosts = u.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++){
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(mx-1,j,m);
            }
        }
    }

    if (pinfo.hasNbr(Side<2>::south())){
        auto ghosts = u.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++){
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,0,m);
            }
        }
    }
    if (pinfo.hasNbr(Side<2>::north())){
        auto ghosts = u.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++){
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,my-1,m);
            }
        }
    }
    applySinglePatch(pinfo,u,f);
}
void fivePoint::applySinglePatch(const PatchInfo<2>& pinfo, 
                                 const PatchView<const double, 2>& u,
                                 const PatchView<double, 2>& f) const 
{
    //const cast since u ghost values have to be modified
    //ThunderEgg doesn't care if ghost values are modified, just don't modify the interior values.

    //fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);

    int mfields = u.getEnd()[2] + 1;
    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];

#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    
    double dx = pinfo.spacings[0];
    double dy = pinfo.spacings[1];


    //if physical boundary
    if (!pinfo.hasNbr(Side<2>::west())){
        auto ghosts = u.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++){
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(0,j,m);
            }
        }
    }
    if (!pinfo.hasNbr(Side<2>::east())){
        auto ghosts = u.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++){
            for(int j = 0; j < my; j++){
                ghosts(j,m) = -u(mx-1,j,m);
            }
        }
    }

    if (!pinfo.hasNbr(Side<2>::south())){
        auto ghosts = u.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++){
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,0,m);
            }
        }
    }
    if (!pinfo.hasNbr(Side<2>::north())){
        auto ghosts = u.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++){
            for(int i = 0; i < mx; i++){
                ghosts(i,m) = -u(i,my-1,m);
            }
        }
    }

    double dx2 = dx*dx;
    double dy2 = dy*dy;

#if 1
    /* Five-point Laplacian */
    for(int m = 0; m < mfields; m++)
        for(int j = 0; j < my; j++)
            for(int i = 0; i < mx; i++)
            {
                double uij = u(i,j,m);
                f(i,j,m) = (u(i+1,j,m) - 2*uij + u(i-1,j,m))/dx2 + 
                             (u(i,j+1,m) - 2*uij + u(i,j-1,m))/dy2;
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


void fivePoint::modifyRHSForInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                       const PatchView<const double,2>& u, 
                                                       const PatchView<double,2>& f) const 
{
#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    

    int mfields = u.getEnd()[2]+1;
    int mx = pinfo.ns[0]; 
    int my = pinfo.ns[1];

    double dx = pinfo.spacings[0];
    double dx2 = dx*dx;

    double dy = pinfo.spacings[1];
    double dy2 = dy*dy;

    for(int m = 0; m < mfields; m++)
    {
        for(int j = 0; j < my; j++)
        {
            /* bool hasNbr(Side<D> s) */
            if (pinfo.hasNbr(Side<2>::west()))
                f(0,j,m) += -(u(-1,j,m)+u(0,j,m))/dx2;

            if (pinfo.hasNbr(Side<2>::east()))
                f(mx-1,j,m) += -(u(mx-1,j,m)+u(mx,j,m))/dx2;
        }

        for(int i = 0; i < mx; i++)
        {
            if (pinfo.hasNbr(Side<2>::south()))
                f(i,0,m) += -(u(i,-1,m)+u(i,0,m))/dy2;

            if (pinfo.hasNbr(Side<2>::north()))
                f(i,my-1,m) += -(u(i,my-1,m)+u(i,my,m))/dy2;
        }
    }
}
 

void fc2d_thunderegg_fivepoint_solve(fclaw_global_t *glob) 
{
    // get needed options
    fclaw_clawpatch_options_t *clawpatch_opt =
    fclaw_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fc2d_thunderegg_options_t *mg_opt = fc2d_thunderegg_get_options(glob);
  
    GhostFillingType fill_type = GhostFillingType::Faces;
#if 0  
    fc2d_thunderegg_vtable_t *mg_vt = fc2d_thunderegg_vt(glob);
#endif  

    // create thunderegg vector for eqn 0
    Vector<2> f = fc2d_thunderegg_get_vector(glob,RHS);

    // get patch size
    array<int, 2> ns = {clawpatch_opt->mx, clawpatch_opt->my};
    int mbc = clawpatch_opt->mbc;

    // get p4est structure
    fclaw_domain_t *domain = glob->domain;
    p4est_wrap_t *wrap = (p4est_wrap_t *) fclaw_domain_get_2d_domain(domain)->pp;

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

    // ghost filler
    BiLinearGhostFiller ghost_filler(te_domain, fill_type);

    // patch operator
    fivePoint op(te_domain,ghost_filler);

    // set the patch solver
    Iterative::BiCGStab<2> p_bicg;
    p_bicg.setTolerance(mg_opt->patch_iter_tol);
    p_bicg.setMaxIterations(mg_opt->patch_iter_max_it);
    Iterative::BiCGStab<2> p_cg;
    p_cg.setTolerance(mg_opt->patch_iter_tol);
    p_cg.setMaxIterations(mg_opt->patch_iter_max_it);
    shared_ptr<PatchSolver<2>>  solver;

    bitset<4> neumann_bitset;
    for(int i=0;i<4;i++){
        neumann_bitset[i]= mg_opt->boundary_conditions[i] == 2;
    }

    switch (mg_opt->patch_solver)
    {
        case BICG:
            solver = make_shared<Iterative::PatchSolver<2>>(p_bicg, op);
            break;
        case CG:
            solver = make_shared<Iterative::PatchSolver<2>>(p_cg, op);
            break;
#ifdef THUNDEREGG_FFTW_ENABLED
        case FFT:
            /* This ignores the five point operator defined above and just uses the 
               ThunderEgg operator 'Poisson'. */
            solver = make_shared<Poisson::FFTWPatchSolver<2>>(op,neumann_bitset);
            break;
#endif
        default:
            fclaw_global_essentialf("thunderegg_fivepoint : No valid patch solver specified\n");
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
        Domain<2> prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
            fivePoint patch_operator(curr_domain, ghost_filler);

            //smoother
            unique_ptr<GMG::Smoother<2>> smoother;
            switch (mg_opt->patch_solver)
            {
                case BICG:
                    smoother.reset(new Iterative::PatchSolver<2>(p_bicg, patch_operator));
                    break;
                case CG:
                    smoother.reset(new Iterative::PatchSolver<2>(p_cg, patch_operator));
                    break;
#ifdef THUNDEREGG_FFTW_ENABLED
                case FFT:
                    smoother.reset(new Poisson::FFTWPatchSolver<2>(patch_operator, neumann_bitset));
                    break;
#endif
                default:
                    fclaw_global_essentialf("thunderegg_fivepoint : No valid " \
                                            "patch solver specified\n");
                    exit(0);            
            }


            //restrictor
            GMG::LinearRestrictor<2> restrictor(curr_domain, 
                                                next_domain);

            //interpolator
            GMG::DirectInterpolator<2> interpolator(curr_domain, 
                                                    prev_domain);

            builder.addIntermediateLevel(patch_operator, *smoother, restrictor, 
                                         interpolator);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        BiLinearGhostFiller ghost_filler(curr_domain, fill_type);
        fivePoint patch_operator(curr_domain, ghost_filler);

        //smoot
        unique_ptr<GMG::Smoother<2>> smoother;
        switch (mg_opt->patch_solver)
        {
            case BICG:
                smoother.reset(new Iterative::PatchSolver<2>(p_bicg, patch_operator));
                break;
            case CG:
                smoother.reset(new Iterative::PatchSolver<2>(p_cg, patch_operator));
                break;
#ifdef THUNDEREGG_FFTW_ENABLED
            case FFT:
                smoother.reset(new Poisson::FFTWPatchSolver<2>(patch_operator, neumann_bitset));
                break;
#endif
            default:
                fclaw_global_essentialf("thunderegg_fivepoint : No valid " \
                                        "patch solver specified\n");
                exit(0);            
        }


        //interpolator
        GMG::DirectInterpolator<2> interpolator(curr_domain, prev_domain);

        builder.addCoarsestLevel(patch_operator, *smoother, interpolator);

        M = builder.getCycle();
    }

    // solve
    Vector<2> u = f.getZeroClone();

    Iterative::BiCGStab<2> iter_solver;
    iter_solver.setMaxIterations(mg_opt->max_it);
    iter_solver.setTolerance(mg_opt->tol);
    bool prt_output = mg_opt->verbosity_level > 0 && glob->mpirank == 0;
    int its = iter_solver.solve(op, u, f, M.get(),prt_output);

    fclaw_global_productionf("Iterations: %i\n", its);    

    /* Solution is copied to right hand side */
    fc2d_thunderegg_store_vector(glob, RHS, u);

}

