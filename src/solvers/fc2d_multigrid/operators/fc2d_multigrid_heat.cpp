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

#include "operators/fc2d_multigrid_heat.h"

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

#include <ThunderEgg/BiCGStab.h>
#include <ThunderEgg/BiCGStabPatchSolver.h>
#include <ThunderEgg/VarPoisson/StarPatchOperator.h>
#include <ThunderEgg/Poisson/FFTWPatchSolver.h>
#include <ThunderEgg/GMG/LinearRestrictor.h>
#include <ThunderEgg/GMG/DirectInterpolator.h>
#include <ThunderEgg/P4estDomGen.h>
#include <ThunderEgg/GMG/CycleBuilder.h>
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/ValVectorGenerator.h>

using namespace std;
using namespace ThunderEgg;
using namespace ThunderEgg::VarPoisson;

class heat : public PatchOperator<2>
{
public:
    heat(std::shared_ptr<const Domain<2>>      domain,
              std::shared_ptr<const GhostFiller<2>> ghost_filler, double dt);

    void applySinglePatch(std::shared_ptr<const PatchInfo<2>> pinfo, 
                          const std::vector<LocalData<2>>& us,
                          std::vector<LocalData<2>>& fs) const override;
    void addGhostToRHS(std::shared_ptr<const PatchInfo<2>> pinfo, 
                       const std::vector<LocalData<2>>& us,
                       std::vector<LocalData<2>>& fs) const override;

    double dt;  

};


heat::heat(std::shared_ptr<const Domain<2>>      domain,
                     std::shared_ptr<const GhostFiller<2>> ghost_filler, double dt) 
                     : PatchOperator<2>(domain,ghost_filler)
{
    this->dt = dt;
}


void heat::applySinglePatch(std::shared_ptr<const PatchInfo<2>> pinfo, 
                                 const std::vector<LocalData<2>>& us,
                                 std::vector<LocalData<2>>& fs) const 
{
    //const cast since u ghost values have to be modified
    //ThunderEgg doesn't care if ghost values are modified, just don't modify the interior values.

    //fc2d_multigrid_options_t *mg_opt = fc2d_multigrid_get_options(glob);

    int mfields = us.size();
    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];

#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    
    double dx = pinfo->spacings[0];
    double dy = pinfo->spacings[1];


    for(int m = 0; m < mfields; m++)
    {
        LocalData<2>& u = const_cast<LocalData<2>&>(us[m]);
        LocalData<2>& f = fs[m];

        //if physical boundary
        if (!pinfo->hasNbr(Side<2>::west())){
            auto ghosts = u.getGhostSliceOnSide(Side<2>::west(),1);
            for(int j = 0; j < my; j++){
                ghosts[{j}] = u[{0,j}];
            }
        }
        if (!pinfo->hasNbr(Side<2>::east())){
            auto ghosts = u.getGhostSliceOnSide(Side<2>::east(),1);
            for(int j = 0; j < my; j++){
                ghosts[{j}] = u[{mx-1,j}];
            }
        }

        if (!pinfo->hasNbr(Side<2>::south())){
            auto ghosts = u.getGhostSliceOnSide(Side<2>::south(),1);
            for(int i = 0; i < mx; i++){
                ghosts[{i}] = u[{i,0}];
            }
        }
        if (!pinfo->hasNbr(Side<2>::north())){
            auto ghosts = u.getGhostSliceOnSide(Side<2>::north(),1);
            for(int i = 0; i < mx; i++){
                ghosts[{i}] = u[{i,my-1}];
            }
        }

        double dx2 = dx*dx;
        double dy2 = dy*dy;
        double dt = this->dt;
        double lambda = -1/dt;

#if 1
        /* Five-point Laplacian */
        for(int j = 0; j < my; j++)
            for(int i = 0; i < mx; i++)
            {
                double uij = u[{i,j}];
                double lap = (u[{i+1,j}] - 2*uij + u[{i-1,j}])/dx2 + 
                             (u[{i,j+1}] - 2*uij + u[{i,j-1}])/dy2;
                /* Forward Euler method */

                //f[{i,j}] = uij - dt*lap;
                f[{i,j}] = lap + lambda*uij;
            }
    
#else

        /* Five-point Laplacian : Slightly slower than then above.*/
        for(int j = 0; j < my; j++)
            for(int i = 0; i < mx; i++)
            {
                double uij = u[{i,j}];
                double flux[4];
                flux[0] = (uij - u[{i-1,j}]);
                flux[1] = (u[{i+1,j}] - uij);
                flux[2] = (uij - u[{i,j-1}]);
                flux[3] = (u[{i,j+1}] - uij);;
                f[{i,j}] = (flux[1]-flux[0])/dx2 + (flux[3] - flux[2])/dy2;
            }
#endif
    }
    
}


void heat::addGhostToRHS(std::shared_ptr<const PatchInfo<2>> pinfo, 
                              const std::vector<LocalData<2>>& us, 
                              std::vector<LocalData<2>>& fs) const 
{
#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    

    int mfields = us.size();
    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];

    double dx = pinfo->spacings[0];
    double dx2 = dx*dx;

    double dy = pinfo->spacings[1];
    double dy2 = dy*dy;

    //double dt = this->dt;

    for(int m = 0; m < mfields; m++)
    {
        const LocalData<2>& u = us[m];
        LocalData<2>& f = fs[m];
        for(int j = 0; j < my; j++)
        {
            /* bool hasNbr(Side<D> s) */
            if (pinfo->hasNbr(Side<2>::west()))
            {
                f[{0,j}] += -(u[{-1,j}]+u[{0,j}])/dx2;
                //f[{0,j}] += dt*(u[{-1,j}]+u[{0,j}])/dx2;
            }

            if (pinfo->hasNbr(Side<2>::east()))
            {                
                f[{mx-1,j}] += -(u[{mx-1,j}]+u[{mx,j}])/dx2;
                //f[{mx-1,j}] += dt*(u[{mx-1,j}]+u[{mx,j}])/dx2;
            }
        }

        for(int i = 0; i < mx; i++)
        {
            if (pinfo->hasNbr(Side<2>::south()))
            {
                f[{i,0}] += -(u[{i,-1}]+u[{i,0}])/dy2;
                //f[{i,0}] += dt*(u[{i,-1}]+u[{i,0}])/dy2;                
            }                

            if (pinfo->hasNbr(Side<2>::north()))
            {
                f[{i,my-1}] += -(u[{i,my-1}]+u[{i,my}])/dy2;
                //f[{i,my-1}] += dt*(u[{i,my-1}]+u[{i,my}])/dy2;
            }
        }
    }
}
 

void fc2d_multigrid_heat_solve(fclaw2d_global_t *glob) 
{
    // get needed options
    fclaw2d_clawpatch_options_t *clawpatch_opt =
    fclaw2d_clawpatch_get_options(glob);
    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    fc2d_multigrid_options_t *mg_opt = fc2d_multigrid_get_options(glob);
  
#if 0  
    fc2d_multigrid_vtable_t *mg_vt = fc2d_multigrid_vt();
#endif  

    // create thunderegg vector for eqn 0
    shared_ptr<Vector<2>> f = make_shared<fc2d_multigrid_vector>(glob);

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

    // ghost filler
    auto ghost_filler = make_shared<BiLinearGhostFiller>(te_domain);

    // patch operator
    double dt = glob->curr_dt;
    auto op = make_shared<heat>(te_domain,ghost_filler,dt);

    // set the patch solver
    shared_ptr<PatchSolver<2>>  solver;
    switch (mg_opt->patch_solver)
    {
        case BICG:
            solver = make_shared<BiCGStabPatchSolver<2>>(op,
                                                         mg_opt->patch_bcgs_tol,
                                                         mg_opt->patch_bcgs_max_it);
            break;
        case FFT:
            /* This ignores the five point operator defined above and just uses the 
               ThunderEgg operator 'Poisson'. */
            solver = make_shared<Poisson::FFTWPatchSolver<2>>(op);
            break;
        default:
            fclaw_global_essentialf("multigrid_heat : No valid patch solver specified\n");
            exit(0);            
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
        auto prev_domain = curr_domain;
        curr_domain = next_domain;
        while(domain_gen.hasCoarserDomain())
        {
            next_domain = domain_gen.getCoarserDomain();

            //operator
            auto ghost_filler = make_shared<BiLinearGhostFiller>(curr_domain);
            patch_operator = make_shared<heat>(curr_domain, ghost_filler,dt);

            //smoother
            shared_ptr<GMG::Smoother<2>> smoother;
            switch (mg_opt->patch_solver)
            {
                case BICG:
                    smoother = make_shared<BiCGStabPatchSolver<2>>(patch_operator,
                                                                   mg_opt->patch_bcgs_tol,
                                                                   mg_opt->patch_bcgs_max_it);
                    break;
                case FFT:
                    smoother = make_shared<Poisson::FFTWPatchSolver<2>>(patch_operator);
                    break;
                default:
                    fclaw_global_essentialf("multigrid_heat : No valid " \
                                            "patch solver specified\n");
                    exit(0);            
            }


            //restrictor
            auto restrictor = make_shared<GMG::LinearRestrictor<2>>(curr_domain, 
                                                            next_domain, clawpatch_opt->rhs_fields);

            //interpolator
            auto interpolator = make_shared<GMG::DirectInterpolator<2>>(curr_domain, 
                                                                prev_domain, clawpatch_opt->rhs_fields);

            //vector generator
            vg = make_shared<ValVectorGenerator<2>>(curr_domain, clawpatch_opt->rhs_fields);

            builder.addIntermediateLevel(patch_operator, smoother, restrictor, 
                                         interpolator, vg);

            prev_domain = curr_domain;
            curr_domain = next_domain;
        }

        //add coarsest level

        //operator
        auto ghost_filler = make_shared<BiLinearGhostFiller>(curr_domain);
        patch_operator = make_shared<heat>(curr_domain, ghost_filler,dt);

        //smoother
        switch (mg_opt->patch_solver)
        {
            case BICG:
                smoother = make_shared<BiCGStabPatchSolver<2>>(patch_operator,
                                                               mg_opt->patch_bcgs_tol,
                                                               mg_opt->patch_bcgs_max_it);
                break;
            case FFT:
                smoother = make_shared<Poisson::FFTWPatchSolver<2>>(patch_operator);
                break;
            default:
                fclaw_global_essentialf("multigrid_heat : No valid " \
                                        "patch solver specified\n");
                exit(0);            
        }


        //interpolator
        auto interpolator = make_shared<GMG::DirectInterpolator<2>>(curr_domain, prev_domain, clawpatch_opt->rhs_fields);

        //vector generator
        vg = make_shared<ValVectorGenerator<2>>(curr_domain, clawpatch_opt->rhs_fields);

        builder.addCoarsestLevel(patch_operator, smoother, interpolator, vg);

        M = builder.getCycle();
    }

    // solve
    auto vg = make_shared<ValVectorGenerator<2>>(te_domain, clawpatch_opt->rhs_fields);
    shared_ptr<Vector<2>> u = vg->getNewVector();

    int its = BiCGStab<2>::solve(vg, A, u, f, M, mg_opt->max_it, mg_opt->tol);

    fclaw_global_productionf("Iterations: %i\n", its);    

    /* Solution is copied to right hand side */
    f->copy(u);

#if 0    
    fclaw_global_productionf("f-2norm:   %24.16f\n", f->twoNorm());
    fclaw_global_productionf("f-infnorm: %24.16f\n", f->infNorm());
    fclaw_global_productionf("u-2norm:   %24.16f\n", u->twoNorm());
    fclaw_global_productionf("u-infnorm: %24.16f\n\n", u->infNorm());

    // copy solution into rhs
    fclaw_global_productionf("Checking if copy function works:\n");
    fclaw_global_productionf("fcopy-2norm:   %24.16f\n", f->twoNorm());
    fclaw_global_productionf("fcopy-infnorm: %24.16f\n\n", f->infNorm());
#endif    
}

