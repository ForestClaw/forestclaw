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
//using namespace Thunderegg::Schur;
using namespace Thunderegg::VarPoisson;

double beta_coeff(const std::array<double,2>& coord)
{
    return 1.0;
}  



class fivePoint : public PatchOperator<2>
{
public:
    fivePoint(std::shared_ptr<const Domain<2>>      domain,
              std::shared_ptr<const GhostFiller<2>> ghost_filler);

    void applySinglePatch(std::shared_ptr<const PatchInfo<2>> pinfo, const LocalData<2> u,
                          LocalData<2> f) const;
    void addGhostToRHS(std::shared_ptr<const PatchInfo<2>> pinfo, const LocalData<2> u,
                       LocalData<2> f) const;


    class Generator
    {
        private:
        /**
         * @brief generator for ghost fillers
         */
        std::function<std::shared_ptr<const GhostFiller<2>>(
        std::shared_ptr<const GMG::Level<2>> level)>
        filler_gen;
        /**
         * @brief Generated operators are stored here.
         */
        std::map<std::shared_ptr<const Domain<2>>, std::shared_ptr<const fivePoint>>
        generated_operators;

        public:
        /**
         * @brief Construct a new fivePoint 2enerator
         *
         * @param finest_op the finest star pach operator
         * @param filler_gen returns a GhostFiller for a given level
         */
        Generator(std::shared_ptr<const fivePoint> finest_op,
                  std::function<
                  std::shared_ptr<const GhostFiller<2>>(std::shared_ptr<const GMG::Level<2>> level)>
                  filler_gen)
        {
            generated_operators[finest_op->domain] = finest_op;
            this->filler_gen                       = filler_gen;
        }
        /**
         * @brief Return a fivePoint 2or a given level
         *
         * @param level the level in GMG
         * @return std::shared_ptr<const fivePoint> the operator
         */
        std::shared_ptr<const fivePoint>
        operator()(std::shared_ptr<const GMG::Level<2>> level)
        {
            auto &coarser_op = generated_operators[level->getDomain()];
            if (coarser_op != nullptr) {
                return coarser_op;
            }

            std::shared_ptr<const Domain<2>> finer_domain = level->getFiner()->getDomain();
            auto                             finer_op     = generated_operators[finer_domain];
            //auto new_coeffs = PetscVector<2>::GetNewVector(level->getDomain());
            //level->getFiner()->getRestrictor().restrict(new_coeffs, finer_op->coeffs);
            coarser_op.reset(
            new fivePoint(level->getDomain(), filler_gen(level)));
            return coarser_op;
        }
    };
};


#if 0
/* Store beta as member in fivePoint */
void fivePoint::fivePoint>(beta_vec,te_domain,ghost_filler);
#endif


fivePoint::fivePoint(std::shared_ptr<const Domain<2>>      domain,
                     std::shared_ptr<const GhostFiller<2>> ghost_filler)
{
    if (domain->getNumGhostCells() < 1) {
        throw 88;
    }
    this->domain       = domain;
    this->ghost_filler = ghost_filler;
}


void fivePoint::applySinglePatch(std::shared_ptr<const PatchInfo<2>> pinfo, 
                                 const LocalData<2> u,
                                 LocalData<2> f) const 
{

    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];
#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
    double dy = pinfo->spacings[1];
#endif    
    double dx = pinfo->spacings[0];


    for(int i = 0; i < mx; i++)
    {
        for(int j = 0; j < my; j++)
        {
            f[{i,j}] = (u[{i+1,j}] + u[{i-1,j}] + u[{i,j+1}] + u[{i,j-1}] - 4*u[{i,j}])/(dx*dx);
        }
    }
}

#if 0
void fivePoint::apply(std::shared_ptr<const Vector<2>> u, std::shared_ptr<Vector<2>> f) const 
{
    //ghost_filler->fillGhost(u);
    for (auto pinfo : domain->getPatchInfoVector()) {
        applySinglePatch(pinfo, u->getLocalData(pinfo->local_index),
                         f->getLocalData(pinfo->local_index));
    }
}
#endif


void fivePoint::addGhostToRHS(std::shared_ptr<const PatchInfo<2>> pinfo, const LocalData<2> u,
                             LocalData<2> f) const 
{
    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];
#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    
    double dx = pinfo->spacings[0];
    double dx2 = dx*dx;

    double dy = pinfo->spacings[1];
    double dy2 = dy*dy;

    for(int j = 0; j < my; j++)
    {
        /* bool hasNbr(Side<D> s) */
        if (pinfo->hasNbr(Side<2>::east))
            f[{0,j}] += -u[{-1,j}]/dx2;

        if (pinfo->hasNbr(Side<2>::west))
            f[{mx-1,j}] += -u[{mx,j}]/dx2;
    }

    for(int i = 0; i < mx; i++)
    {
        if (pinfo->hasNbr(Side<2>::south))
            f[{i,0}] += -u[{i,-1}]/dy2;

        if (pinfo->hasNbr(Side<2>::north))
            f[{i,my-1}] += -u[{i,my}]/dy2;
    }
}
 


void fc2d_multigrid_solve(fclaw2d_global_t *glob) 
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
    shared_ptr<Vector<2>> f(new fc2d_multigrid_vector(glob, 0));

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
        return mg_opt->boundary_conditions[s.toInt()] == 2;
    };

    // generates levels of patches for GMG
    shared_ptr<P4estDomGen> domain_gen(new P4estDomGen(wrap->p4est, ns, 1,inf, bmf));

    // get finest level
    shared_ptr<Domain<2>> te_domain = domain_gen->getFinestDomain();

    // define operators for problems

#if 0
    // get beta function
    auto beta_func = [&](const std::array<double,2>& coord){
        double beta;
        double grad[2];
        mg_vt->fort_beta(&coord[0],&coord[1],&beta,grad);
        return beta;
    };
#endif  

    // create vector for beta
    auto beta_vec = ValVector<2>::GetNewVector(te_domain);
    DomainTools<2>::setValuesWithGhost(te_domain, beta_vec, beta_coeff);

    // ghost filler
    auto ghost_filler = make_shared<BiLinearGhostFiller>(te_domain);
    // patch operator
    auto op = make_shared<fivePoint>(te_domain,ghost_filler);

    // set the patch solver
    auto  solver = make_shared<BiCGStabPatchSolver<2>>(te_domain, ghost_filler, op,
                                                       mg_opt->patch_bcgs_tol,
                                                       mg_opt->patch_bcgs_max_it);

    // create matrix
    shared_ptr<Operator<2>> A = op;

    // create gmg preconditioner
    shared_ptr<Operator<2>> M;
#if 1  
    if(mg_opt->mg_prec)
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

        //generators for levels
        BiLinearGhostFiller::Generator filler_gen(ghost_filler);
        fivePoint::Generator   op_gen(op, filler_gen);
        BiCGStabPatchSolver<2>::Generator smooth_gen(solver, filler_gen, op_gen);
        GMG::LinearRestrictor<2>::Generator        restrictor_gen;
		GMG::DrctIntp<2>::Generator       interpolator_gen;
        M = GMG::CycleFactory<2>::getCycle(copts, domain_gen, restrictor_gen, interpolator_gen,
                                           smooth_gen, op_gen);
    }
#endif

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

