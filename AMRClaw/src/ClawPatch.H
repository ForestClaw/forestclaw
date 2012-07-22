#ifndef __CLAWPATCH__
#define __CLAWPATCH__

#include <iostream>
#include <vector>


using std::ifstream;
using std::ios;
using std::vector;


#include "amr_utils.H"


// The following are the parameters needed by Clawpack.  Most of these are included
// as member data items in the class ClawPatch, defined below.  Several are not used
// in the current implemenation, but are put in for future use (e.g. aux and maux).
//
//  c    meqn is the number of equations in the system of
//  c         conservation laws.
//  c
//  c    mwaves is the number of waves that result from the
//  c           solution of each Riemann problem.  Often mwaves = meqn but
//  c           for some problems these may be different, e.g. for the Euler
//  c           equations meqn = 5 but mwaves =  since there are only 3
//  c           distinct wave speeds.
//  c
//  c    mbc is the number of "ghost cells" that must be added on to each
//  c       side of the domain to handle boundary conditions.  The cells
//  c       actually in the physical domain are labelled from 1 to mx in x,
//  c       from 1 to my in y, and from 1 to mz in z.
//  c       The arrays are dimensioned actually indexed
//  c       from 1-mbc to mx+mbc, from 1-mbc to my+mbc and from 1-mbc
//  c       to mz+mbc.
//  c       For the methods currently implemented, mbc = 2 should be used.
//  c       If the user implements another method that has a larger stencil and
//  c       hence requires more ghost cells, a larger value of mbc could be used.
//  c       q is extended from the physical domain to the ghost cells by the
//  c       user-supplied routine bc2.
//  c
//  c    method(1:7) = array of values specifying the numerical method to use
//  c                  and also indicating whether source terms, capacity
//  c                  function, auxiliary variables are present in the equation.
//
//  c
//  c         method(1) = 0 if fixed size time steps are to be taken.
//  c                       In this case, dt = dtv(1) in all steps.
//  c                   = 1 if variable time steps are to be used.
//  c
//  c         method(2) = 1 if only first order increment waves are to be used.
//  c                   = 2 if second order correction terms are to be added, with
//  c                       a flux limiter as specified by mthlim.
//  c
//  c         method(3) <  0 Gives dimensional splitting using Godunov
//  c                        splitting, i.e. formally first order
//  c                        accurate.
//  c                      0 Gives the Donor cell method. No transverse
//  c                        propagation of neither the increment wave
//  c                        nor the correction wave.
//  c                   = 10 Transverse propagation of the increment wave
//  c                        as in 2D. Note that method (2,10) is
//  c                        unconditionally unstable.
//  c                   = 11 Corner transport upwind of the increment
//  c                        wave. Note that method (2,11) also is
//  c                        unconditionally unstable.
//  c                   = 20 Both the increment wave and the correction
//  c                        wave propagate as in the 2D case. Only to
//  c                        be used with method(2) = 2.
//  c                   = 21 Corner transport upwind of the increment wave,
//  c                        and the correction wave propagates as in 2D.
//  c                        Only to be used with method(2) = 2.
//  c                   = 22 3D propagation of both the increment wave and
//  c                        the correction wave. Only to be used with
//  c                        method(2) = 2.
//  c
//  c         method(4) = 0 to suppress printing
//  c                   = 1 to print dt and Courant number every time step
//  c
//  c         method(5) = 0 if there is no source term psi.  In this case
//  c                       the subroutine src2 is never called so a dummy
//  c                       parameter can be given.
//  c                   = 1 if there is a source term.  In this case
//  c                       the subroutine src2 must be provided.
//  c
//  c         method(6) = 0 if there is no capacity function capa.
//  c                   = mcapa > 0 if there is a capacity function.  In this case
//  c                       aux(i,j,k,mcapa) is the capacity of cell (i,j,k)
//  c                       and you must also specify method(7) .ge. mcapa
//  c                       and set aux.
//  c
//  c         method(7) = 0 if there is no aux array used.
//  c                   = maux > 0  if there are maux auxiliary variables.
//  c
//  c         The recommended choice of methods for most problems is
//  c            method(1) = 1,  method(2) = 2,  method(3) = 22.
//  c
//  c    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
//  c                     in each wave family mw.  Often the same value will be used
//  c                     for each value of mw, but in some cases it may be
//  c                     desirable to use different limiters.  For example,
//  c                     for the Euler equations the superbee limiter might be
//  c                     used for the contact discontinuity (mw=2) while another
//  c                     limiter is used for the nonlinear waves.  Several limiters
//  c                     are built in and others can be added by modifying the
//  c                     subroutine philim.
//  c
//  c        mthlim(mw) = 0 for no limiter
//  c                   = 1 for minmod
//  c                   = 2 for superbee
//  c                   = 3 for van Leer
//  c                   = 4 for monotonized centered
//  c
//  c    aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 1-mbc:maxmy+mbc,maux)
//  c        Array of auxiliary variables that are used in specifying the problem.
//  c        If method(7) = 0 then there are no auxiliary variables and aux
//  c                         can be a dummy variable.
//  c        If method(7) = maux > 0 then there are maux auxiliary variables
//  c                         and aux must be dimensioned as above.
//  c
//  c        Capacity functions are one particular form of auxiliary variable.
//  c        These arise in some applications, e.g. the
//  c        determinant of the Jacobian if a mapped grid is used, or a density
//  c        or porosity function in some advection problems.
//  c        See Clawpack Note # 5 for examples.
//  c
//  c        If method(6) = 0 then there is no capacity function.
//  c        If method(6) = mcapa > 0  then there is a capacity function and
//  c            capa(i,j,k), the "capacity" of the (i,j,k) cell, is assumed to be
//  c            stored in aux(i,j,k,mcapa).
//  c            In this case we require method(7).ge.mcapa.
//  c

extern "C" {

    void set_common_levels_(const int& maxlevel,
                            const int& a_level,
                            const int& refratio);

#if CH_SPACEDIM==2
    void qinit_(const int& maxmx,
                 const int& maxmy,
                 const int& meqn,
                 const int& mbc,
                 const int& mx,
                 const int& my,
                 const Real& xlower,
                 const Real& ylower,
                 const Real& dx,
                 const Real& dy,
                 Real q[],
                 const int& maux,
                 Real aux[]);


        //       subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
        //      &                  maux,aux)

    void setaux_(const int& maxmx,
                 const int& maxmy,
                 const int& mbc,
                 const int& mx,
                 const int& my,
                 const Real& xlower,
                 const Real& ylower,
                 const Real& dx,
                 const Real& dy,
                 const int& maux,
                 Real aux[]);

    void setprob_();

    void bc2_(const int& maxmx,
              const int& maxmy,
              const int& meqn,
              const int& mbc,
              const int& mx,
              const int& my,
              const Real& xlower,
              const Real& ylower,
              const Real& dx,
              const Real& dy,
              const Real q[],
              const int& maux,
              const Real aux[],
              const Real& t,
              const Real& dt,
              const int mthbc[]);

    void clawpatch2_(const int& maxm,
                     const int& meqn,
                     const int& maux,
                     const int& mbc,
                     const int method[],
                     const int mthlim[],
                     const int& mcapa,
                     const int& mwaves,
                     const int& mx,
                     const int& my,
                     Real qold[],
                     Real auxold[],
                     const Real& dx,
                     const Real& dy,
                     const Real& dt,
                     const Real& cfl,
                     Real work[],
                     const int& mwork,
                     Real qold_coarse[],
                     Real auxold_coarse[],
                     Real qadd_x[],
                     Real qadd_y[],
                     const int auxtype_int[],
                     const Real& xlower,
                     const Real& ylower,
                     int* intersectsBoundaries,
                     const int& level,
                     const int mthbc[],
                     const Real& t,
                     const int& mxc,
                     const int& myc,
                     Real fp[],
                     Real fm[],
                     Real gp[],
                     Real gm[],
                     Real fp_chombo[],
                     Real fm_chombo[],
                     Real gp_chombo[],
                     Real gm_chombo[],
                     Real fpc_chombo[],
                     Real fmc_chombo[],
                     Real gpc_chombo[],
                     Real gmc_chombo[]);


    void estimate_error2_(const int& mx,
                         const int& my,
                         const int& mbc,
                         const int& meqn,
                         const Real& xlower,
                         const Real& ylower,
                         const Real& dx,
                         const Real& dy,
                         const Real& t,
                         const int& level,
                         const int isBoundary[],
                         const Real& a_refineThresh,
                         const Real *q,
                         Real em[]);


    void write_qfile_(const int& maxmx,
                      const int& maxmy,
                      const int& meqn,
                      const int& mbc,
                      const int& mx,
                      const int& my,
                      const Real& xlower,
                      const Real& ylower,
                      const Real& dx,
                      const Real& dy,
                      Real q[],
                      const int& iframe,
                      const int& patch_num,
                      const int& level);

#elif CH_SPACEDIM==3
    void qinit_(const int& maxmx, const int& maxmy, const int& maxmz,
               const int& meqn,const int& mbc,
               const int& mx, const int& my, const int& mz,
               const Real& xlower, const Real& ylower, const Real& zlower,
               const Real& dx, const Real& dy, const Real& dz,
               Real q[], const int& maux, Real aux[]);

    void setaux_(const int& maxmx,
                 const int& maxmy,
                 const int& maxmz,
                 const int& mbc,
                 const int& mx,
                 const int& my,
                 const int& mz,
                 const Real& xlower,
                 const Real& ylower,
                 const Real& zlower,
                 const Real& dx,
                 const Real& dy,
                 const Real& dz,
                 const int& maux,
                 Real aux[]);

    void setprob_();

    void bc3_(const int& maxmx,
              const int& maxmy,
              const int& maxmz,
              const int& meqn,
              const int& mbc,
              const int& mx,
              const int& my,
              const int& mz,
              const Real& xlower,
              const Real& ylower,
              const Real& zlower,
              const Real& dx,
              const Real& dy,
              const Real& dz,
              const Real q[],
              const int& maux,
              const Real aux[],
              const Real& t,
              const Real& dt,
              const int mthbc[]);

    void clawpatch3_(const int& maxm,
                     const int& meqn,
                     const int& maux,
                     const int& mbc,
                     int method[],
                     int mthlim[],
                     const int& mcapa,
                     const int& mwaves,
                     const int& mx,
                     const int& my,
                     const int& mz,
                     Real qold[],
                     Real aux[],
                     const Real& dx,
                     const Real& dy,
                     const Real& dz,
                     const Real& dt,
                     const Real& cfl,
                     Real work[],
                     const int& mwork,
                     Real qold_coarse[],
                     Real auxold_coarse[],
                     Real qadd_x[],
                     Real qadd_y[],
                     Real qadd_z[],
                     const int auxtype_int[],
                     const Real& xlower,
                     const Real& ylower,
                     const Real& zlower,
                     int* intersectsBoundaries,
                     const int& level,
                     const int mthbc[],
                     const Real& t,
                     const int& mxc,
                     const int& myc,
                     const int& mzc,
                     Real fp[],
                     Real fm[],
                     Real gp[],
                     Real gm[],
                     Real hm[],
                     Real hp[],
                     Real fp_chombo[],
                     Real fm_chombo[],
                     Real gp_chombo[],
                     Real gm_chombo[],
                     Real hp_chombo[],
                     Real hm_chombo[],
                     Real fpc_chombo[],
                     Real fmc_chombo[],
                     Real gpc_chombo[],
                     Real gmc_chombo[],
                     Real hpc_chombo[],
                     Real hmc_chombo[]);


    void estimate_error3_(const int& mx,
                         const int& my,
                         const int& mz,
                         const int& mbc,
                         const int& meqn,
                         const Real& xlower,
                         const Real& ylower,
                         const Real& zlower,
                         const Real& dx,
                         const Real& dy,
                         const Real& dz,
                         const Real& t,
                         const int& level,
                         const int isBoundary[],
                         const Real& a_refineThresh,
                         const Real* q,
                         Real em[]);

#endif



}

class ClawPatch
{

public :
    ClawPatch();
    ~ClawPatch();

    // Do I need this?
    // void define(ProblemDomain& a_domain, const Real& a_dx);

#if CH_SPACEDIM == 2
    void define(const Real&  a_xlower,
                const Real&  a_ylower,
                const Real&  a_xupper,
                const Real&  a_yupper,
                const global_parms *a_gparms);
#else
    void define(const Real& xlower,
                const Real& ylower,
                const Real& zlower,
                const Real& xupper,
                const Real& yupper,
                const Real& zupper,
                const global_parms *a_gparms);
#endif

    bool isDefined();

    int get_mx() const;
    int get_my() const;
#if CH_SPACEDIM == 3
    int get_mz() const;
#endif

    Real get_xlower() const;
    Real get_ylower() const;
#if CH_SPACEDIM == 3
    Real get_zlower() const;
#endif

    Real get_xupper() const;
    Real get_yupper() const;
#if CH_SPACEDIM == 3
    Real get_zupper() const;
#endif

    //    void initialize(LevelData<FArrayBox>& a_phi, LevelData<FArrayBox>& a_aux, Real dx);
    // void initialize(FArrayBox& a_phi,FArrayBox& a_aux, const Box& a_box, const global_parms
    // &gparms);
    void initialize();

    void setAuxArray(const int& a_maxlevel,
                     const int& a_refratio,
                     const int& a_level);

    void write_patch_data(const int& a_iframe,
                          const int& a_patch_num,
                          const int& a_level);

    void edge_exchange_step1(const int& iside,
                             const int& refratio,
                             ClawPatch *neighbor_cp[]);

    void corner_exchange_step1(const int& icorner,
                               const int& refratio,
                               ClawPatch *corner_cp);

    void set_physbc(const int a_intersects_bc[],
                    const int a_mthbc[],
                    const Real& t,
                    const Real& dt);



    Real ClawPatchIntegrator(FArrayBox& a_phiPatch,
                             FArrayBox& a_auxPatch,
                             FArrayBox a_fluxp[],
                             FArrayBox a_fluxmc[],
                             FArrayBox a_fluxpc[],
                             FArrayBox a_fluxm[],
                             FArrayBox& a_phiCoarse,
                             FArrayBox& a_auxCoarse,
                             FArrayBox a_qadd[],
                             const Box& a_patchBox,
                             const ProblemDomain& a_domain,
                             const Real& a_time,
                             const Real& a_dt,
                             const Real& a_dx,
                             const int& refRatio,const int& a_level,
                             const global_parms &gparms);

    void estimateError(const FArrayBox& a_phiPatch,
                       const Box& a_patchBox,
                       const ProblemDomain& a_domain,
                       const Real& a_time,
                       const Real& a_dt,
                       const Real& a_dx,
                       const int& a_level,
                       const int isBoundary[],
                       const Real& a_refineThresh,
                       FArrayBox& a_error_measure,
                       const global_parms &gparms);

protected :

    int  m_mx;
    int  m_my;
#if CH_SPACEDIM == 3
    int m_mz;
#endif

    Real  m_dx;
    Real  m_dy;
#if CH_SPACEDIM == 3
    Real m_dz;
#endif

    Real m_xlower;
    Real m_ylower;
#if CH_SPACEDIM == 3
    Real m_zlower;
#endif

    Real m_xupper;
    Real m_yupper;
#if CH_SPACEDIM == 3
    Real m_zupper;
#endif

    // These are global variables, but are needed to dimensions the patch
    // and auxilliary data correctly, and so are useful to have locally.
    int m_meqn;
    int m_maux;
    int m_mbc;

    // Here is where the data is stored.
    FArrayBox m_griddata;
    FArrayBox m_auxarray;

    bool m_isDefined;

}; // End class ClawPatch

#endif
