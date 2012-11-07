#ifndef __CLAWPATCH__
#define __CLAWPATCH__

#include <iostream>
#include <vector>

using namespace std;

#include "amr_utils.H"
#include "fclaw_defs.H"

// Think about making ClawPatch members all static

extern "C" {

    void set_common_levels_(const int& maxlevel,const int& a_level,
                            const int& refratio);

    // ----------------------------------------------------------------------------------
    // Initialization routines
    // ----------------------------------------------------------------------------------
    void setprob_();

    void set_maptype_();
    void set_block_(const int * a_blockno);

    void qinit_(const int& maxmx, const int& maxmy, const int& meqn,
                const int& mbc, const int& mx, const int& my,
                const Real& xlower, const Real& ylower, const Real& dx,
                const Real& dy, Real q[], const int& maux, Real aux[]);

    void qinit_mapped_(const int& mx, const int& my, const int& meqn,
                       const int& mbc, const Real& xlower, const Real& ylower,
                       const Real& dx, const Real& dy, Real xp[], Real yp[], Real zp[],
                       Real q[], const int& maux, Real aux[], const int& blockno);

    void setaux_(const int& maxmx, const int& maxmy, const int& mbc,
                 const int& mx, const int& my, const Real& xlower,
                 const Real& ylower, const Real& dx, const Real& dy,
                 const int& maux, Real aux[]);

    void setaux_mapped_(const int& mx, const int& my, const int& mbc,
                        const Real& dx, const Real& dy,
                        Real xp[], Real yp[], Real zp[],
                        Real xd[], Real yd[], Real zd[],
                        Real area[], const int& maux, Real aux[]);

    void b4step2_mapped_(const int& mx, const int& my, const int& mbc,
                         const int& meqn, Real q[], const Real& dx, const Real& dy,
                         Real xp[], Real yp[], Real zp[],
                         Real xd[], Real yd[], Real zd[], const Real& t,
                         const Real& dt, const int& maux, Real aux[]);

    void average_face_ghost_mapped_(const int& mx, const int& my, const int& mbc,
                                    const int& meqn, Real qfine[],Real qcoarse[],
                                    const Real auxcoarse[], const Real auxfine[],
                                    const int& maux, const int& idir, const int& iside,
                                    const int& num_neighbors, const int& refratio,
                                    const int& igrid);

    void mb_average_face_ghost_(const int& mx, const int& my, const int& mbc,
                                const int& meqn, Real qfine[],Real qcoarse[],
                                const Real auxcoarse[], const Real auxfine[], const int& maux,
                                const int& idir, const int& iside, const int& num_neighbors,
                                const int& refratio, const int& igrid);

    void fixcapaq2_(const int& mx, const int& my, const int& mbc, const int& meqn,
                    Real qcoarse[], Real qfine[], Real auxcoarse[], Real auxfine[],
                    const int& maux, const int& p4est_refineFactor,
                    const int& refratio, const int& igrid);

    // ----------------------------------------------------------------------------------
    // Internal boundary conditions
    // ----------------------------------------------------------------------------------
    void exchange_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                              Real qthis[],Real qneighbor[], const int& a_idir);

    void average_face_ghost_(const int& mx, const int& my, const int& mbc,
                             const int& meqn, Real qfine[],Real qcoarse[],
                             const int& idir, const int& iside, const int& num_neighbors,
                             const int& refratio, const int& igrid);

    void interpolate_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                                 Real qthis[],Real qcoarse[], const int& idir, const int& iside,
                                 const int& num_neighbors,const int& refratio, const int& igrid);


    void exchange_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                const int& meqn, Real this_q[],Real neighbor_q[],
                                const int& a_corner);

    void average_corner_ghost_(const int& mx, const int& my, const int& mbc,
                               const int& meqn, const int& a_refratio, Real this_q[],
                               Real neighbor_q[],const int& a_corner);

    void interpolate_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                   const int& meqn, const int& a_refratio, Real this_q[],
                                   Real neighbor_q[], const int& a_corner);

    void mb_exchange_face_ghost_(const int& mx, const int& my, const int& mbc, const int& meqn,
                                 Real qthis[], Real qneighbor[], const int& iface,
                                 const int& iblock);

    void mb_exchange_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                   const int& meqn,
                                   Real qthis[], Real qneighbor[], const int& icorner,
                                   int bdry[], const int& iblock);

    void mb_exchange_block_corner_ghost_(const int& mx, const int& my,
                                         const int& mbc, const int& meqn,
                                         Real qthis[], Real qneighbor[], const int& icorner,
                                         const int& iblock);

    void mb_interpolate_face_ghost_(const int& mx, const int& my, const int& mbc,
                                    const int& meqn,
                                    Real qthis[],Real qcoarse[], const int& idir,
                                    const int& iside,
                                    const int& num_neighbors,const int& refratio,
                                    const int& igrid);

    void mb_average_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                  const int& meqn,
                                  const int& refratio, Real qcoarse[],  Real qfine[],
                                  Real auxcoarse[], Real auxfine[], const int& maux,
                                  const int& icorner, int intersects_block[]);

    // Averaging at block boundaries between coarse and fine grids.
    void  mb_average_block_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                         const int& meqn,const int& refratio, Real qcoarse[],
                                         Real qfine[],Real auxcoarse[], Real auxfine[],
                                         const int& maux, const int& a_coarse_corner,
                                         const int& blockno);


    void mb_interpolate_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                      const int& meqn, const int& refratio,
                                      Real qcoarse[],  Real qfine[], const int& icorner,
                                      int intersects_block[]);

    // Averaging at block boundaries between coarse and fine grids.
    void  mb_interpolate_block_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                             const int& meqn,const int& refratio,
                                             Real qcoarse[],
                                             Real qfine[],const int& a_coarse_corner,
                                             const int& blockno);


    // ----------------------------------------------------------------------------------
    // Physical boundary conditions
    // ----------------------------------------------------------------------------------
    void bc2_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
              const int& mx, const int& my, const Real& xlower, const Real& ylower,
              const Real& dx, const Real& dy, const Real q[], const int& maux,
              const Real aux[], const Real& t, const Real& dt, const int mthbc[]);

    void set_phys_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                const int& meqn, Real q[],const int& icorner,
                                const Real &t, const Real& dt, const int mthbc[]);

    void exchange_phys_corner_ghost_(const int& mx, const int& my, const int& mbc,
                                     const int& meqn, Real qthis[],Real qneighbor[],
                                     const int& icorner, const int& iside);

    // ----------------------------------------------------------------------------------
    // Time stepping
    // ----------------------------------------------------------------------------------
    void clawpatch2_(const int& maxm, const int& meqn, const int& maux, const int& mbc,
                     const int method[], const int mthlim[], const int& mcapa,
                     const int& mwaves, const int& mx, const int& my, Real qold[],
                     Real auxold[], const Real& dx, const Real& dy, const Real& dt,
                     const Real& cfl, Real work[], const int& mwork, const Real& xlower,
                     const Real& ylower, const int& level, const Real& t, Real fp[],
                     Real fm[], Real gp[], Real gm[]);

    // ----------------------------------------------------------------------------------
    // Tagging for refinement/coarsening
    // ----------------------------------------------------------------------------------
    void tag_for_refinement_(const int& mx,const int& my,const int& mbc,
                             const int& meqn, const Real& xlower,
                             const Real& ylower, const Real& dx, const Real& dy,
                             Real q[], const int& init_flag, int& tag_patch);

    void tag_for_coarsening_(const int& mx,const int& my,const int& mbc,
                             const int& meqn,const Real& xlower,
                             const Real& ylower, const Real& dx, const Real& dy,
                             Real qcoarsened[], int& tag_patch);

    void interpolate_to_fine_patch_(const int& mx,const int& my,const int& mbc,
                                    const int& meqn, Real qcoarse[], Real qfine[],
                                    const int& p4est_refineFactor,
                                    const int& refratio, const int& igrid);

    void average_to_coarse_patch_(const int& mx,const int& my,const int& mbc,
                                  const int& meqn, Real qcoarse[],Real qfine[],
                                  const int& p4est_refineFactor,
                                  const int& refratio, const int& igrid);

    void average_to_coarse_mapped_(const int& mx,const int& my, const int& mbc,
                                   const int& meqn, Real qcoarse[], Real qfine[],
                                   const Real auxcoarse[], const Real auxfine[],
                                   const int& maux, const int& p4est_refineFactor,
                                   const int& refratio, const int& igrid);

// ----------------------------------------------------------------------------------
// Mapped grids
// ----------------------------------------------------------------------------------
    void setup_mesh_(const int& mx, const int& my, const int& mbc,
                     const Real& xlower, const Real& ylower,
                     const Real& dx, const Real& dy,
                     Real xp[], Real yp[], Real zp[], Real xd[], Real yd[], Real zd[]);

    void compute_area_(const int& mx, const int& my, const int& mbc,
                       const Real& m_dx, const Real& m_dy,const Real& m_xlower,
                       const Real& m_ylower, Real area[],
                       const int& level, const int& maxlevel, const int& refratio);

    // ----------------------------------------------------------------------------------
    // Output and diagnostics
    // ----------------------------------------------------------------------------------
    void write_qfile_(const int& maxmx, const int& maxmy, const int& meqn, const int& mbc,
                      const int& mx, const int& my, const Real& xlower, const Real& ylower,
                      const Real& dx, const Real& dy, Real q[], const int& iframe,
                      const int& patch_idx,
                      const int& level,const int& blockno);

    void compute_sum_(const int& mx, const int& my, const int& mbc, const int& meqn,
                      const Real& dx,
                      const Real& dy, Real q[], Real& sum);


}

class ClawPatch
{
public :

    ClawPatch();
    ~ClawPatch();

    void define(const Real&  a_xlower,
                const Real&  a_ylower,
                const Real&  a_xupper,
                const Real&  a_yupper,
                const int& a_blockno,
                const amr_options_t *a_gparms);

    void copyFrom(ClawPatch *a_cp);

    bool isDefined();

    // ----------------------------------------------------------------------------------
    // Initialization routines
    // ----------------------------------------------------------------------------------
    void initialize();

    void setup_patch(const int& a_level, const int& a_maxlevel,
                     const int& a_refratio);

    void setAuxArray();

    // ----------------------------------------------------------------
    // Time stepping routines
    // ----------------------------------------------------------------

    // Dummy step routine useful for trying out subcycling without
    // anything too terrible happening
    Real step(const Real& a_time,
              const Real& a_dt,
              const int& a_level,
              const amr_options_t &gparms);

    // Actual step routine with no conservative fix-up
    Real step_noqad(const Real& a_time,
                    const Real& a_dt,
                    const int& a_level,
                    const amr_options_t &gparms);


    // Full time step routine with conservative fix-up (not yet implemented)
    Real ClawPatchIntegrator(const Real& a_time,
                             const Real& a_dt,
                             const int& refRatio,const int& a_level,
                             const amr_options_t &gparms);

    void save_step();
    void restore_step();

    void time_interpolate(const int& a_fine_step, const int& a_coarse_step,
                          const int& a_step_inc);

    // ----------------------------------------------------------------
    // Internal boundary conditions
    // ----------------------------------------------------------------
    void exchange_face_ghost(const int& a_idir,ClawPatch *a_neighbor_cp);

    void mb_exchange_face_ghost(const int& a_iface,ClawPatch *a_neighbor_cp);

    void mb_exchange_corner_ghost(const int& a_icorner, bool a_intersects_block[],
                                  ClawPatch *a_neighbor_cp,
                                  const bool& a_is_block_corner);

    // Average finer grid values onto coarser grid ghost cells
    void average_face_ghost(const int& a_idir,
                            const int& a_iside,
                            const int& a_num_neighbors,
                            const int& a_refratio,
                            ClawPatch **neighbor_cp,
                            bool a_time_interp, bool a_block_boundary);

    // Interpoalte coarser values onto fine grid ghost
    void interpolate_face_ghost(const int& a_idir,
                                const int& a_iside,
                                const int& a_num_neighbors,
                                const int& a_refratio,
                                ClawPatch **a_neighbor_cp,
                                bool a_time_interp,
                                bool a_block_boundary);

    // Exchange corner ghost cells with interior neighbor (interior to the domain)
    void exchange_corner_ghost(const int& a_corner, ClawPatch *cp_neighbor);
    void average_corner_ghost(const int& a_corner, const int& a_refratio,
                              ClawPatch *cp_fine,
                              bool a_time_interp);

    void mb_average_corner_ghost(const int& a_corner, const int& a_refratio,
                                 ClawPatch *cp_fine,
                                 bool a_time_interp, bool is_corner_block,
                                 bool intersects_block[]);

    void interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                  ClawPatch *cp_fine,
                                  bool a_time_interp);

    void mb_interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                     ClawPatch *cp_fine,
                                     bool a_time_interp, bool is_block_corner,
                                     bool intersects_block[]);

    // ----------------------------------------------------------------------------------
    // Physical boundary conditions
    // ----------------------------------------------------------------------------------

    // Boundary condition routines that don't require neighbors
    void set_phys_face_ghost(const bool a_intersects_bc[],
                             const int a_mthbc[],
                             const Real& t,
                             const Real& dt);

    // Corners that lie on a physical face - data has to be exchanged
    // in a particular way
    void set_phys_corner_ghost(const int& a_corner, const int a_mthbc[],
                               const Real& t, const Real& dt);

    // Exchange corner ghost with boundary neighbor
    void exchange_phys_face_corner_ghost(const int& a_corner,
                                         const int& a_iface, ClawPatch* cp);

    // ----------------------------------------------------------------
    // Tagging, refining and coarsening
    // ----------------------------------------------------------------
    bool tag_for_refinement(bool init_flag);
    bool tag_for_coarsening(ClawPatch *cp_siblings[],
                            const int& refratio,
                            const int& num_siblings,
                            const int& p4est_refineFactor);

    void interpolate_to_fine_patch(ClawPatch* a_fine,
                                   const int& a_igrid,
                                   const int& a_p4est_refineFactor,
                                   const int& a_refratio);

    void coarsen_from_fine_family(ClawPatch *a_cp_siblings[],
                                  const int& a_refratio,
                                  const int& a_num_siblings,
                                  const int& a_p4est_refineFactor);


    // ----------------------------------------------------------------
    // Mapped grids
    // ----------------------------------------------------------------

    void setup_manifold(const int& a_level, const int& a_maxlevel,
                        const int& a_refratio);

    // ----------------------------------------------------------------
    // Output and diagnostics
    // ----------------------------------------------------------------
    void write_patch_data(const int& a_iframe,
                          const int& a_patch_num,
                          const int& a_level);
    Real compute_sum();

    void dump();
    void dump_last();
    void dump_time_interp();
    void dump_auxarray();

protected :

    int  m_mx;
    int  m_my;
#if FCLAW_SPACEDIM == 3
    int m_mz;
#endif

    Real  m_dx;
    Real  m_dy;
#if FCLAW_SPACEDIM == 3
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

    bool m_manifold;
    bool m_mapped;
    int m_blockno; // This determines which mapping to use for this patch

    // Here is where the data is stored.
    FArrayBox m_griddata;
    FArrayBox m_griddata_time_interp;
    FArrayBox m_griddata_last;  // stores previous time step for interpolation in time
    FArrayBox m_griddata_save;

    FArrayBox m_auxarray;

    FArrayBox m_xp;
    FArrayBox m_yp;
    FArrayBox m_zp;

    FArrayBox m_xd;
    FArrayBox m_yd;
    FArrayBox m_zd;
    FArrayBox m_area;

    bool m_isDefined;

}; // End class ClawPatch

#endif
