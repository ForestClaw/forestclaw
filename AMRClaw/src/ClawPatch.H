#ifndef __CLAWPATCH_H
#define __CLAWPATCH_H

#include <iostream>
#include <vector>

using namespace std;

#include "amr_options.h"
#include "clawpack_fort.H"
#include "fclaw_defs.H"

class ClawPatch
{
public :

    ClawPatch();
    ~ClawPatch();

    void define(const double&  a_xlower,
                const double&  a_ylower,
                const double&  a_xupper,
                const double&  a_yupper,
                const int& a_blockno,
                const amr_options_t *a_gparms);

    void copyFrom(ClawPatch *a_cp);

    double* current_data_ptr();
    bool isDefined();  // Do I even use this?

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
    double step_dummy(const double& a_time,
                      const double& a_dt,
                      const int& a_level,
                      const amr_options_t &gparms);

    // Actual step routine with no conservative fix-up
    double step_waveprop(const double& a_time,
                         const double& a_dt,
                         const int& a_level,
                         const amr_options_t *gparms);

    void step_diffusion(const double& a_time,
                        const double& a_dt,
                        const int& a_level,
                        const amr_options_t *gparms);


    // Full time step routine with conservative fix-up (not yet implemented)
    double ClawPatchIntegrator(const double& a_time,
                             const double& a_dt,
                             const int& refRatio,const int& a_level,
                             const amr_options_t &gparms);

    void save_step();
    void restore_step();

    /*
    void time_interpolate(const int& a_fine_step, const int& a_coarse_step,
                          const int& a_step_inc);
    */

    void time_interpolate(const double& alpha);

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
                             const double& t,
                             const double& dt);

    // Corners that lie on a physical face - data has to be exchanged
    // in a particular way
    void set_phys_corner_ghost(const int& a_corner, const int a_mthbc[],
                               const double& t, const double& dt);

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
    double compute_sum();

    void dump();
    void dump_last();
    void dump_time_interp();
    void dump_auxarray();

    double dx();
    double dy();



#if FCLAW_SPACEDIM == 3
    double dz();
#endif


protected :

    int  m_mx;
    int  m_my;
#if FCLAW_SPACEDIM == 3
    int m_mz;
#endif

    double  m_dx;
    double  m_dy;
#if FCLAW_SPACEDIM == 3
    double m_dz;
#endif

    double m_xlower;
    double m_ylower;
#if FCLAW_SPACEDIM == 3
    double m_zlower;
#endif

    double m_xupper;
    double m_yupper;
#if FCLAW_SPACEDIM == 3
    double m_zupper;
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
