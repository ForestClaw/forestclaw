#ifndef __CLAWPATCH_H
#define __CLAWPATCH_H

#include <iostream>
#include <vector>

using namespace std;

#include "amr_options.h"
#include "clawpack_fort.H"
#include "fclaw2d_defs.H"

class ClawPatch;

typedef void (*fclaw2d_solver_patch_data_new_t)(void **solverdata);
typedef void (*fclaw2d_solver_patch_data_delete_t)(void **solverdata);

void set_clawpatch(fclaw2d_domain_t* domain, fclaw2d_patch_t *this_patch,
                   int blockno, int patchno);

void delete_clawpatch(fclaw2d_patch_t* this_patch);

ClawPatch* get_clawpatch(fclaw2d_patch_t *patch);


size_t pack_size(fclaw2d_domain_t* domain);
void pack_clawpatch(fclaw2d_patch_t* this_patch,double* data);
void unpack_clawpatch(fclaw2d_domain_t* domain, fclaw2d_patch_t* this_patch,
                      int this_block_idx, int this_patch_idx, double *qdata);

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
                const int& a_level,
                const amr_options_t *a_gparms);

    void copyFrom(ClawPatch *a_cp);

    // ----------------------------------------------------------------
    // Time stepping routines
    // ----------------------------------------------------------------

    void save_step();
    void save_current_step();
    void restore_step();

    void time_interpolate(const double& alpha);

    // ----------------------------------------------------------------
    // Internal boundary conditions
    // ----------------------------------------------------------------
    void exchange_face_ghost(const int& a_idir,ClawPatch *a_neighbor_cp);

    void mb_exchange_face_ghost(const int& a_iface,ClawPatch *a_neighbor_cp);

    void mb_exchange_corner_ghost(const int& a_icorner, fclaw_bool a_intersects_block[],
                                  ClawPatch *a_neighbor_cp,
                                  const fclaw_bool& a_is_block_corner);

    // Average finer grid values onto coarser grid ghost cells
    void average_face_ghost(const int& a_idir,
                            const int& a_iside,
                            const int& a_num_neighbors,
                            const int& a_refratio,
                            ClawPatch **neighbor_cp,
                            fclaw_bool a_time_interp, fclaw_bool a_block_boundary);

    // Interpoalte coarser values onto fine grid ghost
    void interpolate_face_ghost(const int& a_idir,
                                const int& a_iside,
                                const int& a_num_neighbors,
                                const int& a_refratio,
                                ClawPatch **a_neighbor_cp,
                                fclaw_bool a_time_interp,
                                fclaw_bool a_block_boundary);

    // Exchange corner ghost cells with interior neighbor (interior to the domain)
    void exchange_corner_ghost(const int& a_corner, ClawPatch *cp_neighbor);
    void average_corner_ghost(const int& a_corner, const int& a_refratio,
                              ClawPatch *cp_fine,
                              fclaw_bool a_time_interp);

    void mb_average_corner_ghost(const int& a_corner, const int& a_refratio,
                                 ClawPatch *cp_fine,
                                 fclaw_bool a_time_interp, fclaw_bool is_corner_block,
                                 fclaw_bool intersects_block[]);

    void interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                  ClawPatch *cp_fine,
                                  fclaw_bool a_time_interp);

    void mb_interpolate_corner_ghost(const int& a_corner, const int& a_refratio,
                                     ClawPatch *cp_fine,
                                     fclaw_bool a_time_interp, fclaw_bool is_block_corner,
                                     fclaw_bool intersects_block[]);

    // ----------------------------------------------------------------------------------
    // Physical boundary conditions
    // ----------------------------------------------------------------------------------

    // Boundary condition routines that don't require neighbors
    void set_phys_face_ghost(const fclaw_bool a_intersects_bc[],
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
    /*
    fclaw_bool tag_for_refinement(fclaw_bool init_flag);
    fclaw_bool tag_for_coarsening(ClawPatch *cp_siblings[],
                            const int& refratio,
                            const int& num_siblings,
                            const int& p4est_refineFactor);
    */

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

    void setup_manifold(const int& a_level, const amr_options_t *gparms);

    // ----------------------------------------------------------------
    // Output and diagnostics
    // ----------------------------------------------------------------
    void write_patch_data(const int& a_iframe,
                          const int& a_patch_num,
                          const int& a_level);
    double compute_sum();

    void dump(int mq);
    void dump_last();
    void dump_time_interp();

    // ----------------------------------------------------------------
    // Access functions
    // ----------------------------------------------------------------
    double dx();
    double dy();

#if FCLAW_SPACEDIM == 3
    double dz();
#endif

    double *xp();
    double *yp();
    double *zp();
    double *xd();
    double *yd();
    double *zd();

    double *area();

    double *xface_normals();
    double *yface_normals();
    double *xface_tangents();
    double *yface_tangents();
    double *surf_normals();
    double *curvature();
    double *edge_lengths();

    double xlower();
    double ylower();

    double* q();


    /* ---------------------------------------------------------------
       This isn't very object oriented, and I will think about a better
       way to do this without deriving classes
    */

    /* ----------------------------------------------------------------
       Wave propagation routines
       Only defined if we are really using the wave prop algorithm.
       The new/delete functions are static because we want to define them
       ---------------------------------------------------------------- */
    void* waveprop_patch_data();
    void set_waveprop_patch_data(void* solver_data);
    static fclaw2d_solver_patch_data_new_t f_waveprop_patch_data_new;
    static fclaw2d_solver_patch_data_delete_t f_waveprop_patch_data_delete;

    // ----------------------------------------------------------------
    // Wave propagation used for ManyClaw
    // Only defined if we are really using ManyClaw
    // The new/delete functions are static because we want to define them
    // ----------------------------------------------------------------
    void* manyclaw_patch_data();
    void set_manyclaw_patch_data(void* solver_data);
    static fclaw2d_solver_patch_data_new_t f_manyclaw_patch_data_new;
    static fclaw2d_solver_patch_data_delete_t f_manyclaw_patch_data_delete;

    // ----------------------------------------------------------------
    // User data
    // ----------------------------------------------------------------
    void* user_patch_data();
    void set_user_patch_data(void* user_data);
    static fclaw2d_solver_patch_data_new_t f_user_patch_data_new;
    static fclaw2d_solver_patch_data_delete_t f_user_patch_data_delete;


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

    /* These are global variables, but are needed to dimensions the patch, so we keep
       them here.*/
    int m_meqn;
    int m_mbc;

    fclaw_bool m_manifold;
    // fclaw_bool m_mapped;
    int m_blockno; /* This determines which mapping to use for this patch */

    /* Current time step */
    FArrayBox m_griddata;

    /* Time interpolated grid data */
    FArrayBox m_griddata_time_interp;

    /* Last time step, needed to get time interpolated data */
    FArrayBox m_griddata_last;

    /* Save last time step in case we need to retake a time step */
    FArrayBox m_griddata_save;

    FArrayBox m_xp;
    FArrayBox m_yp;
    FArrayBox m_zp;

    FArrayBox m_xd;
    FArrayBox m_yd;
    FArrayBox m_zd;

    FArrayBox m_area;
    FArrayBox m_xface_normals;
    FArrayBox m_yface_normals;
    FArrayBox m_xface_tangents;
    FArrayBox m_yface_tangents;
    FArrayBox m_surf_normals;
    FArrayBox m_edge_lengths;

    FArrayBox m_curvature;  // ???

    /* Data associated with the wave propagation algorithm (ClawPack) */
    void *m_waveprop_patch_data;
    /* Data associated with the wave propagation algorithm (ClawPack) */
    void *m_manyclaw_patch_data;

    // Add user data (for coefficients, for example)
    void * m_user_patch_data;
}; // End class ClawPatch

#endif
