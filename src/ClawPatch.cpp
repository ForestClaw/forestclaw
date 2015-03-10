#include "ClawPatch.H"
#include "fclaw2d_typedefs.h"
#include "fclaw2d_timeinterp.H"

fclaw_app_t* ClawPatch::app;

ClawPatch::ClawPatch()
{
    m_package_data_ptr = fclaw_package_data_new();
}

ClawPatch::~ClawPatch()
{
    fclaw_package_patch_data_destroy(ClawPatch::app,
                                     m_package_data_ptr);

    fclaw_package_data_destroy(m_package_data_ptr);
}


void ClawPatch::define(const double&  a_xlower,
                       const double&  a_ylower,
                       const double&  a_xupper,
                       const double&  a_yupper,
                       const int& a_blockno,
                       const int& a_level,
                       const amr_options_t* gparms)
{
    m_mx = gparms->mx;
    m_my = gparms->my;
    m_mbc = gparms->mbc;
    m_blockno = a_blockno;
    m_meqn = gparms->meqn;

    m_manifold = gparms->manifold;

    if (m_manifold)
    {
        m_xlower = a_xlower;
        m_ylower = a_ylower;
        m_xupper = a_xupper;
        m_yupper = a_yupper;
    }
    else
    {
        double ax = gparms->ax;
        double bx = gparms->bx;
        double ay = gparms->ay;
        double by = gparms->by;
        m_xlower = ax + (bx - ax)*a_xlower;
        m_xupper = ax + (bx - ax)*a_xupper;
        m_ylower = ay + (by - ay)*a_ylower;
        m_yupper = ay + (by - ay)*a_yupper;
    }

    m_dx = (m_xupper - m_xlower)/m_mx;
    m_dy = (m_yupper - m_ylower)/m_my;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-m_mbc;
    }
    ur[0] = m_mx + m_mbc;
    ur[1] = m_my + m_mbc;
    Box box(ll,ur);

    // This will destroy any existing memory n m_griddata.
    m_griddata.define(box, m_meqn);
    m_griddata_last.define(box, m_meqn);
    m_griddata_save.define(box, m_meqn);
    m_griddata_time_interpolated.define(box, m_meqn);

    if (m_manifold)
    {
        setup_manifold(a_level,gparms);
    }

    fclaw_package_patch_data_new(ClawPatch::app,m_package_data_ptr);

}

void ClawPatch::copyFrom(ClawPatch *a_cp)
{
    m_griddata = a_cp->m_griddata;
}

// This is used by level_step.
double* ClawPatch::q()
{
    return m_griddata.dataPtr();
}

FArrayBox ClawPatch::newGrid()
{
    /* Create a grid based on the size of the existing grids */
    Box b = m_griddata.box();
    int fields = m_griddata.fields();
    FArrayBox A;
    A.define(b,fields);
    return A;
}

Box ClawPatch::dataBox()
{
    return m_griddata.box();
}

Box ClawPatch::areaBox()
{
    return m_area.box();
}

Box ClawPatch::edgeBox()
{
    return m_edge_lengths.box();
}

Box ClawPatch::nodeBox()
{
    return m_xp.box();
}


/* Return a pointer to either time interpolated data or regular grid data */
double* ClawPatch::q_time_sync(fclaw_bool time_interp)
{
    if (time_interp)
        return m_griddata_time_interpolated.dataPtr();
    else
        return m_griddata.dataPtr();
}

double* ClawPatch::q_time_interp()
{
    return m_griddata_time_interpolated.dataPtr();
}

void ClawPatch::save_current_step()
{
    m_griddata_last = m_griddata; // Copy for time interpolation
}

double ClawPatch::mx()
{
    return m_mx;
}
double ClawPatch::my()
{
    return m_my;
}
double ClawPatch::mbc()
{
    return m_mbc;
}
double ClawPatch::meqn()
{
    return m_meqn;
}

double ClawPatch::dx()
{
    return m_dx;
}

double ClawPatch::dy()
{
    return m_dy;
}

double ClawPatch::xlower()
{
    return m_xlower;
}

double ClawPatch::ylower()
{
    return m_ylower;
}

double ClawPatch::xupper()
{
    return m_xupper;
}

double ClawPatch::yupper()
{
    return m_yupper;
}

double* ClawPatch::xp()
{
    return m_xp.dataPtr();
}

double* ClawPatch::yp()
{
    return m_yp.dataPtr();
}
double* ClawPatch::zp()
{
    return m_zp.dataPtr();
}

double* ClawPatch::xd()
{
    return m_xd.dataPtr();
}
double* ClawPatch::yd()
{
    return m_yd.dataPtr();
}
double* ClawPatch::zd()
{
    return m_zd.dataPtr();
}

double* ClawPatch::area()
{
    return m_area.dataPtr();
}

double* ClawPatch::xface_normals()
{
    return m_xface_normals.dataPtr();
}

double* ClawPatch::yface_normals()
{
    return m_yface_normals.dataPtr();
}

double* ClawPatch::xface_tangents()
{
    return m_xface_tangents.dataPtr();
}

double* ClawPatch::yface_tangents()
{
    return m_yface_tangents.dataPtr();
}

double* ClawPatch::surf_normals()
{
    return m_surf_normals.dataPtr();
}

double* ClawPatch::edge_lengths()
{
    return m_edge_lengths.dataPtr();
}

double* ClawPatch::curvature()
{
    return m_curvature.dataPtr();
}

int* ClawPatch::block_corner_count()
{
    return &m_block_corner_count[0];
}

/* ----------------------------------------------------
   Solver data and functions
   ---------------------------------------------------*/
// Wave propagation algorithms

void* ClawPatch::clawpack_patch_data(int id)
{
    return fclaw_package_get_data(m_package_data_ptr,id);
    // return m_clawpack_patch_data;
}

/* ----------------------------------------------------------------
   Time stepping routines
   ---------------------------------------------------------------- */


void ClawPatch::save_step()
{
    // Store a backup in case the CFL number is too large doesn't work out.
    m_griddata_save = m_griddata;
}

void ClawPatch::restore_step()
{
    m_griddata = m_griddata_save;
}

void ClawPatch::pack_griddata(double *q)
{
    m_griddata.copyToMemory(q);
}

void ClawPatch::unpack_griddata(double *q)
{
    m_griddata.copyFromMemory(q);
}

void ClawPatch::unpack_griddata_time_interpolated(double *q)
{
    m_griddata_time_interpolated.copyFromMemory(q);
}


/* ----------------------------------------------------------------
   Exchange/average/interpolate
   ---------------------------------------------------------------- */

void ClawPatch::exchange_face_ghost(const int& a_iface,
                                    ClawPatch *neighbor_cp,
                                    fclaw2d_transform_data_t* transform_data)
{
    double *qthis = m_griddata.dataPtr();
    double *qneighbor = neighbor_cp->m_griddata.dataPtr();
    exchange_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qthis,qneighbor,a_iface,
                         &transform_data);
}


void ClawPatch::exchange_corner_ghost(const int& a_corner, ClawPatch *cp_corner,
                                      fclaw2d_transform_data_t* transform_data)
{
    double *qthis = m_griddata.dataPtr();
    double *qcorner = cp_corner->m_griddata.dataPtr();

    exchange_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner, a_corner,
                           &transform_data);

}

/* ----------------------------------------------------------------
   Multi-level operations
   ---------------------------------------------------------------- */
void ClawPatch::average_face_ghost(const int& a_idir,
                                   const int& a_iface_coarse,
                                   const int& a_p4est_refineFactor,
                                   const int& a_refratio,
                                   ClawPatch *neighbor_cp,
                                   fclaw_bool a_time_interp,
                                   const int& igrid,
                                   fclaw2d_transform_data_t* transform_data)
{
    double *qcoarse = q_time_sync(a_time_interp);

    double *qfine = neighbor_cp->m_griddata.dataPtr();

    /* These will be empty for non-manifolds cases */
    double *areacoarse = m_area.dataPtr();
    double *areafine = neighbor_cp->m_area.dataPtr();

    int manifold = m_manifold ? 1 : 0;
    average_face_ghost_(m_mx,m_my,m_mbc,m_meqn,
                        qcoarse,qfine,
                        areacoarse, areafine,
                        a_idir,a_iface_coarse,
                        a_p4est_refineFactor,a_refratio,igrid,
                        manifold, &transform_data);
}

void ClawPatch::interpolate_face_ghost(const int& a_idir,
                                       const int& a_iside,
                                       const int& a_p4est_refineFactor,
                                       const int& a_refratio,
                                       ClawPatch *neighbor_cp,
                                       fclaw_bool a_time_interp,
                                       const int& igrid,
                                       fclaw2d_transform_data_t* transform_data)
{
    double *qcoarse = q_time_sync(a_time_interp);
    double *qfine = neighbor_cp->m_griddata.dataPtr();

    interpolate_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iside,
                            a_p4est_refineFactor,a_refratio,igrid,
                            &transform_data);
}

//
void ClawPatch::average_corner_ghost(const int& a_coarse_corner,
                                     const int& a_refratio,
                                     ClawPatch *cp_corner,
                                     fclaw_bool a_time_interp,
                                     fclaw2d_transform_data_t* transform_data)
{

    double *qcoarse = q_time_sync(a_time_interp);
    double *qfine = cp_corner->m_griddata.dataPtr();
    double *areacoarse = this->m_area.dataPtr();
    double *areafine = cp_corner->m_area.dataPtr();

    int manifold = m_manifold ? 1 : 0;
    average_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, a_refratio,
                          qcoarse, qfine,
                          areacoarse, areafine,
                          manifold,
                          a_coarse_corner,&transform_data);
}



void ClawPatch::interpolate_corner_ghost(const int& a_coarse_corner,
                                         const int& a_refratio,
                                         ClawPatch *cp_corner,
                                         fclaw_bool a_time_interp,
                                         fclaw2d_transform_data_t* transform_data)

{
    double *qcoarse = q_time_sync(a_time_interp);

    // qcorner is the finer level.
    double *qfine = cp_corner->m_griddata.dataPtr();

    interpolate_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                              a_refratio, qcoarse, qfine,
                              a_coarse_corner,&transform_data);
}


/* ----------------------------------------------------------------
   Special case : Pillow grid ghost exchanges/average/interpolate
   ---------------------------------------------------------------- */

void ClawPatch::mb_exchange_block_corner_ghost(const int& a_corner,
                                               ClawPatch *cp_corner)
{
    double *qthis = m_griddata.dataPtr();
    double *qcorner = cp_corner->m_griddata.dataPtr();


    mb_exchange_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner,
                                    a_corner, m_blockno);

}

// internal corners only a block boundaries.
void ClawPatch::mb_average_block_corner_ghost(const int& a_coarse_corner,
                                              const int& a_refratio,
                                              ClawPatch *cp_corner,
                                              fclaw_bool a_time_interp)
{
    // 'this' is the finer grid; 'cp_corner' is the coarser grid.
    double *qcoarse = q_time_sync(a_time_interp);


    double *areacoarse = this->m_area.dataPtr();
    double *areafine = cp_corner->m_area.dataPtr();
    double *qfine = cp_corner->m_griddata.dataPtr();

    mb_average_block_corner_ghost_(m_mx,m_my,m_mbc,m_meqn,
                                   a_refratio,qcoarse,qfine,
                                   areacoarse,areafine,
                                   a_coarse_corner,m_blockno);
}


void ClawPatch::mb_interpolate_block_corner_ghost(const int& a_coarse_corner,
                                                  const int& a_refratio,
                                                  ClawPatch *cp_corner,
                                                  fclaw_bool a_time_interp)

{
    double *qcoarse = q_time_sync(a_time_interp);

    /* qcorner is the finer level. */
    double *qfine = cp_corner->m_griddata.dataPtr();

    mb_interpolate_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                                       a_refratio, qcoarse, qfine,
                                       a_coarse_corner, m_blockno);
}



/* ----------------------------------------------------------------
   Phyisical boundary conditions
   ---------------------------------------------------------------- */

void ClawPatch::set_phys_corner_ghost(const int& a_corner, const int a_mthbc[],
                                      const double& t, const double& dt)
{
    double *q = m_griddata.dataPtr();

    /* This routine doesn't do anything ... */
    set_phys_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, q, a_corner, t, dt, a_mthbc);
}

void ClawPatch::exchange_phys_face_corner_ghost(const int& a_corner, const int& a_side,
                                                ClawPatch* cp)
{
    double *this_q = m_griddata.dataPtr();
    double *neighbor_q = cp->m_griddata.dataPtr();

    exchange_phys_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, this_q, neighbor_q,
                                a_corner, a_side);
}

/* ----------------------------------------------------------------
   Time interpolation
   ---------------------------------------------------------------- */
void ClawPatch::setup_for_time_interpolation(const double& alpha)
{
    /* Store time interpolated data that will be use in coarse grid
       exchanges */
    double *qlast = m_griddata_last.dataPtr();
    double *qcurr = m_griddata.dataPtr();
    double *qinterp = m_griddata_time_interpolated.dataPtr();
    int size = m_griddata.size();

    /* Copy everything first to fill in ghost cells of qinterp */
    memcpy(qinterp,qlast,size*sizeof(double));

    /* Do interpolation only on interior, since ghost cells in qcurr
       are invalid and will lead to floating point exceptions */
    TIMEINTERP_INTERIOR(&m_mx,&m_my,&m_mbc,&m_meqn,qcurr,qlast,qinterp,&alpha);
}


/* ----------------------------------------------------------------
   Mapped grids
   ---------------------------------------------------------------- */

void ClawPatch::set_block_corner_count(const int icorner, const int block_corner_count)
{
    m_block_corner_count[icorner] = block_corner_count;
}

void ClawPatch::setup_manifold(const int& level, const amr_options_t *gparms)
{
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int maxlevel = gparms->maxlevel;
    int refratio = gparms->refratio;

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 1;
    ur[1] = my + mbc + 1;

    Box box_p(ll,ur);   /* Store cell centered values here */

    /* Mesh cell centers of physical mesh */
    m_xp.define(box_p,1);
    m_yp.define(box_p,1);
    m_zp.define(box_p,1);
    m_surf_normals.define(box_p,3);
    m_curvature.define(box_p,3);

    // Compute area of the mesh cell.
    m_area.define(box_p,1);

    /* Node centered values */
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -mbc;
    }
    ur[0] = mx + mbc + 2;
    ur[1] = my + mbc + 2;
    Box box_d(ll,ur);

    m_xd.define(box_d,1);
    m_yd.define(box_d,1);
    m_zd.define(box_d,1);

    /* Face centered values */
    m_xface_normals.define(box_d,3);
    m_yface_normals.define(box_d,3);
    m_xface_tangents.define(box_d,3);
    m_yface_tangents.define(box_d,3);
    m_edge_lengths.define(box_d,2);


    /* Get pointers to pass to mesh routine */
    double *xp = m_xp.dataPtr();
    double *yp = m_yp.dataPtr();
    double *zp = m_zp.dataPtr();
    double *xd = m_xd.dataPtr();
    double *yd = m_yd.dataPtr();
    double *zd = m_zd.dataPtr();
    double *area = m_area.dataPtr();

    double *xnormals = m_xface_normals.dataPtr();
    double *ynormals = m_yface_normals.dataPtr();
    double *xtangents = m_xface_tangents.dataPtr();
    double *ytangents = m_yface_tangents.dataPtr();
    double *surfnormals = m_surf_normals.dataPtr();
    double *curvature = m_curvature.dataPtr();
    double *edge_lengths = m_edge_lengths.dataPtr();

    /* Compute centers and corners of mesh cell */
    setup_mesh_(mx,my,mbc,m_xlower,m_ylower,m_dx,m_dy,m_blockno,
                xp,yp,zp,xd,yd,zd);

    /* The level and the refratio is needed here to compute
       areas on coarser meshes based on areas of the finest
       level meshes. */
    compute_area_(mx, my, mbc, m_dx, m_dy,m_xlower, m_ylower,
                  m_blockno, area, level, maxlevel, refratio);

    compute_normals_(mx,my,mbc,xp,yp,zp,xd,yd,zd,
                     xnormals,ynormals);

    compute_tangents_(mx,my,mbc,xd,yd,zd,xtangents,ytangents,edge_lengths);

    compute_surf_normals_(mx,my,mbc,xnormals,ynormals,edge_lengths,
                          curvature, surfnormals,area);
}


/* ----------------------------------------------------------------
   Output and diagnostics
   ---------------------------------------------------------------- */

int ClawPatch::size()
{
    /* Use this to create new data */
    return m_griddata.size();
}
