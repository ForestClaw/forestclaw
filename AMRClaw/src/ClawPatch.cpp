#include "ClawPatch.H"
#include "amr_utils.H"

// This constructors includes all of parameters that are patch independent.
// All of this could also be in some sort of "set_params" function...
ClawPatch::ClawPatch()
{
    m_isDefined = false;
}

ClawPatch::~ClawPatch()
{
    // delete m_corners_set;
}


#if CH_SPACEDIM == 2
void ClawPatch::define(const Real&  a_xlower,
                       const Real&  a_ylower,
                       const Real&  a_xupper,
                       const Real&  a_yupper,
                       const global_parms* a_gparms)
{
    m_mx = a_gparms->m_mx_leaf;
    m_my = a_gparms->m_my_leaf;
    m_mbc = a_gparms->m_mbc;

    m_xlower = a_xlower;
    m_ylower = a_ylower;
    m_xupper = a_xupper;
    m_yupper = a_yupper;

    m_dx = (a_xupper - a_xlower)/m_mx;
    m_dy = (a_yupper - a_ylower)/m_my;

    m_meqn = a_gparms->m_meqn;
    m_maux = a_gparms->m_maux;

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
    m_griddata_time_interp.define(box, m_meqn);
    if (m_maux > 0)
    {
        m_auxarray.define(box,m_maux);
    }

    m_mapped = a_gparms->m_mapped;
    m_manifold = a_gparms->m_manifold;

    m_isDefined = true;
}
#else
void ClawPatch::define(const Real& a_xlower,
                       const Real& a_ylower,
                       const Real& a_zlower,
                       const Real& a_xupper,
                       const Real& a_yupper,
                       const Real& a_zupper,
                       const global_parms* a_gparms)
{
    m_mx = a_gparms->m_mx_leaf;
    m_my = a_gparms->m_my_leaf;
    m_mz = a_gparms->m_mz_leaf;
    m_mbc = a_gparms->m_mbc;

    m_xlower = a_xlower;
    m_ylower = a_ylower;
    m_xupper = a_xupper;
    m_yupper = a_yupper;
    m_zlower = a_zlower;
    m_zupper = a_zupper;

    m_dx = (a_xupper - a_xlower)/m_mx;
    m_dy = (a_yupper - a_ylower)/m_my;
    m_dz = (a_zupper - a_zlower)/m_mz;


    // Set box for grid data.  Use local indexing for now.
    // Note that box doesn't know anything about the number of boundary conditions
    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-m_mbc;
    }
    ur[0] = m_mx + m_mbc;
    ur[1] = m_my + m_mbc;
    ur[2] = m_mz + m_mbc;
    Box box(ll,ur);

    m_meqn = a_gparms->m_meqn;
    m_maux = a_gparms->m_maux;
    m_mapped = a_gparms->m_mapped;
    m_manifold = a_gparms->m_manifold;

    m_griddata.define(box, m_meqn);
    m_griddata_last.define(box,m_meqn);
    if (m_maux > 0)
    {
        m_auxarray.define(box,m_maux);
    }

    m_isDefined = true;
}
#endif

void ClawPatch::copyFrom(ClawPatch *a_cp)
{
    m_mx = a_cp->m_mx;
    m_my = a_cp->m_my;
    m_mbc = a_cp->m_mbc;
    m_meqn = a_cp->m_meqn;
    m_maux = a_cp->m_maux;
    m_mapped = a_cp->m_mapped;
    m_manifold = a_cp->m_manifold;

    m_xlower = a_cp->m_xlower;
    m_ylower = a_cp->m_ylower;
    m_xupper = a_cp->m_xupper;
    m_yupper = a_cp->m_yupper;

    m_dx = a_cp->m_dx;
    m_dy = a_cp->m_dy;

    m_griddata = a_cp->m_griddata;
    m_auxarray = a_cp->m_auxarray;

    Box box = m_griddata.box();
    m_griddata_time_interp.define(box,m_meqn);

    // m_griddata_save and m_griddata_last will get allocated when we set them equal to
    // current time steps.
}


bool ClawPatch::isDefined()
{
    return m_isDefined;
}

// ----------------------------------------------------------------
// Time stepping routines
// ----------------------------------------------------------------

void ClawPatch::setup_patch(const int& level, const int& maxlevel, const int& refratio)
{
    if (m_maux > 0)
    {
        if (m_manifold)
        {
            setup_manifold(level, maxlevel, refratio);
        }
        setAuxArray();
    }
    // Anything else to do?
}


void ClawPatch::initialize()
{
    Real* q = m_griddata.dataPtr();
    Real* aux = m_auxarray.dataPtr();

#if CH_SPACEDIM == 2
    if (m_manifold)
    {
        qinit_mapped_(m_mx, m_my, m_meqn, m_mbc,
               m_xp.dataPtr(), m_yp.dataPtr(), m_zp.dataPtr(), q, m_maux, aux);
    }
    else
    {
        qinit_(m_mx,m_my,m_meqn,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,q,m_maux,aux);
    }
#elif CH_SPACEDIM == 3
    qinit_(mx,my,mz,m_meqn,m_mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,m_maux,aux);
#endif
}

void ClawPatch::setAuxArray()
{
    Real* aux = m_auxarray.dataPtr();

#if CH_SPACEDIM == 2
    if (m_manifold)
    {
        setaux_mapped_(m_mx,m_my,m_mbc,m_dx,m_dy,
                       m_xp.dataPtr(),m_yp.dataPtr(),m_zp.dataPtr(),
                       m_xd.dataPtr(),m_yd.dataPtr(),m_zd.dataPtr(),
                       m_area.dataPtr(), m_maux,aux);
    }
    else
    {
        setaux_(m_mx,m_my,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,m_maux,aux);
    }
#elif CH_SPACEDIM == 3
    setaux_(m_mx,m_my,m_mz,m_mbc,m_mx,m_my,m_mz,m_xlower,m_ylower,m_zlower,m_dx,m_dy,m_dz,m_maux,aux);
#endif
}

Real ClawPatch::step(const Real& a_time,
                     const Real& a_dt,
                     const int& a_level,
                     const global_parms& gparms)
{
    Real maxwavespeed = 1; // Making this up...
    Real cfl_grid = a_dt/m_dx*maxwavespeed; //
    return cfl_grid;
}

#if CH_SPACEDIM == 2

Real ClawPatch::step_noqad(const Real& a_time,
                           const Real& a_dt,
                           const int& a_level,
                           const global_parms& gparms)
{
    Real* qold = m_griddata.dataPtr();
    Real* aux = m_auxarray.dataPtr();

    // Mysterious bug in this call when mx=my=4 and levels are (2,3).
    // m_griddata_last.dataPtr() == NULL, for some reason.
    m_griddata_last = m_griddata; // Copy for time interpolation

    // We also call a 'b4step2' in clawpath2, below.  But it won't
    // do anything in the mapped case.
    if (m_manifold)
    {
        /*
        b4step2_mapped_(m_mx,m_my, m_mbc,m_meqn,qold, m_dx,m_dy,
                        m_xp.dataPtr(), m_yp.dataPtr(), m_zp.dataPtr(),
                        m_xd.dataPtr(), m_yd.dataPtr(), m_zd.dataPtr(),
                        a_time, a_dt, m_maux, aux);
        */
    }

    int maxm = max(m_mx,m_my);

    Real cflgrid;

    int mwork = (maxm+2*m_mbc)*(12*m_meqn + (m_meqn+1)*gparms.m_mwaves + 3*m_maux + 2);
    Real* work = new Real[mwork];

    Real* fp = new Real[m_meqn*(m_mx+2*m_mbc)*(m_my+2*m_mbc)];
    Real* fm = new Real[m_meqn*(m_mx+2*m_mbc)*(m_my+2*m_mbc)];
    Real* gp = new Real[m_meqn*(m_mx+2*m_mbc)*(m_my+2*m_mbc)];
    Real* gm = new Real[m_meqn*(m_mx+2*m_mbc)*(m_my+2*m_mbc)];

    clawpatch2_(maxm, m_meqn, m_maux, m_mbc, gparms.m_method,
                gparms.m_mthlim, gparms.m_mcapa, gparms.m_mwaves, m_mx, m_my, qold,
                aux, m_dx, m_dy, a_dt, cflgrid, work, mwork, m_xlower, m_ylower,a_level,
                a_time, fp, fm, gp, gm);

    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}
#endif


Real ClawPatch::ClawPatchIntegrator(const Real& a_time,
                                    const Real& a_dt,
                                    const int& a_refRatio,
                                    const int& a_level,
                                    const global_parms& gparms)
{

    // Real dt = a_dt;

  // Data for step2 or step3.  This is overwritten by updated values.
  // Real* qold = m_griddata.dataPtr();
  // Real* aux = m_auxarray.dataPtr();

  // int maxm = max(m_mx,m_my);
#if CH_SPACEDIM == 3
  maxm = max(maxm,mz);
#endif

  // set common block for level
  set_common_levels_(gparms.m_maxlevel,a_level,gparms.m_refratio);


  Real cflgrid = 1;


  /*
#if CH_SPACEDIM==2
  int mwork = (maxm+2*gparms.m_mbc)*(12*gparms.m_meqn + (gparms.m_meqn+1)*gparms.m_mwaves + 3*gparms.m_maux + 2);
  Real* work = new Real[mwork];

  Real* fp = new Real[gparms.m_meqn*(m_mx+2*gparms.m_mbc)*(m_my+2*gparms.m_mbc)];
  Real* fm = new Real[gparms.m_meqn*(m_mx+2*gparms.m_mbc)*(m_my+2*gparms.m_mbc)];
  Real* gp = new Real[gparms.m_meqn*(m_mx+2*gparms.m_mbc)*(m_my+2*gparms.m_mbc)];
  Real* gm = new Real[gparms.m_meqn*(m_mx+2*gparms.m_mbc)*(m_my+2*gparms.m_mbc)];

  Real* fp_chombo = a_fluxp[0].dataPtr();
  Real* fm_chombo = a_fluxm[0].dataPtr();
  Real* fpc_chombo = a_fluxpc[0].dataPtr();
  Real* fmc_chombo = a_fluxmc[0].dataPtr();

  Real* gp_chombo = a_fluxp[1].dataPtr();
  Real* gm_chombo = a_fluxm[1].dataPtr();
  Real* gpc_chombo = a_fluxpc[1].dataPtr();
  Real* gmc_chombo = a_fluxmc[1].dataPtr();

  Real* qadd_x = a_qadd[0].dataPtr();
  Real* qadd_y = a_qadd[1].dataPtr();

  int m_auxtype_int[10];  // dummy for now;  fix!
  clawpatch2_(maxm, gparms.m_meqn, gparms.m_maux, gparms.m_mbc, gparms.m_method,
              gparms.m_mthlim,
              gparms.m_mcapa, gparms.m_mwaves, m_mx, m_my, qold, aux,
              m_dx, m_dy, dt, cflgrid, work, mwork, qold_coarse, auxold_coarse,
              qadd_x, qadd_y, m_auxtype_int, m_xlower, m_ylower,
              intersectsBoundary,a_level,
              gparms.m_mthbc, a_time, mxc, myc, fp, fm, gp, gm,
              fp_chombo,fm_chombo,gp_chombo,gm_chombo,
              fpc_chombo,fmc_chombo,gpc_chombo,gmc_chombo);

#elif CH_SPACEDIM==3
  int mwork = (maxm+2*m_mbc)*(46*m_meqn + (m_meqn+1)*m_mwaves + 9*m_maux + 3);
  Real* work = new Real[mwork];

  Real* fp_chombo = a_fluxp[0].dataPtr();
  Real* fm_chombo = a_fluxm[0].dataPtr();
  Real* fpc_chombo = a_fluxpc[0].dataPtr();
  Real* fmc_chombo = a_fluxmc[0].dataPtr();

  Real* gp_chombo = a_fluxp[1].dataPtr();
  Real* gm_chombo = a_fluxm[1].dataPtr();
  Real* gpc_chombo = a_fluxpc[1].dataPtr();
  Real* gmc_chombo = a_fluxmc[1].dataPtr();

  Real* hp_chombo = a_fluxp[2].dataPtr();
  Real* hm_chombo = a_fluxm[2].dataPtr();
  Real* hpc_chombo = a_fluxpc[2].dataPtr();
  Real* hmc_chombo = a_fluxmc[2].dataPtr();

  Real* fp = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];
  Real* fm = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];

  Real* gp = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];
  Real* gm = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];

  Real* hp = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];
  Real* hm = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)*(mz+2*m_mbc)];

  Real* qadd_x = a_qadd[0].dataPtr();
  Real* qadd_y = a_qadd[1].dataPtr();
  Real* qadd_z = a_qadd[2].dataPtr();

  clawpatch3_(maxm, m_meqn, m_maux, m_mbc, m_method, m_mthlim,
              m_mcapa, m_mwaves, mx, my, mz, qold, aux,
              dx, dy, dz, dt, cflgrid, work, mwork, qold_coarse, auxold_coarse,
              qadd_x, qadd_y, qadd_z, m_auxtype_int, xlower, ylower, zlower,
              intersectsBoundary,a_level,
              m_mthbc, a_time, mxc, myc, mzc, fp, fm, gp, gm, hp, hm,
              fp_chombo,fm_chombo,gp_chombo,gm_chombo,hp_chombo,hm_chombo,
              fpc_chombo, fmc_chombo, gpc_chombo, gmc_chombo,hpc_chombo,hmc_chombo);

#endif

  delete [] fp;
  delete [] fm;
  delete [] gp;
  delete [] gm;
#if CH_SPACEDIM == 3
  delete [] hp;
  delete [] hm;
#endif

  delete [] intersectsBoundary;
  delete [] work;

  // Real maxWaveSpeed = (dx/dt)*cflgrid;
  // return maxWaveSpeed;
  */

  return cflgrid;
}


void ClawPatch::save_step()
{
    // Store a backup in case the CFL number is too large doesn't work out.
    m_griddata_save = m_griddata;
}

void ClawPatch::restore_step()
{
    m_griddata = m_griddata_save;
}


void ClawPatch::time_interpolate(const int& a_fine_step, const int& a_coarse_step,
                                 const int& a_refratio)
{
    Real alpha = Real(a_fine_step)/Real(a_refratio);

    Real *qlast = m_griddata_last.dataPtr();
    Real *qcurr = m_griddata.dataPtr();
    Real *qtimeinterp = m_griddata_time_interp.dataPtr();
    int size = m_griddata.size();

    // This works, even when we have a system (meqn > 1).  Note that all ghost values
    // will be interpolated.
    for(int i = 0; i < size; i++)
    {
        // There is surely a BLAS routine that does this...
        qtimeinterp[i] = qlast[i] + alpha*(qcurr[i] - qlast[i]);
    }
}


// ----------------------------------------------------------------
// Single level exchanges
// ----------------------------------------------------------------

void ClawPatch::exchange_face_ghost(const int& a_idir, ClawPatch *neighbor_cp)
{
    Real *qthis = m_griddata.dataPtr();
    Real *qneighbor = neighbor_cp->m_griddata.dataPtr();
    exchange_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qthis,qneighbor,a_idir);
}

void ClawPatch::exchange_corner_ghost(const int& a_corner, ClawPatch *cp_corner)
{
    Real *qthis = m_griddata.dataPtr();
    Real *qcorner = cp_corner->m_griddata.dataPtr();

    exchange_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner, a_corner);

}

void ClawPatch::set_phys_face_ghost(const bool a_intersects_bc[], const int a_mthbc[],
                                    const Real& t, const Real& dt)
{
    Real *q = m_griddata.dataPtr();
    Real *aux = m_auxarray.dataPtr();

    // Set a local copy of mthbc that can be used for a patch.
    int mthbc[2*SpaceDim];
    for(int i = 0; i < 2*SpaceDim; i++)
    {
        if (a_intersects_bc[i])
        {
            mthbc[i] = a_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }
    bc2_(m_mx,m_my,m_meqn,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,q,m_maux,aux,t,dt,mthbc);
}


void ClawPatch::set_phys_corner_ghost(const int& a_corner, const int a_mthbc[],
                                      const Real& t, const Real& dt)
{
    Real *q = m_griddata.dataPtr();

    // No code yet
    set_phys_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, q, a_corner, t, dt, a_mthbc);
}

void ClawPatch::exchange_phys_face_corner_ghost(const int& a_corner, const int& a_side,
                                                ClawPatch* cp)
{
    Real *this_q = m_griddata.dataPtr();
    Real *neighbor_q = cp->m_griddata.dataPtr();

    exchange_phys_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, this_q, neighbor_q,
                                a_corner, a_side);
}


// ----------------------------------------------------------------
// Multi-level operations
// ----------------------------------------------------------------
void ClawPatch::average_face_ghost(const int& a_idir,
                                   const int& a_iface_coarse,
                                   const int& a_p4est_refineFactor,
                                   const int& a_refratio,
                                   ClawPatch **neighbor_cp,
                                   bool a_time_interp)
{
    Real *qcoarse;
    if (a_time_interp)
    {
        qcoarse = m_griddata_time_interp.dataPtr();
    }
    else
    {
        qcoarse = m_griddata.dataPtr();
    }
    for(int igrid = 0; igrid < a_p4est_refineFactor; igrid++)
    {
        Real *qfine = neighbor_cp[igrid]->m_griddata.dataPtr();
        if (m_manifold)
        {
            Real *auxcoarse = m_auxarray.dataPtr();
            Real *auxfine = neighbor_cp[igrid]->m_auxarray.dataPtr();
            average_face_ghost_mapped_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,
                                       auxcoarse, auxfine, m_maux,
                                       a_idir,a_iface_coarse,
                                       a_p4est_refineFactor,a_refratio,igrid);
        }
        else
        {
            average_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iface_coarse,
                                a_p4est_refineFactor,a_refratio,igrid);
        }
    }
}

void ClawPatch::interpolate_face_ghost(const int& a_idir,
                                       const int& a_iside,
                                       const int& a_p4est_refineFactor,
                                       const int& a_refratio,
                                       ClawPatch **neighbor_cp,
                                       bool a_time_interp)
{
    Real *qcoarse;
    if (a_time_interp)
    {
        qcoarse = m_griddata_time_interp.dataPtr();
    }
    else
    {
        qcoarse = m_griddata.dataPtr();
    }

    for(int ir = 0; ir < a_p4est_refineFactor; ir++)
    {
        Real *qfine = neighbor_cp[ir]->m_griddata.dataPtr();
        int igrid = ir; // indicates which grid we are averaging from.
        interpolate_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iside,
                                a_p4est_refineFactor,a_refratio,igrid);
    }
}

//
void ClawPatch::average_corner_ghost(const int& a_coarse_corner, const int& a_refratio,
                                     ClawPatch *cp_corner, bool a_time_interp)
{
    // 'this' is the finer grid; 'cp_corner' is the coarser grid.
    Real *qcoarse;
    if (a_time_interp)
    {
        qcoarse = m_griddata_time_interp.dataPtr();
    }
    else
    {
        qcoarse = m_griddata.dataPtr();
    }

    Real *qfine = cp_corner->m_griddata.dataPtr();

    average_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, a_refratio, qcoarse, qfine, a_coarse_corner);

}

void ClawPatch::interpolate_corner_ghost(const int& a_coarse_corner, const int& a_refratio,
                                         ClawPatch *cp_corner, bool a_time_interp)
{
    Real *qcoarse;
    if (a_time_interp)
    {
        qcoarse = m_griddata_time_interp.dataPtr();
    }
    else
    {
        qcoarse = m_griddata.dataPtr();
    }

    // qcorner is the finer level.
    Real *qfine = cp_corner->m_griddata.dataPtr();

    interpolate_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, a_refratio, qcoarse, qfine, a_coarse_corner);
}


// ----------------------------------------------------------------
// Tagging, refining and coarsening
// ----------------------------------------------------------------

void ClawPatch::interpolate_to_fine_patch(ClawPatch* a_fine,
                                          const int& a_igrid,
                                          const int& a_p4est_refineFactor,
                                          const int& a_refratio)
{
    Real *qcoarse = this->m_griddata.dataPtr();
    Real *qfine = a_fine->m_griddata.dataPtr();

    // Use linear interpolation with limiters.
    interpolate_to_fine_patch_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_p4est_refineFactor,
                               a_refratio,a_igrid);
    if (m_manifold)
    {
        /* Doesn't quite work
        Real *auxcoarse = m_griddata.dataPtr();
        Real *auxfine = a_fine->m_griddata.dataPtr();
        fixcapaq2_(m_mx, m_my, m_mbc, m_meqn,qcoarse, qfine, auxcoarse, auxfine,
                   m_maux, a_p4est_refineFactor, a_refratio, a_igrid);
        */
    }
}

void ClawPatch::coarsen_from_fine_patch(ClawPatch * a_fine, const int& a_igrid,
                                   const int& a_p4est_refineFactor, const int& a_refratio)
{
    Real *qfine = a_fine->m_griddata.dataPtr();
    Real *qcoarse = this->m_griddata.dataPtr();

    if (m_manifold)
    {
        Real *auxcoarse = m_auxarray.dataPtr();
        Real *auxfine = a_fine->m_auxarray.dataPtr();
        average_to_coarse_mapped_(m_mx, m_my, m_mbc, m_meqn, qcoarse, qfine,
                                  auxcoarse, auxfine, m_maux, p4est_refineFactor,
                                  a_refratio, a_igrid);
    }
    else
    {
        average_to_coarse_patch_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,p4est_refineFactor,a_refratio,a_igrid);
    }
}

bool ClawPatch::tag_for_refinement(bool a_init_flag)
{
    Real *q = m_griddata.dataPtr();
    int tag_patch;  // == 0 or 1
    int iflag;

    // I am sure there is a way to pass a boolean to Fortran ...
    if (a_init_flag)
    {
        iflag = 1;
    }
    else
    {
        iflag = 0;
    }
    tag_for_refinement_(m_mx,m_my,m_mbc,m_meqn,m_xlower,m_ylower,m_dx, m_dy,q,iflag,tag_patch);
    return tag_patch == 1;
}

bool ClawPatch::tag_for_coarsening(const int& a_refratio)
{
    Real *q = m_griddata.dataPtr();
    int tag_patch;
    int mx2 = m_mx/a_refratio;
    int my2 = m_my/a_refratio;
    Real *qcoarse = new Real[(mx2 + 2*m_mbc)*(my2 + 2*m_mbc)*m_meqn];
    coarsen_(m_mx,m_my,m_mbc,m_meqn,a_refratio, q, qcoarse);

    // it would be nice to know where the fine grid is in relation to the underlying coarse
    // quadrant.
    Real dx2 = a_refratio*m_dx;
    Real dy2 = a_refratio*m_dy;
    int iflag = 0;
    tag_for_refinement_(mx2,my2,m_mbc,m_meqn,m_xlower,m_ylower,dx2,dy2,qcoarse,iflag,tag_patch);
    delete [] qcoarse;
    return tag_patch == 0;
}


// ----------------------------------------------------------------
// Mapped grids
// ----------------------------------------------------------------

void ClawPatch::setup_manifold(const int& a_level, const int& a_maxlevel, const int& a_refratio)
{

    // This is a generic routine that sets all things related to the mapping.
    set_maptype_();
    // Do we really ever use the "box"?
    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -m_mbc;
    }
    ur[0] = m_mx + m_mbc + 1;
    ur[1] = m_my + m_mbc + 1;

    Box box_p(ll,ur);

    // Mesh cell centers of physical mesh
    m_xp.define(box_p,1);
    m_yp.define(box_p,1);
    m_zp.define(box_p,1);

    // Compute area of the mesh cell.
    m_area.define(box_p,1);

    // Mesh cell corners of physical mesh
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = -m_mbc;
    }
    ur[0] = m_mx + m_mbc + 2;
    ur[1] = m_my + m_mbc + 2;
    Box box_d(ll,ur);

    m_xd.define(box_d,1);
    m_yd.define(box_d,1);
    m_zd.define(box_d,1);

    // Compute centers and corners of mesh cell
    setup_mesh_(m_mx,m_my,m_mbc,m_xlower,m_ylower,m_dx,m_dy,
                m_xp.dataPtr(),m_yp.dataPtr(),m_zp.dataPtr(),
                m_xd.dataPtr(),m_yd.dataPtr(),m_zd.dataPtr());


    compute_area_(m_mx, m_my, m_mbc, m_dx, m_dy,m_xlower, m_ylower,
                  m_area.dataPtr(), a_level, a_maxlevel, a_refratio);
}


// ----------------------------------------------------------------
// Output and diagnostics
// ----------------------------------------------------------------


void ClawPatch::write_patch_data(const int& a_iframe, const int& a_patch_num, const int& a_level)
{
    Real *q = m_griddata.dataPtr();
    write_qfile_(m_mx,m_my,m_meqn,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,q,a_iframe,a_patch_num,a_level);
}

Real ClawPatch::compute_sum()
{
    Real *q = m_griddata.dataPtr();
    Real sum;
    compute_sum_(m_mx,m_my,m_mbc,m_meqn,m_dx, m_dy, q,sum);
    return sum;
}

void ClawPatch::dump()
{
    Real *q;
    q = m_griddata.dataPtr();
    int k = 0;
    for(int j = 1-m_mbc; j <= m_my+m_mbc; j++)
    {
        for(int i = 1-m_mbc; i <= m_mx+m_mbc; i++)
        {
            printf("q[%2d,%2d] = %24.16e\n",i,j,q[k]);
            k++;
        }
        printf("\n");
    }
}

void ClawPatch::dump_last()
{
    Real *q;
    q = m_griddata_last.dataPtr();
    int k = 0;
    for(int j = 1-m_mbc; j <= m_my+m_mbc; j++)
    {
        for(int i = 1-m_mbc; i <= m_mx+m_mbc; i++)
        {
            printf("q[%2d,%2d] = %24.16e\n",i,j,q[k]);
            k++;
        }
        printf("\n");
    }
}

void ClawPatch::dump_time_interp()
{
    Real *q;
    q = m_griddata_time_interp.dataPtr();
    int k = 0;
    for(int j = 1-m_mbc; j <= m_my+m_mbc; j++)
    {
        for(int i = 1-m_mbc; i <= m_mx+m_mbc; i++)
        {
            printf("q[%2d,%2d] = %24.16e\n",i,j,q[k]);
            k++;
        }
        printf("\n");
    }
}
