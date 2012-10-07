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


void ClawPatch::define(const Real&  a_xlower,
                       const Real&  a_ylower,
                       const Real&  a_xupper,
                       const Real&  a_yupper,
                       const int& a_blockno,
                       const global_parms* a_gparms)
{
    m_mx = a_gparms->m_mx_leaf;
    m_my = a_gparms->m_my_leaf;
    m_mbc = a_gparms->m_mbc;
    m_blockno = a_blockno;

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
    m_griddata_last.define(box, m_meqn);
    m_griddata_save.define(box, m_meqn);
    m_griddata_time_interp.define(box, m_meqn);
    if (m_maux > 0)
    {
        m_auxarray.define(box,m_maux);
    }


    m_mapped = a_gparms->m_mapped;
    m_manifold = a_gparms->m_manifold;

    m_isDefined = true;
}

void ClawPatch::copyFrom(ClawPatch *a_cp)
{
    m_griddata = a_cp->m_griddata;
    /*
    // This should not be needed, as these various states will all be
    // recreated in the next time step.
    m_griddata_last = a_cp->m_griddata_last;
    m_griddata_save = a_cp->m_griddata_save;
    m_griddata_time_interp = a_cp->m_griddata_time_interp;
    */
}


bool ClawPatch::isDefined()
{
    return m_isDefined;
}

// ----------------------------------------------------------------
// Time stepping routines
// ----------------------------------------------------------------

void ClawPatch::setup_patch(const int& a_level, const int& a_maxlevel,
                            const int& a_refratio)
{
    if (m_maux > 0)
    {
        if (m_manifold)
        {
            setup_manifold(a_level, a_maxlevel, a_refratio);
        }
        setAuxArray();
    }
}


void ClawPatch::initialize()
{
    Real* q = m_griddata.dataPtr();
    Real* aux = m_auxarray.dataPtr();

    set_block_(&m_blockno);

    if (m_manifold)
    {
        qinit_mapped_(m_mx, m_my, m_meqn, m_mbc,m_xlower, m_ylower, m_dx, m_dy,
                      m_xp.dataPtr(), m_yp.dataPtr(), m_zp.dataPtr(),
                      q, m_maux, aux, m_blockno);
    }
    else
    {
        qinit_(m_mx,m_my,m_meqn,m_mbc,m_mx,m_my,m_xlower,m_ylower,
               m_dx,m_dy,q,m_maux,aux);
    }
}

void ClawPatch::setAuxArray()
{

    set_block_(&m_blockno);
    Real* aux = m_auxarray.dataPtr();

    if (m_manifold)
    {
        setaux_mapped_(m_mx,m_my,m_mbc,m_dx,m_dy,
                       m_xp.dataPtr(),m_yp.dataPtr(),m_zp.dataPtr(),
                       m_xd.dataPtr(),m_yd.dataPtr(),m_zd.dataPtr(),
                       m_area.dataPtr(), m_maux,aux);
    }
    else
    {
        setaux_(m_mx,m_my,m_mbc,m_mx,m_my,m_xlower,m_ylower,
                m_dx,m_dy,m_maux,aux);
    }
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
    set_block_(&m_blockno);

    Real* qold = m_griddata.dataPtr();
    Real* aux = m_auxarray.dataPtr();

    // Mysterious bug in this call when mx=my=4 and levels are (2,3).
    // m_griddata_last.dataPtr() == NULL, for some reason.
    m_griddata_last = m_griddata; // Copy for time interpolation

    // We also call a 'b4step2' in clawpath2, below.  But it won't
    // do anything in the mapped case.
    if (m_manifold)
    {

        b4step2_mapped_(m_mx,m_my, m_mbc,m_meqn,qold, m_dx,m_dy,
                        m_xp.dataPtr(), m_yp.dataPtr(), m_zp.dataPtr(),
                        m_xd.dataPtr(), m_yd.dataPtr(), m_zd.dataPtr(),
                        a_time, a_dt, m_maux, aux);

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
  int mwork = (maxm+2*gparms.m_mbc)*(12*gparms.m_meqn +
  (gparms.m_meqn+1)*gparms.m_mwaves + 3*gparms.m_maux + 2);
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

  delete [] fp;
  delete [] fm;
  delete [] gp;
  delete [] gm;

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
    set_block_(&m_blockno);
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

void ClawPatch::mb_exchange_face_ghost(const int& a_iface, ClawPatch *neighbor_cp)
{
    Real *qthis = m_griddata.dataPtr();
    Real *qneighbor = neighbor_cp->m_griddata.dataPtr();
    mb_exchange_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qthis,qneighbor,a_iface,m_blockno);
}

void ClawPatch::exchange_corner_ghost(const int& a_corner, ClawPatch *cp_corner)
{
    Real *qthis = m_griddata.dataPtr();
    Real *qcorner = cp_corner->m_griddata.dataPtr();

    exchange_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner, a_corner);

}

void ClawPatch::mb_exchange_corner_ghost(const int& a_corner, bool a_intersects_block[],
                                         ClawPatch *cp_corner, const bool& a_is_block_corner)
{
    Real *qthis = m_griddata.dataPtr();
    Real *qcorner = cp_corner->m_griddata.dataPtr();

    if (a_is_block_corner)
    {
        // We know we are at a block corner, which is handled differently than a corner that is
        // only at an edge, but not at a corner.
        mb_exchange_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner,
                                        a_corner, m_blockno);
    }
    else
    {
        int numfaces = 2*SpaceDim;
        int bdry[numfaces];
        for(int m = 0; m < numfaces; m++)
        {
            bdry[m] = a_intersects_block[m] ? 1 : 0;
        }
        mb_exchange_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, qthis, qcorner,
                                  a_corner, bdry, m_blockno);

    }
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
                                   bool a_time_interp,
                                   bool a_block_boundary)
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
            if (a_block_boundary)
            {
                mb_average_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,
                                       auxcoarse, auxfine, m_maux,
                                       a_idir,a_iface_coarse,
                                       a_p4est_refineFactor,a_refratio,igrid);
            }
            else
            {
                average_face_ghost_mapped_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,
                                           auxcoarse, auxfine, m_maux,
                                           a_idir,a_iface_coarse,
                                           a_p4est_refineFactor,a_refratio,igrid);
            }
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
                                       bool a_time_interp,
                                       bool a_block_boundary)
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
        if (a_block_boundary)
        {
            mb_interpolate_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iside,
                                              a_p4est_refineFactor,a_refratio,igrid);
        }
        else
        {
            interpolate_face_ghost_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,a_idir,a_iside,
                                    a_p4est_refineFactor,a_refratio,igrid);
        }
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


// internal corners only a block boundaries.
void ClawPatch::mb_average_corner_ghost(const int& a_coarse_corner, const int& a_refratio,
                                        ClawPatch *cp_corner, bool a_time_interp,bool is_block_corner,
                                        bool intersects_block[])
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

    Real *auxcoarse = this->m_auxarray.dataPtr();
    Real *auxfine = cp_corner->m_auxarray.dataPtr();
    Real *qfine = cp_corner->m_griddata.dataPtr();

    if (is_block_corner)
    {
        mb_average_block_corner_ghost_(m_mx,m_my,m_mbc,m_meqn,a_refratio,qcoarse,qfine,
                                       auxcoarse,auxfine,m_maux,a_coarse_corner,m_blockno);
    }
    else
    {
        int block_bdry[4];
        for (int m = 0; m < 4; m++)
        {
            block_bdry[m] = intersects_block[m] ? 1 : 0;
        }
        mb_average_corner_ghost_(m_mx, m_my, m_mbc, m_meqn, a_refratio, qcoarse, qfine,
                                 auxcoarse, auxfine, m_maux, a_coarse_corner, block_bdry);
    }
}


void ClawPatch::mb_interpolate_corner_ghost(const int& a_coarse_corner,
                                            const int& a_refratio,
                                            ClawPatch *cp_corner,
                                            bool a_time_interp, bool is_block_corner,
                                            bool intersects_block[])

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

    if (is_block_corner)
    {
        // This doesn't do anything right now.
        mb_interpolate_block_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                                           a_refratio, qcoarse, qfine,
                                           a_coarse_corner, m_blockno);
    }
    else
    {
        int bdry[4];
        for(int m = 0; m < 4; m++)
        {
            bdry[m] = intersects_block[m] ? 1 : 0;
        }
        mb_interpolate_corner_ghost_(m_mx, m_my, m_mbc, m_meqn,
                                     a_refratio, qcoarse, qfine,
                                     a_coarse_corner, bdry);
    }

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

void ClawPatch::coarsen_from_fine_family(ClawPatch *a_cp_siblings[],
                                         const int& a_refratio,
                                         const int& a_num_siblings,
                                         const int& a_p4est_refineFactor)
{
    Real *qcoarse = m_griddata.dataPtr();
    for(int igrid = 0; igrid < a_num_siblings; igrid++)
    {
        Real *qfine = a_cp_siblings[igrid]->m_griddata.dataPtr();
        if (m_manifold)
        {
            Real *auxcoarse = m_auxarray.dataPtr();
            Real *auxfine = a_cp_siblings[igrid]->m_auxarray.dataPtr();
            average_to_coarse_mapped_(m_mx, m_my, m_mbc, m_meqn, qcoarse, qfine,
                                      auxcoarse, auxfine, m_maux,
                                      a_p4est_refineFactor,
                                      a_refratio, igrid);
        }
        else
        {
            average_to_coarse_patch_(m_mx,m_my,m_mbc,m_meqn,qcoarse,qfine,
                                     a_p4est_refineFactor,a_refratio,igrid);
        }
    }
}

bool ClawPatch::tag_for_refinement(bool a_init_flag)
{
    Real *q = m_griddata.dataPtr();
    int tag_patch;  // == 0 or 1
    int iflag = a_init_flag ? 1 : 0;
    tag_for_refinement_(m_mx,m_my,m_mbc,m_meqn,m_xlower,m_ylower,
                        m_dx, m_dy,q,iflag,tag_patch);
    return tag_patch == 1;
}

bool ClawPatch::tag_for_coarsening(ClawPatch *a_cp_siblings[],
                                   const int& a_refratio,
                                   const int& a_num_siblings,
                                   const int& a_p4est_refineFactor)
{
    this->coarsen_from_fine_family(a_cp_siblings,a_refratio,a_num_siblings,
                                   a_p4est_refineFactor);
    int tag_patch;
    Real *qcoarse = m_griddata.dataPtr();
    tag_for_coarsening_(m_mx,m_my,m_mbc,m_meqn,m_xlower,m_ylower,m_dx,m_dy,
                              qcoarse,tag_patch);
    return tag_patch == 0;
}


// ----------------------------------------------------------------
// Mapped grids
// ----------------------------------------------------------------

void ClawPatch::setup_manifold(const int& a_level,
                               const int& a_maxlevel, const int& a_refratio)
{
    // Set fortran common block
    set_block_(&m_blockno);

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
    write_qfile_(m_mx,m_my,m_meqn,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,q,
                 a_iframe,a_patch_num,a_level,m_blockno);
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

void ClawPatch::dump_auxarray()
{
    Real *q;
    q = m_auxarray.dataPtr();
    int k = 0;
    for (int m = 0; m < m_maux; m++)
    {
        for(int j = 1-m_mbc; j <= m_my+m_mbc; j++)
        {
            for(int i = 1-m_mbc; i <= m_mx+m_mbc; i++)
            {
                printf("q[%2d,%2d,%2d] = %24.16e\n",i,j,m,q[k]);
                k++;
            }
            printf("\n");
        }
        printf("\n");
        printf("\n");
    }
}
