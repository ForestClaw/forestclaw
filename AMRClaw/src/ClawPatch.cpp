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

    int ll[SpaceDim];
    int ur[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        ll[idir] = 1-m_mbc;
    }
    ur[0] = m_mx + m_mbc;
    ur[1] = m_my + m_mbc;
    Box box(ll,ur);


    m_meqn = a_gparms->m_meqn;
    m_maux = a_gparms->m_maux;

    int gridsize = (m_mx+2*m_mbc)*(m_my+2*m_mbc);
    m_griddata.define(gridsize*m_meqn,box);
    if (m_maux > 0)
    {
        m_auxarray.define(gridsize*m_maux,box);
    }

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

    int gridsize = (m_mx+2*m_mbc)*(m_my+2*m_mbc)*(m_mz+2*m_mbc);
    m_griddata.define(gridsize*m_meqn,box);
    if (m_maux > 0)
    {
        m_auxarray.define(gridsize*m_maux,box);
    }

    m_isDefined = true;
}
#endif

bool ClawPatch::isDefined()
{
    return m_isDefined;
}


int ClawPatch::get_mx() const
{
    return m_mx;
}

int ClawPatch::get_my() const
{
    return m_my;
}

#if CH_SPACEDIM == 3
int ClawPatch::get_mz() const
{
    return m_mz;
}
#endif

Real ClawPatch::get_xlower() const
{
    return m_xlower;
}

Real ClawPatch::get_ylower() const
{
    return m_ylower;
}

#if CH_SPACEDIM == 3
Real ClawPatch::get_zlower() const
{
    return m_zlower;
}
#endif


Real ClawPatch::get_xupper() const
{
    return m_xupper;
}

Real ClawPatch::get_yupper() const
{
    return m_yupper;
}

#if CH_SPACEDIM==3
Real ClawPatch::get_zupper() const
{
    return m_zupper;
}
#endif


void ClawPatch::initialize()
{
    CH_assert(m_isDefined);

    Real* q = m_griddata.dataPtr();
    Real* aux = m_auxarray.dataPtr();

    //  c     =====================================================
    //         subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
    //       &                   dx,dy,q,maux,aux)
    //  c     =====================================================

    // first two arguments are maxmx/maxmy
    // is it ok to send mx/my for that?  (ndk)
#if CH_SPACEDIM == 2
    qinit_(m_mx,
           m_my,
           m_meqn,
           m_mbc,
           m_mx,
           m_my,
           m_xlower,
           m_ylower,
           m_dx,
           m_dy,
           q,
           m_maux,
           aux);
#elif CH_SPACEDIM == 3
    qinit_(mx,my,mz,m_meqn,m_mbc,mx,my,mz,xlower,ylower,zlower,
           dx,dy,dz,q,m_maux,aux);
#endif
}

void ClawPatch::setAuxArray(const int& a_maxlevel,
                            const int& a_refratio,
                            const int& a_level)
{
    CH_assert(m_isDefined);

    // set common block for level
    set_common_levels_(a_maxlevel,a_level,a_refratio);

    Real* aux = m_auxarray.dataPtr();

    // Default setaux.f routine.  (2d version)
    // c     ============================================
    //       subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
    //      &                  maux,aux)
    // c     ============================================

    // first two arguments are maxmx/maxmy
    // is it ok to send mx/my for that?  (ndk)
    // yes - maxmx, maxmy only used for static dimensioning of arrays in Clawpack (DAC).
#if CH_SPACEDIM == 2
    setaux_(m_mx,m_my,m_mbc,m_mx,m_my,m_xlower,m_ylower,m_dx,m_dy,m_maux,aux);
#elif CH_SPACEDIM == 3
    setaux_(m_mx,m_my,m_mz,m_mbc,m_mx,m_my,m_mz,m_xlower,m_ylower,m_zlower,m_dx,m_dy,m_dz,m_maux,aux);
#endif

}


Real ClawPatch::ClawPatchIntegrator(FArrayBox& a_phiPatch,
                                    FArrayBox& a_auxPatch,
                                    FArrayBox a_fluxp[],
                                    FArrayBox a_fluxm[],
                                    FArrayBox a_fluxpc[],
                                    FArrayBox a_fluxmc[],
                                    FArrayBox& a_phiCoarse,
                                    FArrayBox& a_auxCoarse,
                                    FArrayBox a_qadd[],
                                    const Box& a_patchBox,
                                    const ProblemDomain& a_domain,
                                    const Real& a_time,
                                    const Real& a_dt,
                                    const Real& a_dx,
                                    const int& a_refRatio,
                                    const int& a_level,
                                    const global_parms& gparms)
{

    // Since m_meqn is hardwired in the constructor for the ClawPatch,
    // we should check it here.
    // CH_assert(a_phiPatch.nComp() == m_meqn);

    /*
    Real dx = a_dx;
    Real dy = dx;
    Real xlower, ylower;
    int mx, my;
#if CH_SPACEDIM == 2
    get_clawvars(a_patchBox,dx,dy,xlower,ylower,mx,my);
#elif CH_SPACEDIM == 3
    Real dz = dx;
    Real zlower;
    int mz;
    get_clawvars(a_patchBox,dx,dy,dz,xlower,ylower,zlower,mx,my,mz);
#endif
    */

  // needed for qad
  // const IntVect& v = a_phiCoarse.size();
  int v[2];
  int mxc = v[0] - 2;  // since there are only one layer of ghost cells
  int myc = v[1] - 2;



#if CH_SPACEDIM == 3
  int mzc = v[2] - 2;
#endif

  Real dt = a_dt;

  // Data for step2 or step3.  This is overwritten by updated values.
  Real* qold = a_phiPatch.dataPtr();
  Real* aux = NULL;
  if (gparms.m_maux > 0)
  {
      aux = a_auxPatch.dataPtr();
  }

  // This is needed for qad
  Real* qold_coarse = a_phiCoarse.dataPtr();

  Real* auxold_coarse = NULL;
  if (gparms.m_maux > 0)
  {
      auxold_coarse = a_auxCoarse.dataPtr();
  }

  // This seems to be different from a_patchBox, because phiBox contains
  // ghost cell indexing, which I need to determine if we are near a physical
  // boundary.
  // Box phiBox = a_phiPatch.box();

  int *intersectsBoundary = new int[SpaceDim*2];
  for (int i = 0; i < 2*SpaceDim; i++)
  {
      intersectsBoundary[i] = 0;
  }

  for(int idir = 0; idir < SpaceDim; idir++)
  {
      if (a_domain.isPeriodic(idir))
      {
          // Skip this boundary;  don't set any boundary conditions.
      }
      else
      {
          /*
          // Don't include this for now, since we don't have routines for adjCellLo and
          // adjCellHi
          Box bdryBox[2];
          bdryBox[0] =  adjCellLo(a_domain,idir,m_meqn);
          bdryBox[1] =  adjCellHi(a_domain,idir,m_meqn);
          for(int edge = 0; edge < 2; edge++)
          {
              if (phiBox.intersects(bdryBox[edge]))
              {
                  // Only keep track of non-periodic physical boundaries
                  intersectsBoundary[2*idir + edge] = 1;
              }
          }
          */
      }
  }


  int maxm = Max(m_mx,m_my);
#if CH_SPACEDIM == 3
  maxm = Max(maxm,mz);
#endif

  // set common block for level
  set_common_levels_(gparms.m_maxlevel,a_level,gparms.m_refratio);


  Real cflgrid;


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

  return cflgrid;
}


void ClawPatch::estimateError(const FArrayBox& a_phiPatch,
                               const Box& a_patchBox,
                               const ProblemDomain& a_domain,
                               const Real& a_time,
                               const Real& a_dt,
                               const Real& a_dx,
                               const int& a_level,
                               const int isBoundary[],
                               const Real& a_refineThresh,
                              FArrayBox& a_error_measure,
                              const global_parms &gparms)
{
    /*
      Real dx = a_dx;
      Real dy = dx;
      Real xlower, ylower;
      int mx, my;
      #if CH_SPACEDIM == 2
      get_clawvars(a_patchBox,dx,dy,xlower,ylower,mx,my);
      #elif CH_SPACEDIM == 3
      Real dz = dx;
      Real zlower;
      int mz;
      get_clawvars(a_patchBox,dx,dy,dz,xlower,ylower,zlower,mx,my,mz);
      #endif
    */

    const Real* q = a_phiPatch.dataPtr();
    Real* em = a_error_measure.dataPtr();
    Real tol = a_refineThresh;

#if CH_SPACEDIM == 2
    estimate_error2_(m_mx,m_my,gparms.m_mbc,gparms.m_meqn,m_xlower, m_ylower,m_dx,m_dy,
                     a_time,a_level,isBoundary, tol, q,em);
#else
    estimate_error3_(mx,my,mz,m_mbc,m_meqn,xlower, ylower, zlower,
                     dx,dy,dz,a_time, a_level, isBoundary, tol, q, em);
#endif

}

/*
#if CH_SPACEDIM == 2
void ClawPatch::get_clawvars(const Box& a_box, const Real& dx, const Real& dy,
                             Real& xlower, Real& ylower, int& mx, int& my)
{
    mx = a_box.bigEnd(0) - a_box.smallEnd(0) + 1;
    my = a_box.bigEnd(1) - a_box.smallEnd(1) + 1;

    // Lower edge of computational domain.
    xlower = m_xlower + dx*a_box.smallEnd(0);
    ylower = m_ylower + dy*a_box.smallEnd(1);
}
#elif CH_SPACEDIM == 3

void ClawPatch::get_clawvars(const Box& a_box, const Real& dx, const Real& dy, const Real& dz,
                             Real& xlower, Real& ylower, Real& zlower,
                             int& mx, int& my, int& mz)
{
    mx = a_box.bigEnd(0) - a_box.smallEnd(0) + 1;
    my = a_box.bigEnd(1) - a_box.smallEnd(1) + 1;
    mz = a_box.bigEnd(2) - a_box.smallEnd(2) + 1;

    // Lower edge of computational domain.
    xlower = m_xlower + dx*a_box.smallEnd(0);
    ylower = m_ylower + dy*a_box.smallEnd(1);
    zlower = m_zlower + dz*a_box.smallEnd(2);
}

#endif
*/
