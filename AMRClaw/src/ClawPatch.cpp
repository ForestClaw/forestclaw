#include "ClawPatch.H"
#include "amr_utils.H"

// This constructors includes all of parameters that are patch independent.
// All of this could also be in some sort of "set_params" function...
ClawPatch::ClawPatch()
{
    // This is dimensioned here, rather than in the more logical place,
    // set_method, because other routines, like set_mcapa and set_src
    // will attempt to access this array, possibly before set_method has
    // been called.
    m_method = new int[7];
    for (int i = 0; i < 7; i++)
    {
        m_method[i] = 0;
    }

    // Other arrays, like m_mthbc and m_mthlim, are
    // dimensioned in their 'set' routines.
    m_mthbc = NULL;  // NULL pointers
    m_mthlim = NULL;
    m_refratios = NULL;
    m_mwavesIsSet = false; // needed for proper dimensioning of mthlim
    m_maxlevelIsSet = false;
    m_isDefined = false;
}

ClawPatch::~ClawPatch()
{
    if( m_mthlim      != NULL ) { delete [] m_mthlim     ; m_mthlim      = NULL; }
    if( m_method      != NULL ) { delete [] m_method     ; m_method      = NULL; }
    if( m_mthbc       != NULL ) { delete [] m_mthbc      ; m_mthbc       = NULL; }
    if( m_refratios   != NULL ) { delete [] m_refratios  ; m_refratios   = NULL; }
    if( m_auxtype_int != NULL ) { delete [] m_auxtype_int; m_auxtype_int = NULL; }
}

void ClawPatch::define(const ClawPatch& a_clawPatch)
{
    set_mx(a_clawPatch.m_mx);
    set_my(a_clawPatch.m_my);
    set_method(a_clawPatch.m_method);
    set_mbc(a_clawPatch.m_mbc);
    set_meqn(a_clawPatch.m_meqn);
    set_maux(a_clawPatch.m_maux);
    set_mcapa(a_clawPatch.m_mcapa);
    set_mwaves(a_clawPatch.m_mwaves);
    set_mthlim(a_clawPatch.m_mthlim);
    set_mthbc(a_clawPatch.m_mthbc);
    set_stateNames(a_clawPatch.m_stateNames);
    set_auxtype(a_clawPatch.m_auxtype);
    set_xlower(a_clawPatch.m_xlower);
    set_ylower(a_clawPatch.m_ylower);
    set_xupper(a_clawPatch.m_xupper);
    set_yupper(a_clawPatch.m_yupper);
    set_initial_dt(a_clawPatch.m_initial_dt);
    set_maxlevel(a_clawPatch.m_maxlevel);
    set_refratios(a_clawPatch.m_refratios);
    set_max_cfl(a_clawPatch.m_max_cfl);

#if CH_SPACEDIM == 3
    set_mz(a_clawPatch.m_mz);
    set_zlower(a_clawPatch.m_zlower);
    set_zupper(a_clawPatch.m_zupper);
#endif

    m_isDefined = true;
}

bool ClawPatch::isDefined()
{
    return m_isDefined;
}

void ClawPatch::get_inputParams()
{
    global_parms parms;

    // read data using Fortran file
    inputparms_(parms.m_initial_dt,
                parms.m_tfinal,
                parms.m_max_cfl,
                parms.m_nout,
                parms.m_src_term,
                parms.m_mcapa,
                parms.m_maux,
                parms.m_meqn,
                parms.m_mwaves,
                parms.m_maxmwaves,
                parms.m_mthlim,
                parms.m_mbc,
                parms.m_mthbc,
                parms.m_order);

    // cp->set_mx(domain->mx_leaf);
    // cp->set_my(domain->my_leaf);

    set_src(parms.m_src_term);
    set_mcapa(parms.m_mcapa);
    set_maux(parms.m_maux);
    set_mbc(parms.m_mbc);
    set_meqn(parms.m_meqn);
    set_mwaves(parms.m_mwaves);
    set_initial_dt(parms.m_initial_dt);

    set_max_cfl(parms.m_max_cfl);


    // Do this nonsense because I used std::vector<int> in ClawPatch.
    std::vector<int> order(2,0);
    for(int m = 0; m < 2; m++)
    {
        order[m] = parms.m_order[m];
    }
    set_order(order);

    std::vector<int> mthlim(parms.m_mwaves,0);
    for (int m = 0; m < parms.m_mwaves; m++)
    {
        mthlim[m] = parms.m_mthlim[m];
    }
    set_mthlim(mthlim);

    std::vector<int> mthbc(4,0);
    for (int m = 0; m < 4; m++)
    {
        mthbc[m] = parms.m_mthbc[m];
    }
    set_mthbc(mthbc);

}

void ClawPatch::print_inputParams()
{
    // Killed everything that used "pout()"

  pout() << endl;
  pout() << "CLAWPACK PARAMETERS : " << endl;
  pout() << "Initial dt " << m_initial_dt << endl;
  pout() << "maximum cfl " << m_max_cfl << endl;
  // pout() << "method(1) (fixed time step)       = " << fixedDt << endl;
  pout() << "method(2:3) (order of integration) = " << m_method[1] << " " << m_method[2] << endl;
  pout() << "method(5) (source term splitting) = " << m_method[4] << endl;
  pout() << "method(6) (mcapa) = " << m_method[5] << endl;
  pout() << "method(7) (maux) = " << m_method[6] << endl;
  pout() << endl;

  pout() << "meqn (number of equations) = " << m_meqn << endl;
  pout() << "maux (number of auxiliary variables) = " << m_maux << endl;
  pout() << "mcapa (location of capacity function) = " << m_mcapa << endl;
  pout() << "mwaves (number of waves) = " << m_mwaves << endl;

  pout() << "mthlim(mwaves) (limiters) = ";
  for(int i = 0; i < m_mwaves; i++)
    {
      pout() << m_mthlim[i] << " ";
    }
  pout() << endl << endl;

  pout() << "mbc (number of ghost cells) = " << m_mbc << endl;

  pout() << "Domain values " << endl;
  pout() << "xlower, xupper = " << m_xlower << " " << m_xupper << endl;
  pout() << "ylower, yupper = " << m_ylower << " " << m_yupper << endl;

#if CH_SPACEDIM == 3
  pout() << "zlower, zupper = " << m_zlower << " " << m_zupper << endl;
#endif

  /*
  pout() << "Auxiliary array type : " << endl;
  for (int i = 0; i < m_maux; i++)
    {
      pout() << "  " << m_auxtype[i] << endl;
    }
  pout() << endl;
  */


  pout() << "mthbc(2*dim) (boundary conditions) = ";
  for(int i = 0; i < 2*SpaceDim; i++)
    {
      pout() << m_mthbc[i] << " ";
    }
  pout() << endl << endl;


}

// clawpack parameters
void ClawPatch::set_mx(int a_mx)
{
    m_mx = a_mx;
}


void ClawPatch::set_my(int a_my)
{
    m_my = a_my;
}

#if (CH_SPACEDIM == 3)
void ClawPatch::set_mz(int a_mz)
{
    m_mz = a_mz;
}
#endif

void ClawPatch::set_initial_dt(Real a_initial_dt)
{
    m_initial_dt = a_initial_dt;
}

void ClawPatch::set_mbc(int a_mbc)
{
  m_mbc = a_mbc;
}

void ClawPatch::set_meqn(int a_meqn)
{
  m_meqn = a_meqn;
}

void ClawPatch::set_mwaves(int a_mwaves)
{
  m_mwaves = a_mwaves;
  m_mwavesIsSet = true;
}

void ClawPatch::set_mthlim(std::vector<int> a_mthlim)
{
    if (m_mthlim != NULL)
    {
        delete m_mthlim;
    }
    m_mthlim = new int[a_mthlim.size()];
    for (int i = 0; i < a_mthlim.size(); i++)
    {
        m_mthlim[i] = a_mthlim[i];
    }
}

void ClawPatch::set_mthlim(int* a_mthlim)
{
    // CH_assert(m_mwavesIsSet);
    if (m_mthlim != NULL)
    {
        delete m_mthlim;
    }
    m_mthlim = new int[get_mwaves()];
    for (int i = 0; i < m_mwaves; i++)
    {
        m_mthlim[i] = a_mthlim[i];
    }
}


// method[0] = dt fixed (not used in ChomboClaw)
// method[1] = order of integration of for normal direction
// method[2] = order of integration of transverse direction(s)
// method[3] = verbosity level (not used in ChomboClaw)
// method[4] = source term splitting
// method[5] = mcapa
// method[6] = maux

void ClawPatch::set_order(std::vector<int> a_order)
{
    m_method[1] = a_order[0];
    if (SpaceDim == 2)
    {
        m_method[2] = a_order[1];
    }
    else
    {
        m_method[2] = 10*a_order[1] + a_order[2];
    }
}

void ClawPatch::set_src(int a_src)
{
  m_method[4] = a_src;
}

void ClawPatch::set_mcapa(int a_mcapa)
{
    m_mcapa = a_mcapa;
    m_method[5] = a_mcapa;
}

void ClawPatch::set_maux(int a_maux)
{
    m_maux = a_maux;
    m_method[6] = a_maux;
}

void ClawPatch::set_mthbc(std::vector<int> a_mthbc)
{
    if (m_mthbc != NULL)
    {
        delete m_mthbc;
    }
    m_mthbc = new int[2*SpaceDim];
    for(int i = 0; i < 2*SpaceDim; i++)
    {
        m_mthbc[i] = a_mthbc[i];
    }
}

void ClawPatch::set_mthbc(int* a_mthbc)
{
    if (m_mthbc != NULL)
    {
        delete m_mthbc;
    }
    m_mthbc = new int[2*SpaceDim];
    for (int i = 0; i < 2*SpaceDim; i++)
    {
        m_mthbc[i] = a_mthbc[i];
    }
}


void ClawPatch::set_stateNames(std::vector<std::string> a_stateNames)
{
  m_stateNames = a_stateNames;
}

void ClawPatch::set_auxtype(std::vector<std::string> a_auxtype)
{
    m_auxtype = a_auxtype;
    int maux = a_auxtype.size();

    // This is used to communicate with fortran routines, which don't know about
    // std::vector<string>...
    m_auxtype_int = new int[maux];
    for(int m = 0; m < maux; m++)
    {
        if (a_auxtype[m] == "xleft")
        {
            m_auxtype_int[m] = 1;
        }
        else if (a_auxtype[m] == "yleft")
        {
            m_auxtype_int[m] = 2;
        }
        else if (a_auxtype[m] == "zleft")
        {
            m_auxtype_int[m] = 3;
        }
        else if (a_auxtype[m] == "center")
        {
            m_auxtype_int[m] = 4;
        }
        else if (a_auxtype[m] == "capacity")
        {
            m_auxtype_int[m] = 5;
        }
    }
}

// Private, only to be called by methods in this class.
void ClawPatch::set_method(std::vector<int> a_method)
{
    for (int i = 0; i < 7; i++)
    {
        m_method[i] = a_method[i];
    }
}

void ClawPatch::set_method(int *a_method)
{
    for (int i = 0; i < 7; i++)
    {
        m_method[i] = a_method[i];
    }
}


void ClawPatch::set_xlower(Real a_xlower)
{
    m_xlower = a_xlower;
}

void ClawPatch::set_ylower(Real a_ylower)
{
    m_ylower = a_ylower;
}

#if CH_SPACEDIM == 3
void ClawPatch::set_zlower(Real a_zlower)
{
    m_zlower = a_zlower;
}
#endif


void ClawPatch::set_xupper(Real a_xupper)
{
    m_xupper = a_xupper;
}

void ClawPatch::set_yupper(Real a_yupper)
{
    m_yupper = a_yupper;
}

#if CH_SPACEDIM == 3
void ClawPatch::set_zupper(Real a_zupper)
{
    m_zupper = a_zupper;
}
#endif



void ClawPatch::set_refratios(std::vector<int> a_refratios)
{
    if (m_refratios != NULL)
    {
        delete m_refratios;
    }
    m_refratios = new int[a_refratios.size()];
    for (int i = 0; i < a_refratios.size(); i++)
    {
        m_refratios[i] = a_refratios[i];
    }
}

void ClawPatch::set_refratios(int* a_refratios)
{
    // CH_assert(m_maxlevelIsSet);
    if (m_refratios != NULL)
    {
        delete m_refratios;
    }
    m_refratios = new int[m_maxlevel];
    for (int i = 0; i < m_maxlevel; i++)
    {
        m_refratios[i] = a_refratios[i];
    }
}


void ClawPatch::set_maxlevel(const int& a_maxlevel)
{
    m_maxlevel = a_maxlevel;
    m_maxlevelIsSet = true;
}

int* ClawPatch::get_refratios() const
{
    return m_refratios;
}

int ClawPatch::get_maxlevel() const
{
    return m_maxlevel;
}


Real ClawPatch::get_max_cfl() const
{
    return m_max_cfl;
}

void ClawPatch::set_max_cfl(const Real& a_max_cfl)
{
    m_max_cfl = a_max_cfl;
}


int ClawPatch::get_mx() const
{
    return m_mx;
}

int ClawPatch::get_my() const
{
    return m_my;
}

Real ClawPatch::get_initial_dt() const
{
    return m_initial_dt;
}

#if CH_SPACEDIM == 3
int ClawPatch::get_mz() const
{
    return m_mz;
}
#endif

int ClawPatch::get_meqn() const
{
  return m_meqn;
}

int ClawPatch::get_mbc() const
{
  return m_mbc;
}

int ClawPatch::get_maux() const
{
    return m_maux;
}

int ClawPatch::get_mcapa() const
{
    return m_mcapa;
}

int ClawPatch::get_mwaves() const
{
    CH_assert(m_mwavesIsSet);
    return m_mwaves;
}

int* ClawPatch::get_mthbc() const
{
    return m_mthbc;
}

std::vector<std::string> ClawPatch::get_stateNames() const
{
  return m_stateNames;
}

std::vector<std::string> ClawPatch::get_auxtype() const
{
  return m_auxtype;
}

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


//void ClawPatch::initialize(LevelData<FArrayBox>& a_phi, LevelData<FArrayBox>& a_aux, Real a_dx)
void ClawPatch::initialize(FArrayBox& a_phi, FArrayBox& a_aux, const Box& a_box, Real a_dx)
{
    CH_assert(m_isDefined);

    int mx, my;
    Real xlower, ylower;
    Real dx = a_dx;
    Real dy = a_dx;
#if CH_SPACEDIM == 2
    get_clawvars(a_box, dx, dy, xlower, ylower, mx, my);
#elif CH_SPACEDIM == 3
    Real dz = a_dx;
    int mz;
    Real zlower;
    get_clawvars(a_box, dx,dy, dz, xlower, ylower, zlower, mx, my, mz);
#endif

    Real* q = a_phi.dataPtr();
    Real* aux = a_aux.dataPtr();

    //  c     =====================================================
    //         subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
    //       &                   dx,dy,q,maux,aux)
    //  c     =====================================================

    // first two arguments are maxmx/maxmy
    // is it ok to send mx/my for that?  (ndk)
#if CH_SPACEDIM == 2
    qinit_(mx,my,m_meqn,m_mbc,mx,my,xlower,ylower,dx,dy,q,m_maux,aux);
#elif CH_SPACEDIM == 3
    qinit_(mx,my,mz,m_meqn,m_mbc,mx,my,mz,xlower,ylower,zlower,
           dx,dy,dz,q,m_maux,aux);
#endif
}

void ClawPatch::setAuxArray(FArrayBox& a_aux, const Box& a_box, Real a_dx,
                            const int& a_level)
{
    CH_assert(m_isDefined);

    int mx, my;
    Real xlower, ylower;
    Real dx = a_dx;
    Real dy = a_dx;

#if CH_SPACEDIM == 2
    get_clawvars(a_box, dx, dy, xlower, ylower, mx, my);
#elif CH_SPACEDIM == 3
    Real dz = a_dx;
    int mz;
    Real zlower;
    get_clawvars(a_box, dx, dy, dz, xlower, ylower, zlower, mx, my, mz);
#endif

    // set common block for level
    set_common_levels_(m_maxlevel,a_level,m_refratios);

    Real* aux = a_aux.dataPtr();

    // Default setaux.f routine.  (2d version)
    // c     ============================================
    //       subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
    //      &                  maux,aux)
    // c     ============================================

    // first two arguments are maxmx/maxmy
    // is it ok to send mx/my for that?  (ndk)
    // yes - maxmx, maxmy only used for static dimensioning of arrays in Clawpack (DAC).
#if CH_SPACEDIM == 2
    setaux_(mx,my,m_mbc,mx,my,xlower,ylower,dx,dy,m_maux,aux);
#elif CH_SPACEDIM == 3
    setaux_(mx,my,mz,m_mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,m_maux,aux);
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
                                    const int& a_level)
{

  // Since m_meqn is hardwired in the constructor for the ClawPatch,
  // we should check it here.
    // CH_assert(a_phiPatch.nComp() == m_meqn);

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
  if (m_maux > 0)
  {
      aux = a_auxPatch.dataPtr();
  }

  // This is needed for qad
  Real* qold_coarse = a_phiCoarse.dataPtr();

  Real* auxold_coarse = NULL;
  if (m_maux > 0)
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


  int maxm = Max(mx,my);
#if CH_SPACEDIM == 3
  maxm = Max(maxm,mz);
#endif

  // set common block for level
  set_common_levels_(m_maxlevel,a_level,m_refratios);


  Real cflgrid;


#if CH_SPACEDIM==2
  int mwork = (maxm+2*m_mbc)*(12*m_meqn + (m_meqn+1)*m_mwaves + 3*m_maux + 2);
  Real* work = new Real[mwork];

  Real* fp = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)];
  Real* fm = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)];
  Real* gp = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)];
  Real* gm = new Real[m_meqn*(mx+2*m_mbc)*(my+2*m_mbc)];

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

  clawpatch2_(maxm, m_meqn, m_maux, m_mbc, m_method, m_mthlim,
              m_mcapa, m_mwaves, mx, my, qold, aux,
              dx, dy, dt, cflgrid, work, mwork, qold_coarse, auxold_coarse,
              qadd_x, qadd_y, m_auxtype_int, xlower, ylower,
              intersectsBoundary,a_level,
              m_mthbc, a_time, mxc, myc, fp, fm, gp, gm,
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
                               FArrayBox& a_error_measure)
{
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

  const Real* q = a_phiPatch.dataPtr();
  Real* em = a_error_measure.dataPtr();
  Real tol = a_refineThresh;

#if CH_SPACEDIM == 2
  estimate_error2_(mx,my,m_mbc,m_meqn,xlower, ylower,dx,dy,a_time,a_level,
                  isBoundary, tol, q,em);
#else
  estimate_error3_(mx,my,mz,m_mbc,m_meqn,xlower, ylower, zlower,
                 dx,dy,dz,a_time, a_level, isBoundary, tol, q, em);
#endif

}

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
