#include "phasefield_patch_operator.h"

using namespace ThunderEgg;

double phasefield::lambda{999};

void phasefield::setLambda(double lambda){
    phasefield::lambda = lambda;
}
double phasefield::getLambda(){
    return lambda;
}

phasefield::phasefield(fclaw2d_global_t *glob,
                       const Vector<2>& phi_n,
                       const Domain<2>& domain,
                       const GhostFiller<2>& ghost_filler) 
    : phasefield(fc2d_thunderegg_get_options(glob),phasefield_get_options(glob),phi_n,domain,ghost_filler) 
{
    //this just calls the other constructor
}
phasefield::phasefield(const fc2d_thunderegg_options *mg_opt,const phasefield_options* phase_opt,
                       const Vector<2>& phi_n_in,
                       const Domain<2>& domain,
                       const GhostFiller<2>& ghost_filler) 
    : PatchOperator<2>(domain,ghost_filler),
      phi_n(phi_n_in),
      phase_opt(phase_opt)

{
    /* User should call 'fc2d_thunderegg_phasefield_set_lambda' before calling elliptic solve */
    FCLAW_ASSERT(phasefield::lambda <= 0);

    /* Get scale needed to apply homogeneous boundary conditions. 
       For Dirichlet (bctype=1) :   scalar is -1 
       For Neumann (bctype=2)   :   scalar is 1 */
    for(int m = 0; m < 4; m++)
    {
        s[m] = 2*mg_opt->boundary_conditions[m] - 3;
    }
}


phasefield* phasefield::clone() const{
    return new phasefield(*this);
}
void phasefield::applySinglePatchWithInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                                const PatchView<const double,2>& us,
                                                                const PatchView<double,2>& fs) const 
{
    int mfields = us.getEnd()[2]+1;
    int mx =pinfo.ns[0]; 
    int my =pinfo.ns[1];

    if (pinfo.hasNbr(Side<2>::west()))
    {
        auto ghosts = us.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++)
            for(int j = 0; j < my; j++)
                ghosts(j,m) = -us(0,j,m);
    }

    if (pinfo.hasNbr(Side<2>::east()))
    {
        auto ghosts = us.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++)
            for(int j = 0; j < my; j++)
                ghosts(j,m) = -us(mx-1,j,m);
    }

    if (pinfo.hasNbr(Side<2>::south()))
    {
        auto ghosts = us.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++)
            for(int i = 0; i < mx; i++)
                ghosts(i,m) = -us(i,0,m);
    }

    if (pinfo.hasNbr(Side<2>::north()))
    {
        auto ghosts = us.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++)
            for(int i = 0; i < mx; i++)
                ghosts(i,m) = -us(i,my-1,m);
    }
    applySinglePatch(pinfo,us,fs);
}
void phasefield::applySinglePatch(const PatchInfo<2>& pinfo, 
                                  const PatchView<const double,2>& us,
                                  const PatchView<double,2>& fs) const 
{
    //const cast since u ghost values have to be modified
    //ThunderEgg doesn't care if ghost values are modified, just don't modify the interior values.

    int mfields = us.getEnd()[2]+1;
    int mx =pinfo.ns[0]; 
    int my =pinfo.ns[1];

#if 0    
    int mbc =pinfo.num_ghost_cells;
    double xlower =pinfo.starts[0];
    double ylower =pinfo.starts[1];
#endif    

    /* Apply boundary conditions */
    if (!pinfo.hasNbr(Side<2>::west()))
    {
        /* Physical boundary */
        auto ghosts = us.getGhostSliceOn(Side<2>::west(),{0});
        for(int m = 0; m < mfields; m++)
            for(int j = 0; j < my; j++)
                ghosts(j,m) = s[0]*us(0,j,m);
    }

    if (!pinfo.hasNbr(Side<2>::east()))
    {
        /* Physical boundary */
        auto ghosts = us.getGhostSliceOn(Side<2>::east(),{0});
        for(int m = 0; m < mfields; m++)
            for(int j = 0; j < my; j++)
                ghosts(j,m) = s[1]*us(mx-1,j,m);            
    } 

    if (!pinfo.hasNbr(Side<2>::south()))
    {
        /* Physical boundary */
        auto ghosts = us.getGhostSliceOn(Side<2>::south(),{0});
        for(int m = 0; m < mfields; m++)
            for(int i = 0; i < mx; i++)
                ghosts(i,m) = s[2]*us(i,0,m);
    }

    if (!pinfo.hasNbr(Side<2>::north()))
    {
        /* Physical boundary */
        auto ghosts = us.getGhostSliceOn(Side<2>::north(),{0});
        for(int m = 0; m < mfields; m++)
            for(int i = 0; i < mx; i++)
                ghosts(i,m) = s[3]*us(i,my-1,m);
    }

#if 1
    /* Five-point Laplacian - not anisotropic yet */
    ComponentView<const double,2> u = us.getComponentView(0);
    ComponentView<double,2> Au = fs.getComponentView(0);

    ComponentView<const double,2> phi = us.getComponentView(1);
    ComponentView<double,2> Aphi = fs.getComponentView(1);

    double dx =pinfo.spacings[0];
    double dy =pinfo.spacings[1];
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    /* Check already done at construction, but this is a double check */
    FCLAW_ASSERT(lambda <= 0);

    double xi = phase_opt->xi;
    double m = phase_opt->m;
    double S = phase_opt->S;
    double alpha = phase_opt->alpha;
    double beta = xi*xi/m;

    double c1 = 30/S;
    double c2 = 30*alpha*S*xi;

    /*  For the isotropic case 

              T(\theta) = xi^2; 

        For the anisotropic case

              T(\theta) = \eta(\theta)[\eta(\theta), -\eta'(\theta);
                                       \eta'(\theta), \eta(\theta)]

        In the anisotropic case, we have to compute T inside the loop.
    */
    double  T = xi*xi;

    /* Get local view into phi_n : component 1 */
    ComponentView<const double,2> pn = phi_n.getComponentView(1, pinfo.local_index);

    for(int j = 0; j < my; j++)
        for(int i = 0; i < mx; i++)
        {
            double uij   = u(i,j);
            double lap_u = (u(i+1,j) - 2*uij + u(i-1,j))/dx2 + 
                           (u(i,j+1) - 2*uij + u(i,j-1))/dy2;

            double phi_ij = phi(i,j);
            double lap_phi = (phi(i+1,j) - 2*phi_ij + phi(i-1,j))/dx2 + 
                             (phi(i,j+1) - 2*phi_ij + phi(i,j-1))/dy2;

            /* Use phi_n from previous time step */
            double pn_ij = pn(i,j);
            double g0 = pn_ij*(1-pn_ij);
            double g = g0*g0;
            double S1 = c1*g;
            double S2 = c2*g;

            Au(i,j)   = lap_u + lambda*(uij + S1*phi_ij);
            //f_phi(i,j) = T*lap_phi + lambda*(1/lambda*S2*uij + beta*phi_ij);
            Aphi(i,j) = T*lap_phi + S2*uij + lambda*beta*phi_ij;
        }
    
#else

        /* Five-point Laplacian : Slightly slower than then above.*/
        for(int j = 0; j < my; j++)
            for(int i = 0; i < mx; i++)
            {
                double uij = u(i,j);
                double flux[4];
                flux[0] = (uij - u(i-1,j));
                flux[1] = (u(i+1,j) - uij);
                flux[2] = (uij - u(i,j-1));
                flux[3] = (u(i,j+1) - uij);;
                f(i,j) = (flux[1]-flux[0])/dx2 + (flux[3] - flux[2])/dy2;
            }
#endif    
    
}


void phasefield::modifyRHSForInternalBoundaryConditions(const PatchInfo<2>& pinfo, 
                                                        const PatchView<const double,2>& us, 
                                                        const PatchView<double,2>& Aus) const 
{
    int mfields = us.getEnd()[2]+1;
    int mx =pinfo.ns[0]; 
    int my =pinfo.ns[1];


#if 0

    double dx =pinfo.spacings[0];
    double dx2 = dx*dx;

    double dy =pinfo.spacings[1];
    double dy2 = dy*dy;

    for(int m = 0; m < mfields; m++)
    {
        const LocalData<2>& u = us[m];
        LocalData<2>& Au = Aus[m];
        for(int j = 0; j < my; j++)
        {
            /* bool hasNbr(Side<D> s) */
            if pinfo.hasNbr(Side<2>::west()))
                Au(0,j) += -(u(-1,j)+u(0,j))/dx2;

            if pinfo.hasNbr(Side<2>::east()))
                Au(mx-1,j) += -(u(mx-1,j)+u(mx,j))/dx2;
        }

        for(int i = 0; i < mx; i++)
        {
            if pinfo.hasNbr(Side<2>::south()))
                Au(i,0) += -(u(i,-1)+u(i,0))/dy2;

            if pinfo.hasNbr(Side<2>::north()))
                Au(i,my-1) += -(u(i,my-1)+u(i,my))/dy2;
        }
    }
#else       
    PatchArray<2> new_us(pinfo.ns,mfields,1);
    PatchArray<2> new_Aus(pinfo.ns,mfields,1);
    if(pinfo.hasNbr(Side<2>::west())){
        for(int field=0; field<mfields; field++){
            for(int j=0; j < my; j++){
                new_us(-1,j,field)=us(-1,j,field)+us(0,j,field);
            }
        }
    }
    if(pinfo.hasNbr(Side<2>::east())){
        for(int field=0; field<mfields; field++){
            for(int j=0; j < my; j++){
                new_us(my,j,field)=us(my,j,field)+us(my-1,j,field);
            }
        }
    }
    if(pinfo.hasNbr(Side<2>::south())){
        for(int field=0; field<mfields; field++){
            for(int i=0; i < mx; i++){
                new_us(i,-1,field)=us(i,-1,field)+us(i,0,field);
            }
        }
    }
    if(pinfo.hasNbr(Side<2>::north())){
        for(int field=0; field<mfields; field++){
            for(int i=0; i < mx; i++){
                new_us(i,mx,field)=us(i,mx,field)+us(i,mx-1,field);
            }
        }
    }
    applySinglePatch(pinfo,new_us.getView(),new_Aus.getView());
    //modify rhs
    for(int field=0; field<mfields; field++){
        for(int j=0; j < my; j++){
            for(int i=0; i < mx; i++){
                Aus(i,j,field)-=new_Aus(i,j,field);
            }
        }
    }
#endif    
}

const int *phasefield::getS() const{
    return s;
}