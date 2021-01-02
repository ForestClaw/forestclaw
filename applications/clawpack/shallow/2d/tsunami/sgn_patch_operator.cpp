#include "sgn_patch_operator.h"
#include <ThunderEgg/ValVector.h>

using namespace ThunderEgg;

sgn::sgn(fclaw2d_global_t *glob,
               std::shared_ptr<const Vector<2>> q_n,
               std::shared_ptr<const Domain<2>> domain,
               std::shared_ptr<const GhostFiller<2>> ghost_filler) 
    : sgn(fc2d_thunderegg_get_options(glob),tsunami_get_options(glob),q_n,domain,ghost_filler) 
{
    //this just calls the other constructor
}
sgn::sgn(const fc2d_thunderegg_options *mg_opt,const user_options* sgn_opt,
               std::shared_ptr<const Vector<2>> q_n_in,
               std::shared_ptr<const Domain<2>> domain,
               std::shared_ptr<const GhostFiller<2>> ghost_filler) 
    : PatchOperator<2>(domain,ghost_filler),
      q_n(q_n_in),
      sgn_opt(sgn_opt)

{
    /* Get scale needed to apply homogeneous boundary conditions. 
       For Dirichlet (bctype=1) :   scalar is -1 
       For Neumann (bctype=2)   :   scalar is 1 */
    for(int m = 0; m < 4; m++)
    {
        s[m] = 2*mg_opt->boundary_conditions[m] - 3;
    }
}


void sgn::applySinglePatch(std::shared_ptr<const PatchInfo<2>> pinfo, 
                                 const std::vector<LocalData<2>>& us,
                                 std::vector<LocalData<2>>& fs,
                                 bool interior_dirichlet) const 
{
    //const cast since u ghost values have to be modified
    //ThunderEgg doesn't care if ghost values are modified, just don't modify the interior values.

    int mfields = us.size();
    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];

#if 0    
    int mbc = pinfo->num_ghost_cells;
    double xlower = pinfo->starts[0];
    double ylower = pinfo->starts[1];
#endif    

    /* Apply boundary conditions */
    for(int m = 0; m < mfields; m++)
    {
        LocalData<2>& u = const_cast<LocalData<2>&>(us[m]);

        if (!pinfo->hasNbr(Side<2>::west()))
        {
            /* Physical boundary */
            auto ghosts = u.getGhostSliceOnSide(Side<2>::west(),1);
            for(int j = 0; j < my; j++)
                ghosts[{j}] = s[0]*u[{0,j}];
        }
        else if (interior_dirichlet)
        {
            auto ghosts = u.getGhostSliceOnSide(Side<2>::west(),1);
            for(int j = 0; j < my; j++)
                ghosts[{j}] = -u[{0,j}];
        }

        if (!pinfo->hasNbr(Side<2>::east()))
        {
            /* Physical boundary */
            auto ghosts = u.getGhostSliceOnSide(Side<2>::east(),1);
            for(int j = 0; j < my; j++)
                ghosts[{j}] = s[1]*u[{mx-1,j}];            
        } 
        else if (interior_dirichlet)
        {
            auto ghosts = u.getGhostSliceOnSide(Side<2>::east(),1);
            for(int j = 0; j < my; j++)
                ghosts[{j}] = -u[{mx-1,j}];
        }

        if (!pinfo->hasNbr(Side<2>::south()))
        {
            /* Physical boundary */
            auto ghosts = u.getGhostSliceOnSide(Side<2>::south(),1);
            for(int i = 0; i < mx; i++)
                ghosts[{i}] = s[2]*u[{i,0}];
        }
        else if (interior_dirichlet)
        {
            auto ghosts = u.getGhostSliceOnSide(Side<2>::south(),1);
            for(int i = 0; i < mx; i++)
                ghosts[{i}] = -u[{i,0}];
        }

        if (!pinfo->hasNbr(Side<2>::north()))
        {
            /* Physical boundary */
            auto ghosts = u.getGhostSliceOnSide(Side<2>::north(),1);
            for(int i = 0; i < mx; i++)
                ghosts[{i}] = s[3]*u[{i,my-1}];
        }
        else if (interior_dirichlet)
        {
            auto ghosts = u.getGhostSliceOnSide(Side<2>::north(),1);
            for(int i = 0; i < mx; i++)
                ghosts[{i}] = -u[{i,my-1}];
        }
    }

    /* Five-point Laplacian - not anisotropic yet */
    LocalData<2>& Dx = const_cast<LocalData<2>&>(us[0]);
    LocalData<2>& Fx = fs[0];

    LocalData<2>& Dy = const_cast<LocalData<2>&>(us[1]);
    LocalData<2>& Fy = fs[1];

    double dx = pinfo->spacings[0];
    double dy = pinfo->spacings[1];
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    /* We need to discretize this (from Basilisk.fr)

    res.x[] = b.x[] -
      (-alpha_d/3.*(hr3*D.x[1] + hl3*D.x[-1] - 
            (hr3 + hl3)*D.x[])/sq(Delta) +
       hc*(alpha_d*(dxeta*dxzb + hc/2.*d2x(zb)) + 1.)*D.x[] +
       alpha_d*hc*((hc/2.*d2xy(zb) + dxeta*dy(zb))*D.y[] + 
               hc/2.*dy(zb)*dx(D.y) - sq(hc)/3.*d2xy(D.y)
               - hc*dy(D.y)*(dxh + dxzb/2.)));

    Simplfied : 
       -alpha_d/3.*(  hr3*D_{i+1,j} + hl3*D_{i-1,j} - (hr3 + hl3)*D{i,j}  )/sq(Delta)
                     + hc*D_{i,j}
       + alpha_d*hc*(  -sq(hc)/3.*d2xy(D.y) - hc*dy(D.y)*(dxh)  );
    */

    double alpha = sgn_opt->alpha;


    /* Get local view into q_n (component 0) */
    LocalData<2> h = q_n->getLocalData(0,pinfo->local_index);


    for(int j = 0; j < my; j++)
    {
        for(int i = 0; i < mx; i++)
        {

            double hc = h[{i,j}];

            double hl = h[{i-1,j}];
            double hr = h[{i+1,j}];
            double hu = h[{i,j+1}];
            double hd = h[{i,j-1}];

            /* Compute h^3 at edges */
            double hl3 = pow((hl + hc)/2.0,3);
            double hr3 = pow((hr + hc)/2.0,3);
            double hu3 = pow((hu + hc)/2.0,3);
            double hd3 = pow((hd + hc)/2.0,3);

            double dxh = (hr - hl)/(2*dx);
            double dyh = (hu - hd)/(2*dy);

            double hc2 = hc*hc;

            if (hc <= 0)
            {
                fclaw_global_essentialf("sgn_patch_operator : h(i,j) = 0\n");
                exit(0);
            }

            /* Operator without bathymetry */
            double dxh3Dx = (hl3*Dx[{i-1,j}] - (hl3+hr3)*Dx[{i,j}] + hr3*Dx[{i+1,j}])/dx2;
            double dyh3Dy = (hd3*Dy[{i,j-1}] - (hu3+hd3)*Dy[{i,j}] + hu3*Dy[{i,j+1}])/dy2;

            double dxyDx = ((Dx[{i+1,j+1}] - Dx[{i+1,j-1}]) - 
                            (Dx[{i-1,j+1}] - Dx[{i-1,j-1}]))/(4*dx*dy);

            double dxyDy = ((Dy[{i+1,j+1}] - Dy[{i+1,j-1}]) - 
                            (Dy[{i-1,j+1}] - Dy[{i-1,j-1}]))/(4*dx*dy);

            double dxDx = (Dx[{i+1,j}] - Dx[{i-1,j}])/(2*dx);
            double dyDy = (Dy[{i,j+1}] - Dy[{i,j-1}])/(2*dy);            

            Fx[{i,j}]   = -alpha/3.0*dxh3Dx + hc*Dx[{i,j}] - 
                           alpha*hc*(hc*dxh*dyDy + hc2*dxyDy/3.0);
            Fy[{i,j}]   = -alpha/3.0*dyh3Dy + hc*Dy[{i,j}] - 
                           alpha*hc*(hc*dyh*dxDx + hc2*dxyDx/3.0);

            printf("%5d %5d %12.4e %12.4e %12.4e\n",i,j,h[{i,j}],Fx[{i,j}],Fy[{i,j}]);
        }
    }
    printf("Done with patch iteration\n\n");
    
}


void sgn::addGhostToRHS(std::shared_ptr<const PatchInfo<2>> pinfo, 
                              const std::vector<LocalData<2>>& us, 
                              std::vector<LocalData<2>>& Aus) const 
{
    int mfields = us.size();
    int mx = pinfo->ns[0]; 
    int my = pinfo->ns[1];


    ValVector<2> new_u(MPI_COMM_WORLD,pinfo->ns,1,us.size(),1);
    ValVector<2> new_Au(MPI_COMM_WORLD,pinfo->ns,1,us.size(),1);
    auto new_us = new_u.getLocalDatas(0);
    auto new_Aus = new_Au.getLocalDatas(0);
    if(pinfo->hasNbr(Side<2>::west())){
        for(int field=0; field<mfields; field++){
            for(int j=0; j < my; j++){
                new_us[field][{-1,j}]=us[field][{-1,j}]+us[field][{0,j}];
            }
        }
    }
    if(pinfo->hasNbr(Side<2>::east())){
        for(int field=0; field<mfields; field++){
            for(int j=0; j < my; j++){
                new_us[field][{my,j}]=us[field][{my,j}]+us[field][{my-1,j}];
            }
        }
    }
    if(pinfo->hasNbr(Side<2>::south())){
        for(int field=0; field<mfields; field++){
            for(int i=0; i < mx; i++){
                new_us[field][{i,-1}]=us[field][{i,-1}]+us[field][{i,0}];
            }
        }
    }
    if(pinfo->hasNbr(Side<2>::north())){
        for(int field=0; field<mfields; field++){
            for(int i=0; i < mx; i++){
                new_us[field][{i,mx}]=us[field][{i,mx}]+us[field][{i,mx-1}];
            }
        }
    }
    
    /* Call patch operator to set boundary conditions correctly */
    applySinglePatch(pinfo,new_us,new_Aus,false);


    //modify rhs
    for(int field=0; field<mfields; field++){
        for(int j=0; j < my; j++){
            for(int i=0; i < mx; i++){
                Aus[field][{i,j}]-=new_Aus[field][{i,j}];
            }
        }
    }
}

const int *sgn::getS() const{
    return s;
}