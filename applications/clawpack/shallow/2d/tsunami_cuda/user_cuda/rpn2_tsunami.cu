#include "../tsunami_user.h"

#include <fc2d_cudaclaw.h>
//#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */
//#include <cassert>

#include <fc2d_cudaclaw_check.cu>

__constant__ double s_grav;
__constant__ double s_dry_tolerance;
__constant__ double s_sea_level;


void geoclaw_setprob_cuda(double grav, double dry_tolerance, double sea_level)
{
    CHECK(cudaMemcpyToSymbol(s_grav,           &grav,          sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_dry_tolerance,  &dry_tolerance, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_sea_level,      &sea_level,     sizeof(double)));
}


__device__ void rpn2_geoclaw(int idir, int meqn, int mwaves, 
                                 int maux, double ql[], double qr[], 
                                 double auxl[], double auxr[],
                                 double fwave[], double s[], 
                                 double amdq[], double apdq[])
{
#if 0
    assert(mwaves == 3);
    assert(meqn == 3);
#endif    

    int mu = 1+idir;
    int mv = 2-idir;    

    double hl  = ql[0];
    double hr  = qr[0];

#if 0
    assert(hl > 0);
    assert(hr > 0);
#endif    

    double hul = ql[mu];
    double hur = qr[mu];

    double hvl = ql[mv];
    double hvr = qr[mv];

    double bl  = auxl[0];
    double br  = auxr[0];

    double ur = hur/hr;
    double vr = hvr/hr;

    double ul = hul/hl;
    double vl = hvl/hl;

    /* Compute fwaves */
    {
        /* Average h at cell interface  */
        double hbar = (hr + hl)/2;


        /* 1 wave speed of left state */
        double sl = ul - sqrt(s_grav*hl);    

        /* 2 wave speed of right state */
        double sr = ur + sqrt(s_grav*hr);    

        /* Roe speeds */
        double uhat = (sqrt(s_grav*hl)*ul + 
                       sqrt(s_grav*hr)*ur)/(sqrt(s_grav*hr) + sqrt(s_grav*hl));
        double chat = sqrt(s_grav*(hr + hl))/2;

        /* Compute wave speeds by comparing speeds above */
        s[0] = min(sl, uhat - chat);
        s[2] = max(sr, uhat + chat);
        s[1] = (s[0] + s[2])/2;
        
        /* Flux differences */
        double fluxdiff[3];
        double phir = s_grav*hr*hr/2 + hur*hur/hr;
        double phil = s_grav*hl*hl/2 + hul*hul/hl;

        fluxdiff[0] = (hr*ur) - (hl*ul);
        fluxdiff[1] = phir - phil + s_grav*hbar*(br - bl);
        fluxdiff[2] = hr*ur*vr - hl*ul*vl;

        /* Wave strengths */
        double beta[3];
        beta[0] = (s[2]*fluxdiff[0] - fluxdiff[1])/(s[2] - s[0]);
        beta[2] = (fluxdiff[1] - s[0]*fluxdiff[0])/(s[2] - s[0]);
        beta[1] = fluxdiff[2] - beta[0]*vl - beta[2]*vr;

        /* Flux waves = beta*R */
        fwave[0]  = beta[0];
        fwave[mu] = beta[0]*s[0];
        fwave[mv] = beta[0]*vl;

        fwave[3]    = 0;
        fwave[3+mu] = 0;
        fwave[3+mv] = beta[1];

        fwave[6]    = beta[2];
        fwave[6+mu] = beta[2]*s[2];
        fwave[6+mv] = beta[2]*vr;    
    }
    
    for(int mq = 0; mq < meqn; mq++)
    {
        /* z = copysign(x,y) : abs(z) = abs(x); sign(z) = sign(y) 
           if y == 0, returns abs(x) (not 0)  Should be okay here */

        int z[3];
        z[0] = (int) copysign(1.,s[0]); 
        z[1] = (int) copysign(1.,s[1]); 
        z[2] = (int) copysign(1.,s[2]); 

        /* (1-z)/2 : (-1 --> 1);  (1 --> 0) */
        amdq[mq]  = (1-z[0])/2*fwave[mq];
        amdq[mq] += (1-z[1])/2*fwave[meqn + mq];
        amdq[mq] += (1-z[2])/2*fwave[2*meqn + mq];

        /* (z+1)/2 : (-1 --> 0);  (1 --> 1)  */
        apdq[mq]  = (z[0]+1)/2*fwave[mq];
        apdq[mq] += (z[1]+1)/2*fwave[meqn + mq];
        apdq[mq] += (z[2]+1)/2*fwave[2*meqn + mq];
    }
}

__device__ cudaclaw_cuda_rpn2_t geoclaw_rpn2 = rpn2_geoclaw;

void geoclaw_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, geoclaw_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (rpn2_geoclaw): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void rpt2_geoclaw(int idir, int meqn, int mwaves, int maux,
                             double ql[], double qr[], 
                             double aux1[], double aux2[], double aux3[],
                             int imp, double asdq[],
                             double bmasdq[], double bpasdq[])

{
    int mu = 1+idir;
    int mv = 2-idir;    

    double hl  = ql[0];
    double hr  = qr[0];
    double hul = ql[mu];
    double hur = qr[mu];
    double hvl = ql[mv];
    double hvr = qr[mv];

    double ul = hul/hl;
    double vl = hvl/hl;

    double ur = hur/hr;
    double vr = hvr/hr;

#if 0
    double eta = hr  + aux2[0];
    double topo1 = aux1[0];
    double topo3 = aux3[0];
#endif    

    /* Construct Roe variables */
    double hrs = sqrt(hr);
    double hls = sqrt(hl);

    double vhat = (vr*hrs)/(hrs + hls) + (vl*hls)/(hrs + hls);
    double uhat = (ur*hrs)/(hrs+hls) + (ul*hls)/(hrs + hls);
    double hhat = (hr + hl)/2;

    double roe1 = vhat - sqrt(s_grav*hhat);
    double roe3 = vhat + sqrt(s_grav*hhat);

    double s1l = vl - sqrt(s_grav*hl);
    double s3r = vr + sqrt(s_grav*hr);

    double s1 = min(roe1,s1l);
    double s3 = max(roe3,s3r);
    // double s2 = (s1 + s3)/2;
    double s2 = uhat;

    double s[3];
    s[0] = s1;
    s[1] = s2;
    s[2] = s3;

    double delf1 = asdq[0];
    double delf2 = asdq[mu];
    double delf3 = asdq[mv];

    double beta[3];
    beta[0] = (s3*delf1/(s3-s1))-(delf3/(s3-s1));
    beta[1] = -s2*delf1 + delf2;
    beta[2] = (delf3/(s3-s1))-(s1*delf1/(s3-s1));

    double wave[9];
    wave[0]  = 1;
    wave[mu] = s2;
    wave[mv] = s1;

    wave[meqn + 0]  = 0;
    wave[meqn + mu] = 1;
    wave[meqn + mv] = 0;

    wave[2*meqn + 0]  = 1;
    wave[2*meqn + mu] = s2;
    wave[2*meqn + mv] = s3;


    double szm[3], szp[3];
    for(int mw = 0; mw < 3; mw++)
    {
        int z = (int) copysign(1.,s[mw]);
        szm[mw] = (1-z)/2*beta[mw];
        szp[mw] = (1+z)/2*beta[mw];
    }

    for(int mq = 0; mq < meqn; mq++)
    {
        bmasdq[mq]  = szm[0]*waveb[mq];
        bmasdq[mq] += szm[1]*waveb[meqn + mq];
        bmasdq[mq] += szm[2]*waveb[2*meqn + mq];

        bpasdq[mq]  = szp[0]*waveb[mq];
        bpasdq[mq] += szp[1]*waveb[meqn + mq];
        bpasdq[mq] += szp[2]*waveb[2*meqn + mq];
    }
}



#if 0
__device__ void rpt2_geoclaw(int idir, int meqn, int mwaves, int maux,
                                 double ql[], double qr[], 
                                 double aux1[], double aux2[], double aux3[],
                                 int imp, double asdq[],
                                 double bmasdq[], double bpasdq[])
{

    double alpha1, alpha2, alpha3;
    double h, hsqrtr, hsqrtl, hsq2, u,v,a;
    double waveb[9], sb[3];

    int mu,mv;
    int mq;

    mu = 1+idir;
    mv = 2-idir;    

    h = (qr[0] + ql[0])/2.0;
    hsqrtr = sqrt(qr[0]);
    hsqrtl = sqrt(ql[0]);
    hsq2 = hsqrtl + hsqrtr;

    u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
    v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
    a = sqrt(s_grav*h);


    alpha1 = ((v+a)*asdq[0] - asdq[mv])/(2.0*a);
    alpha2 = asdq[mu] - u*asdq[0];
    alpha3 = (-(v-a)*asdq[0] + asdq[mv])/(2.0*a);

    waveb[0]  = alpha1;
    waveb[mu] = alpha1*u;
    waveb[mv] = alpha1*(v-a);
    sb[0] = v - a;

    waveb[meqn + 0]  = 0.0;
    waveb[meqn + mu] = alpha2;
    waveb[meqn + mv] = 0.0;
    sb[1] = v;

    waveb[2*meqn + 0]  = alpha3;
    waveb[2*meqn + mu] = alpha3*u;
    waveb[2*meqn + mv] = alpha3*(v+a);
    sb[2] = v + a;

    for(mq = 0; mq < meqn; mq++)
    {
        /* Loop-unrolling! loop over mwaves=3*/
        bmasdq[mq]  = min(sb[0],0.)*waveb[mq];
        bmasdq[mq] += min(sb[1],0.)*waveb[meqn + mq];
        bmasdq[mq] += min(sb[2],0.)*waveb[2*meqn + mq];

        bpasdq[mq]  = max(sb[0],0.)*waveb[mq];
        bpasdq[mq] += max(sb[1],0.)*waveb[meqn + mq];
        bpasdq[mq] += max(sb[2],0.)*waveb[2*meqn + mq];
    }

#endif


__device__ cudaclaw_cuda_rpt2_t geoclaw_rpt2 = rpt2_geoclaw;

void geoclaw_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, geoclaw_rpt2, 
                                          sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (rpt2_geoclaw): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}

