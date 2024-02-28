#include "../shockbubble_user.h"

#include <fc2d_cudaclaw.h>

#include <fc2d_cudaclaw_check.h>

__constant__ double s_gamma;
__constant__ double s_gamma1;


void setprob_cuda()
{
    double gamma;
    FILE *f = fopen("setprob.data","r");
    fscanf(f,"%lf",&gamma);
    fclose(f);    

    double gamma1 = gamma - 1.0;
    CHECK(cudaMemcpyToSymbol(s_gamma, &gamma, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(s_gamma1, &gamma1, sizeof(double)));
}


__device__ void shockbubble_rpn2euler4(int idir, int meqn, int mwaves, 
                                       int maux, double ql[], double qr[], 
                                       double auxl[], double auxr[],
                                       double wave[], double s[], 
                                       double amdq[], double apdq[])
{
	double rhsqrtl, rhsqrtr, pr, pl, rhsq2, u, v, a;
    double enth, u2v2, g1a2, euv, ainv2;
	double delta[4];
	double a1, a2, a3, a4;

    int mu, mv;
    int mq;

    mu = 1+idir;
    mv = 2-idir;    

    rhsqrtr = sqrt(qr[0]);
    rhsqrtl = sqrt(ql[0]);

    pl = s_gamma1*(ql[3] - 0.5*(ql[1]*ql[1] + ql[2]*ql[2])/ql[0]);
    pr = s_gamma1*(qr[3] - 0.5*(qr[1]*qr[1] + qr[2]*qr[2])/qr[0]);

    rhsq2 = rhsqrtl + rhsqrtr;
    u = (ql[mu]/rhsqrtl + qr[mu]/rhsqrtr) / rhsq2;
    v = (ql[mv]/rhsqrtl + qr[mv]/rhsqrtr) / rhsq2;
    enth = ((ql[3] + pl)/rhsqrtl + (qr[3] + pr)/rhsqrtr) / rhsq2;
    u2v2 = u*u + v*v;
    a2 = s_gamma1*(enth - 0.5*u2v2);
    a = sqrt(a2);
    ainv2 = 0.5*rsqrt(a2);
    g1a2 = s_gamma1 / a2;
    euv = enth - u2v2;

    delta[0] = qr[0 ] - ql[0 ];
    delta[1] = qr[mu] - ql[mu];
    delta[2] = qr[mv] - ql[mv];
    delta[3] = qr[3]  - ql[3];

    a3 = g1a2*(euv*delta[0] + u*delta[1] + v*delta[2] - delta[3]);
    a2 = delta[2] - v*delta[0];
    a4 = (delta[1] + (a - u)*delta[0] - a*a3)*ainv2;
    a1 = delta[0] - a3 - a4;

    /* Acoustic wave */
    wave[0]  = a1;
    wave[mu] = a1*(u-a);
    wave[mv] = a1*v;
    wave[3]  = a1*(enth - u*a);
    s[0] = u-a;
    
#if 0    
    /* Shear wave */
    wave[0 +meqn] = 0.0;
    wave[mu+meqn] = 0.0;
    wave[mv+meqn] = a2;
    wave[3 +meqn] = a2*v;
    s[1] = u;
    
    /* Entropy wave */
    wave[0 +meqn*2] = a3;
    wave[mu+meqn*2] = a3*u;
    wave[mv+meqn*2] = a3*v;
    wave[3 +meqn*2] = a3*0.5*u2v2;
    s[2] = u;
#endif

    /* Entropy  + shear wave */
    wave[0 +meqn] = a3;
    wave[mu+meqn] = a3*u;
    wave[mv+meqn] = a3*v + a2;
    wave[3 +meqn] = a3*0.5*u2v2 + a2*v;
    s[1] = u;


    /* Acoustic */
    wave[0 +meqn*2] = a4;
    wave[mu+meqn*2] = a4*(u + a);
    wave[mv+meqn*2] = a4*v;
    wave[3 +meqn*2] = a4*(enth + u*a);
    s[2] = u + a;
    
    //printf("s[0] = %12.4f; s[1] = %12.4f; s[2] = %12.4f; s[3] = %12.4f\n",
    //s[0], s[1], s[2], s[3]);

    /* Entropy fix not implemented yet */

    for(mq = 0; mq < meqn; mq++)
    {
        /* Loop-unrolling! loop over mwaves=4*/
        amdq[mq]  = min(s[0],0.)*wave[0*meqn + mq];
        amdq[mq] += min(s[1],0.)*wave[1*meqn + mq];
        amdq[mq] += min(s[2],0.)*wave[2*meqn + mq];
        //amdq[mq] += min(s[3],0.)*wave[3*meqn + mq];

        apdq[mq]  = max(s[0],0.)*wave[0*meqn + mq];
        apdq[mq] += max(s[1],0.)*wave[1*meqn + mq];
        apdq[mq] += max(s[2],0.)*wave[2*meqn + mq];
        //apdq[mq] += max(s[3],0.)*wave[3*meqn + mq];

    }
}

__device__ cudaclaw_cuda_rpn2_t shockbubble_rpn2 = shockbubble_rpn2euler4;

void shockbubble_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, shockbubble_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radialdam_rpn2adv): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void shockbubble_rpt2euler4(int idir, int meqn, int mwaves, int maux,
                                      double ql[], double qr[], 
                                      double aux1[], double aux2[], double aux3[],
                                      int imp, double asdq[],
                                      double bmasdq[], double bpasdq[])
{
    int mu = 1+idir;
    int mv = 2-idir;    

    /* First get (u,v), etc at normal interface */
    double rhsqrtr = sqrt(qr[0]);
    double rhsqrtl = sqrt(ql[0]);

    double pl = s_gamma1*(ql[3] - 0.5*(ql[1]*ql[1] + ql[2]*ql[2])/ql[0]);
    double pr = s_gamma1*(qr[3] - 0.5*(qr[1]*qr[1] + qr[2]*qr[2])/qr[0]);

    double rhsq2 = rhsqrtl + rhsqrtr;
    double u = (ql[mu]/rhsqrtl + qr[mu]/rhsqrtr) / rhsq2;
    double v = (ql[mv]/rhsqrtl + qr[mv]/rhsqrtr) / rhsq2;
    double enth = ((ql[3] + pl)/rhsqrtl + (qr[3] + pr)/rhsqrtr) / rhsq2;
    double u2v2 = u*u + v*v;
    double a2 = s_gamma1*(enth - 0.5*u2v2);
    double a = sqrt(a2);
    double ainv2 = 0.5*rsqrt(a2);
    double g1a2 = s_gamma1 / a2;
    double euv = enth - u2v2;

    double delta[4];
    delta[0] = asdq[0];
    delta[1] = asdq[mu];
    delta[2] = asdq[mv];
    delta[3] = asdq[3];
    //a3 = g1a2(i) * (euv(i)*asdq(i,1) + u(i)*asdq(i,mu) + v(i)*asdq(i,mv) - asdq(i,4))
    double a3 = g1a2*(euv*delta[0] + u*delta[1] + v*delta[2] - delta[3]);
    // a2 = asdq(i,mu) - u(i)*asdq(i,1)  /* From fortran */
    a2 = delta[1] - u*delta[0];
    // a4 = (asdq(i,mv) + (a(i)-v(i))*asdq(i,1) - a(i)*a3) / (2.d0*a(i))
    double a4 = (delta[2] + (a - v)*delta[0] - a*a3)*ainv2;

    // a1 = asdq(i,1) - a3 - a4
    double a1 = delta[0] - a3 - a4;

    /* Acoustic wave */
    double wave[12], s[3];
    wave[0]  = a1;
    wave[mu] = a1*u;
    wave[mv] = a1*(v-a);
    wave[3]  = a1*(enth - v*a);
    s[0] = v - a;
    
#if 0
    /* Shear wave */
    wave[0 +meqn] = 0.0;
    wave[mu+meqn] = a2;
    wave[mv+meqn] = 0.0;
    wave[3 +meqn] = a2*u;
    s[1] = v;
    
    /* Entropy wave */
    wave[0 +meqn*2] = a3;
    wave[mu+meqn*2] = a3*u;
    wave[mv+meqn*2] = a3*v;
    wave[3 +meqn*2] = a3*0.5*u2v2;
    s[2] = v;
#endif

#if 1
    /* Entropy and shear waves (combine to save memory) */
    wave[0 +meqn] = a3;
    wave[mu+meqn] = a3*u        + a2;
    wave[mv+meqn] = a3*v;
    wave[3 +meqn] = a3*0.5*u2v2 + a2*u;
    s[1] = v;
#endif    



    /* Acoustic */
    wave[0 +meqn*2] = a4;
    wave[mu+meqn*2] = a4*u;
    wave[mv+meqn*2] = a4*(v+a);
    wave[3 +meqn*2] = a4*(enth + v*a);
    s[2] = v + a;
    
    /* Entropy fix not implemented yet */

    for(int mq = 0; mq < meqn; mq++)
    {
        /* Loop-unrolling! loop over mwaves=4*/
        bmasdq[mq]  = min(s[0],0.)*wave[mq];
        bmasdq[mq] += min(s[1],0.)*wave[meqn + mq];
        bmasdq[mq] += min(s[2],0.)*wave[2*meqn + mq];
        //bmasdq[mq] += min(s[3],0.)*wave[3*meqn + mq];

        bpasdq[mq]  = max(s[0],0.)*wave[mq];
        bpasdq[mq] += max(s[1],0.)*wave[meqn + mq];
        bpasdq[mq] += max(s[2],0.)*wave[2*meqn + mq];
        //bpasdq[mq] += max(s[3],0.)*wave[3*meqn + mq];
    }
}


__device__ cudaclaw_cuda_rpt2_t shockbubble_rpt2 = shockbubble_rpt2euler4;

void shockbubble_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, shockbubble_rpt2, 
                                          sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (shockbubble_rpt2euler4): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void shockbubble_src2_cuda(int meqn, int maux, 
                                      double xlower, double ylower, 
                                      double dx, double dy, 
                                      double q[], double aux[], 
                                      double t, double dt, 
                                      int i, int j)

{
    double dt2 = dt/2;

    double rad = ylower + (j-0.5)*dy;
    double gamma1 = s_gamma1;

    /* Stage 1 */
    double rho = q[0];
    double u = q[1]/rho;
    double v = q[2]/rho;
    double E = q[3];
    double press = gamma1*(E - 0.5*rho*(u*u + v*v));

    double qstar[4];

    qstar[0] = q[0] - dt2*q[2]/rad;
    qstar[1] = q[1] - dt2*rho*u*v/rad;
    qstar[2] = q[2] - dt2*rho*v*v/rad;
    qstar[3] = q[3] - dt2*v*(E + press)/rad;

    /* Stage 2 */
    rho = qstar[0];
    u = qstar[1]/rho;
    v = qstar[2]/rho;
    E = qstar[3];
    press = gamma1*(E - 0.5*rho*(u*u + v*v));
    q[0] = q[0] - dt*qstar[2]/rad;
    q[1] = q[1] - dt*rho*u*v/rad;
    q[2] = q[2] - dt*rho*v*v/rad;
    q[3] = q[3] - dt*v*(E + press)/rad;
}

__device__ cudaclaw_cuda_src2_t shockbubble_src2 = shockbubble_src2_cuda;

void shockbubble_assign_src2(cudaclaw_cuda_src2_t *src2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(src2, shockbubble_src2, 
                                          sizeof(cudaclaw_cuda_src2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (shockbubble_src2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}



