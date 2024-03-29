#include "../radialdam_user.h"

#include <fc2d_cudaclaw.h>
//#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */
//#include <cassert>

#include <fc2d_cudaclaw_check.h>

__constant__ double s_grav;


void radialdam_setprob_cuda(double grav)
{
    CHECK(cudaMemcpyToSymbol(s_grav,  &grav, sizeof(double)));
}


__device__ void radialdam_rpn2shallow(int idir, int meqn, int mwaves, 
    								  int maux, double ql[], double qr[], 
                              	  	  double auxl[], double auxr[],
                              	  	  double wave[], double s[], 
                              	  	  double amdq[], double apdq[])
{
    //assert(mwaves == 3);
    //assert(meqn == 3);

	double h, hsqrtl, hsqrtr, hsq2, u, v, a;
	double delta1, delta2, delta3;
	double a1, a2, a3;

    int mu, mv;
    int mq;

    mu = 1+idir;
    mv = 2-idir;    

    /* This should be true ... */
    //assert(qr[0] > 0);
    //assert(ql[0] > 0);

	h = (qr[0] + ql[0])/2.0;
	hsqrtr = sqrt(qr[0]);
	hsqrtl = sqrt(ql[0]);
	hsq2 = hsqrtl + hsqrtr;

	u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
	v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
	a = sqrt(s_grav*h);
    // a = sqrt(abs(h));

	delta1 = qr[0 ] -ql[0 ];
	delta2 = qr[mu] -ql[mu];
	delta3 = qr[mv] -ql[mv];
	a1 = ( (u+a)*delta1 - delta2 )/(2*a);
	a2 = -v*delta1 + delta3;
	a3 = (-(u-a)*delta1 + delta2 )/(2*a);
    
    wave[0]  = a1;
    wave[mu] = a1*(u-a);
    wave[mv] = a1*v;
    s[0] = u-a;
    
    wave[0 +meqn] = 0.0;
    wave[mu+meqn] = 0.0;
    wave[mv+meqn] = a2;
    s[1] = u;
    
    wave[0 +meqn*2] = a3;
    wave[mu+meqn*2] = a3*(u+a);
    wave[mv+meqn*2] = a3*v;
    s[2] = u+a;
    
    for(mq = 0; mq < meqn; mq++)
    {
        /* Loop-unrolling! loop over mwaves=3*/
        amdq[mq]  = min(s[0],0.)*wave[mq];
        amdq[mq] += min(s[1],0.)*wave[meqn + mq];
        amdq[mq] += min(s[2],0.)*wave[2*meqn + mq];

        apdq[mq]  = max(s[0],0.)*wave[mq];
        apdq[mq] += max(s[1],0.)*wave[meqn + mq];
        apdq[mq] += max(s[2],0.)*wave[2*meqn + mq];

    }
}

__device__ cudaclaw_cuda_rpn2_t radialdam_rpn2 = radialdam_rpn2shallow;

void radialdam_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, radialdam_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radialdam_rpn2adv): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}


__device__ void radialdam_rpt2shallow(int idir, int meqn, int mwaves, int maux,
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
}


__device__ cudaclaw_cuda_rpt2_t radialdam_rpt2 = radialdam_rpt2shallow;

void radialdam_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, radialdam_rpt2, 
                                          sizeof(cudaclaw_cuda_rpt2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (radialdam_rpt2shallow): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}

