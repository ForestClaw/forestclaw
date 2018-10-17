#include "../radialdam_user.h"

#include <fc2d_cudaclaw.h>
#include <fclaw_base.h>  /* Needed for SC_MIN, SC_MAX */
#include <cassert>

static __device__ double grav = 1.0;

__device__ void radialdam_rpn2shallow(int idir, int meqn, int mwaves, 
    								  int maux, double ql[], double qr[], 
                              	  	  double auxl[], double auxr[],
                              	  	  double wave[], double s[], 
                              	  	  double amdq[], double apdq[])
{
    /* wave[mwaves][meqn] */
    /* idir in 0,1 : needed to get correct  */
    assert(mwaves == 3);
    assert(meqn == 3);

	int mu, mv;

	if (idir == 0)
	{
		mu = 1;
		mv = 2;
	}
	else
	{
		mu = 2;
		mv = 1;
	}

	double h, hsqrtl, hsqrtr, hsq2, u, v, a;
	double delta1, delta2, delta3;
	double a1, a2, a3;

    assert(qr[0] > 0);
    assert(ql[0] > 0);

	h = 0.5*(qr[0] + ql[0]);
	hsqrtr = sqrt(qr[0]);
	hsqrtl = sqrt(ql[0]);
	hsq2 = hsqrtl + hsqrtr;

	u = (qr[mu]/hsqrtr + ql[mu]/hsqrtl) / hsq2;
	v = (qr[mv]/hsqrtr + ql[mv]/hsqrtl) / hsq2;
	a = sqrt(grav*h);
    // a = sqrt(abs(h));

	delta1 = qr[0 ] -ql[0 ];
	delta2 = qr[mu] -ql[mu];
	delta3 = qr[mv] -ql[mv];
	a1 = ( (u+a)*delta1 - delta2 )*0.5/a;
	a2 = -v*delta1 + delta3;
	a3 = (-(u-a)*delta1 + delta2 )*0.5/a;
    
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
    
    for(int mq = 0; mq < meqn; mq++)
    {
    	amdq[mq] = 0.0;
    	apdq[mq] = 0.0;
    	for (int mw = 0; mw < mwaves; mw++)
    	{
    		amdq[mq] += SC_MIN(s[mw],0) * wave[mq+mw*meqn];
        	apdq[mq] += SC_MAX(s[mw],0) * wave[mq+mw*meqn]; 
    	    // amdq[mq] += 0;
            // apdq[mq] += 0; 
        }
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