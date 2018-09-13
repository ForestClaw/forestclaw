#include "../fc2d_cudaclaw5.h"
#include "../fc2d_cudaclaw5_fort.h"
#include "cudaclaw5_update_q.h"

void cudaclaw5_step2_wrap(int maxm, int meqn, int maux, int mbc,
                          int method[], int mthlim[], int mcapa, int mwaves, 
                          int mx, int my, double qold[], double aux[],
                          double dx, double dy, double dt, double cfl, 
                          double work[], int mwork, double xlower, 
                          double ylower, int  level, double t, 
                          double fp[], double fm[],
                          double gp[], double gm[], 
                          cudaclaw5_fort_rpn2_t rpn2, 
                          cudaclaw5_fort_rpt2_t rpt2,
                          cudaclaw5_fort_flux2_t flux2, 
                          int block_corner_count[], int* ierror)
{

#if 0
c     # Local variables
      integer i0faddm, i0faddp, i0gaddm, i0gaddp
      integer i0q1d, i0dtdx1, i0dtdy1
      integer i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
      integer i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr, i0auxl
      integer i0auxr

      integer i,j,m

c     Needed by Riemann solvers.  This should be fixed later by a 'context'
c     for a Riemann solver.
      double precision dtcom, dxcom,dycom,tcom
      integer icom, jcom
      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

c     # This should be set to actual time, in case the user wants it
c     # it for some reason in the Riemann solver. 

c     # Set up work arrays (these are not used yet) 

      i0faddm = 1
      i0faddp = i0faddm +   (maxm+2*mbc)*meqn
      i0gaddm = i0faddp +   (maxm+2*mbc)*meqn
      i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
      i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d   +   (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1 +   (maxm+2*mbc)
      i0aux1  = i0dtdy1 +   (maxm+2*mbc)
      i0aux2  = i0aux1  +   (maxm+2*mbc)*maux
      i0aux3  = i0aux2  +   (maxm+2*mbc)*maux      
c
c
      i0next  = i0aux3  + (maxm+2*mbc)*maux    !# next free space
      mused   = i0next - 1                    !# space already used
      mwork1  = mwork - mused           !# remaining space (passed to step2)

      if (mused.gt.mwork) then
         ierror = 1
         return
      endif


c     # Include four additional arguments to avoid need for
c     # global array
c      call cudaclaw5_step2(maxm,maxmx,maxmy,meqn,maux, mbc,
c     &      mx,my, qold,aux,dx,dy,dt,
c     &      cfl,fm,fp,gm,gp,
c     &      work(i0faddm),work(i0faddp),
c     &      work(i0gaddm),work(i0gaddp),
c     &      work(i0q1d),work(i0dtdx1),work(i0dtdy1),
c     &      work(i0aux1),work(i0aux2),work(i0aux3),
c     &      work(i0next),mwork1,rpn2,rpt2,flux2,
c     &      mwaves,mcapa,method,mthlim,block_corner_count,ierror)
#endif

    double dtdx, dtdy;

    *ierror = 0;

    CUDACLAW5_STEP2(&maxm,&meqn,&maux,&mbc,&mx,&my,qold,aux,
                     &dx, &dy, &dt, &cfl, fm,fp,gm,gp,rpn2,rpt2,
                     block_corner_count, ierror);

    /* # update q */
    dtdx = dt/dx;
    dtdy = dt/dy;
#if 1    
    CUDACLAW5_FORT_UPDATE_Q(&meqn, &mx, &my, &mbc, &maux,
                           &dtdx, &dtdy,qold,fp,fm,
                           gp, gm, &mcapa);
#else    
    cudaclaw5_update_q(meqn,mx,my,mbc,dtdx,dtdy,qold,
                       fm,fp,gm,fp,mcapa);
#endif      

}