
ccc     =====================================================
c      subroutine clawpack46_rpn2adv(ixy,maxm,meqn,mwaves,mbc,mx,
c     &  ql,qr,auxl,auxr,wave,s,amdq,apdq)
cc     =====================================================

c
c     =====================================================
      subroutine clawpack46_rpn2bu(ixy,maxm,meqn,mwaves,mbc,mx,q,
     &   qp,qb,qt, wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for Burgers' equation in 2d:
c     #  u_t + 0.5*(u^2)_x + 0.5*(u^2)_y = 0
c     
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension    q(1-mbc:maxm+mbc, meqn)
      dimension    qp(1-mbc:maxm+mbc, meqn)
      dimension    qb(1-mbc:maxm+mbc,1)
      dimension    qt(1-mbc:maxm+mbc,1)


      dimension  qpx(1-mbc:maxm+mbc, meqn)
      dimension    f(1-mbc:maxm+mbc, meqn)
      dimension   fp(1-mbc:maxm+mbc, meqn)
      dimension  fpx(1-mbc:maxm+mbc, meqn)

      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)

      logical efix

c
c
c     # x- and y- Riemann problems are identical, so it doesn't matter if
c     # ixy=1 or 2.
c
      efix = .false.


c     compute point values of the flux     
      do 11 i=2-mbc, mx+mbc-1
         fp(i,1)=0.5d0 * qp(i,1)**2.

  11     continue
c
c     compute 3rd order accurate averaged fluxes          
      do 12 i=3-mbc, mx+mbc-2
         f(i,1)=fp(i,1)+(fp(i-1,1)-2d0*fp(i,1)+fp(i+1,1))/24d0
     &   + (0.5*qb(i,1)**2. - 2.*fp(i,1) + 0.5*qt(i,1)**2.)/24d0 

  12     continue
c
c
      do 30 i=3-mbc,mx+mbc-2
c
c        # Compute the wave and speed
c
         wave(i,1,1) = f(i,1) - f(i-1,1)
         s(i,1) = 0.5*(q(i-1,1)+q(i,1))
         
c        3rd order accurate wavespeed (without limiter)         
c        s(i,1) = (7d0*(q(i,1)+q(i-1,1))-q(i-2,1)-q(i+1,1))/12d0        

c
c
c        # compute left-going and right-going flux differences:
c        ------------------------------------------------------
         amdq(i,1) = -dmin1(sign(1d0,s(i,1)), 0.d0) * wave(i,1,1)
         apdq(i,1) =  dmax1(sign(1d0,s(i,1)), 0.d0) * wave(i,1,1)

         if (efix) then
c           # entropy fix for transonic rarefactions:
            if (q(i-1,1).lt.0.d0 .and. q(i,1).gt.0.d0) then
               amdq(i,1) = - f(i-1,1)  
               apdq(i,1) =   f(i,1)    
               endif
            endif
   30   continue
c
      return
      end
