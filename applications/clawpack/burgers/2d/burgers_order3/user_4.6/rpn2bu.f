
c     =====================================================
      subroutine clawpack46_rpn2bu(ixy,maxm,meqn,mwaves,mbc,mx,ql,
     &   qr,auxl,auxr, wave,s,amdq,apdq)
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
      dimension    ql(1-mbc:maxm+mbc, meqn)
      dimension    qr(1-mbc:maxm+mbc, meqn)
      dimension    auxl(1-mbc:maxm+mbc,2)
      dimension    auxr(1-mbc:maxm+mbc,2)


c      dimension  qpx(1-mbc:maxm+mbc, meqn)
      dimension    f(1-mbc:maxm+mbc, meqn)
c      dimension   fp(1-mbc:maxm+mbc, meqn)
c      dimension  fpx(1-mbc:maxm+mbc, meqn)

      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)

      double precision favgl, favgr

      logical efix

c
c
c     # x- and y- Riemann problems are identical, so it doesn't matter if
c     # ixy=1 or 2.
c
      efix = .false.

      do 30 i=3-mbc,mx+mbc-2
c
c        # Compute the wave and speed
c
         favgr = auxl(i,2)
         favgl = auxr(i-1,2)
         wave(i,1,1) = favgr - favgl
         s(i,1) = 0.5*(qr(i-1,1)+ql(i,1))
         
c        # compute left-going and right-going flux differences:
c        ------------------------------------------------------
         amdq(i,1) = -dmin1(sign(1d0,s(i,1)), 0.d0) * wave(i,1,1)
         apdq(i,1) =  dmax1(sign(1d0,s(i,1)), 0.d0) * wave(i,1,1)

         if (efix) then
c           # entropy fix for transonic rarefactions:
c            if (q(i-1,1).lt.0.d0 .and. q(i,1).gt.0.d0) then
c               amdq(i,1) = - f(i-1,1)  
c               apdq(i,1) =   f(i,1)    
c               endif
c            endif
   30   continue
c
      return
      end
