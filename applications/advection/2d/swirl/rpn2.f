c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the sample scalar equation
c     #  q_t + u*q_x + v*q_y = 0
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q,
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c     # maux=0 and aux arrays are unused in this example.
c
c
      implicit none

      integer ixy, maxm, meqn, mwaves, mbc, mx
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision   s(1-mbc:maxm+mbc, mwaves)
      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, *)
      double precision auxr(1-mbc:maxm+mbc, *)

      integer i
      double precision ubar, vbar

      common /comrp/ ubar,vbar

      do i = 2-mbc, mx+mbc-1
         wave(i,1,1) = ql(i,1) - qr(i-1,1)
         if (ixy .eq. 1) then
            s(i,1) = ubar
         else
            s(i,1) = vbar
         endif

c        # flux differences:
         amdq(i,1) = dmin1(s(i,1), 0.d0) * wave(i,1,1)
         apdq(i,1) = dmax1(s(i,1), 0.d0) * wave(i,1,1)

      enddo

      end
