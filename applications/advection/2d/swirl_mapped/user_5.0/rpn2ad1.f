c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
     &                  auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y = 0
c     # where u and v are a given velocity field.
c
      implicit none

      integer ixy,maux,maxm, meqn, mwaves, mbc, mx
      double precision wave(meqn,mwaves, 1-mbc:maxm+mbc)
      double precision    s(mwaves, 1-mbc:maxm+mbc)
      double precision   ql(meqn, 1-mbc:maxm+mbc)
      double precision   qr(meqn, 1-mbc:maxm+mbc)
      double precision apdq(meqn, 1-mbc:maxm+mbc)
      double precision amdq(meqn, 1-mbc:maxm+mbc)
      double precision auxl(maux, 1-mbc:maxm+mbc)
      double precision auxr(maux,1-mbc:maxm+mbc)

      integer i, m


c     # Set wave, speed, and flux differences:
c     ------------------------------------------

      do i = 2-mbc, mx+mbc
         do m = 1,meqn
            wave(m,1,i) = ql(m,i) - qr(m,i-1)
         enddo
         s(1,i) = auxl(ixy,i)

c        # The flux difference df = s*wave  all goes in the downwind
c        direction:
         do m = 1,meqn
            amdq(m,i) = min(s(1,i), 0.d0) * wave(m,1,i)
            apdq(m,i) = max(s(1,i), 0.d0) * wave(m,1,i)
         enddo
      enddo

      return
      end
