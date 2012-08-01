c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                  auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y = 0
c     # where u and v are a given velocity field.
c
      implicit none

      integer ixy,maxm, meqn, mwaves, mbc, mx
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, *)
      double precision auxr(1-mbc:maxm+mbc, *)

      integer i


c     # Set wave, speed, and flux differences:
c     ------------------------------------------

      do i = 2-mbc, mx+mbc
         wave(i,1,1) = ql(i,1) - qr(i-1,1)
         s(i,1) = auxl(i,ixy)

c        # The flux difference df = s*wave  all goes in the downwind direction:
         amdq(i,1) = dmin1(auxl(i,ixy), 0.d0) * wave(i,1,1)
         apdq(i,1) = dmax1(auxl(i,ixy), 0.d0) * wave(i,1,1)
      enddo

      return
      end
