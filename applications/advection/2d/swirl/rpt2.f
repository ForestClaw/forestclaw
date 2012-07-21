c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the scalar equation
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c
      integer ixy, maxm, meqn, mwaves, mbc, mx, ilr

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision aux1(1-mbc:maxm+mbc,*)
      double precision aux2(1-mbc:maxm+mbc,*)
      double precision aux3(1-mbc:maxm+mbc,*)

      integer i
      double precision ubar, vbar, stran, stranm, stranp

      common /comrp/ ubar,vbar
c
c     # transverse wave speeds have been computed in rpn2
c     # maux=0 and aux arrays are unused in this example.
c
      if (ixy .eq. 1) then
         stran = vbar
      else
         stran = ubar
      endif

      stranm = dmin1(stran, 0.d0)
      stranp = dmax1(stran, 0.d0)

      do  i = 2-mbc, mx+mbc
         bmasdq(i,1) = stranm * asdq(i,1)
         bpasdq(i,1) = stranp * asdq(i,1)
      enddo
c
      return
      end
