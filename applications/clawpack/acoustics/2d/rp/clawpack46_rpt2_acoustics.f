c
c
c     =====================================================
      subroutine clawpack46_rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &      ql,qr,aux1,aux2,aux3,
     &      imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the acoustics equations.
c
c     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
c
      integer ixy, maxm, meqn, mwaves, mbc, mx, imp
      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision aux1(1-mbc:maxm+mbc, *)
      double precision aux2(1-mbc:maxm+mbc, *)
      double precision aux3(1-mbc:maxm+mbc, *)

c
c     # density, bulk modulus, and sound speed, and impedence of medium:
c     # (should be set in setprob.f)
      double precision rho, bulk, cc, zz
      common /cparam/ rho,bulk,cc,zz

      integer mu, mv, i
      double precision a1, a2

c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
      do i = 2-mbc, mx+mbc
         a1 = (-asdq(i,1) + zz*asdq(i,mv)) / (2.d0*zz)
         a2 = ( asdq(i,1) + zz*asdq(i,mv)) / (2.d0*zz)
c
c        # The down-going flux difference bmasdq is the product  -c * wave
c
         bmasdq(i,1)  = cc * a1*zz
         bmasdq(i,mu) = 0.d0
         bmasdq(i,mv) = -cc * a1
c
c        # The up-going flux difference bpasdq is the product  c * wave
c
         bpasdq(i,1)  = cc * a2*zz
         bpasdq(i,mu) = 0.d0
         bpasdq(i,mv) = cc * a2
c
      end do
c
      return
      end
