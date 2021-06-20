      subroutine clawpack46_rpt2_acoustics_vc(ixy,maxm,meqn,
     &      mwaves,mbc,mx,ql,qr,
     &      aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)
      implicit none
c
c     # Riemann solver in the transverse direction for the acoustics equations
c     # with varying material properties rho and kappa
c
c     # auxN(i,1) holds rho
c     # auxN(i,2) holds c
c     #  N = 1 for row below
c     #      2 for this row
c     #      3 for row above
c
c     # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
c
c     # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq
c
      integer ixy, maxm, meqn, mwaves, mbc, mx, imp
      double precision    ql(1-mbc:maxm+mbc, meqn)
      double precision    qr(1-mbc:maxm+mbc, meqn)
      double precision    asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)

      integer i, i1, mu, mv
      double precision c, cm, cp, zm, zz, zp, a1, a2

      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c
      do  i = 2-mbc, mx+mbc
c
c        # imp is used to flag whether wave is going to left or right,
c        # since material properties are different on the two sides
c
         if (imp.eq.1) then
c            # asdq = amdq, moving to left
             i1 = i-1
           else
c            # asdq = apdq, moving to right
             i1 = i
           endif
c
c        # The flux difference asdq is split into downward moving part
c        # traveling at speed -c relative to the medium below and
c        # an upward moving part traveling
c        # at speed +c relative to the medium above.
c
c        # Note that the sum of these parts does not give all of asdq
c        # since there is also reflection at the interfaces which decreases
c        # the flux.
c
c        # sound speed in each row of cells:
         cm = aux1(i1,2)
         c = aux2(i1,2)
         cp = aux3(i1,2)
c
c        # impedances:
         zm = aux1(i1,1)*aux1(i1,2)
         zz = aux2(i1,1)*aux2(i1,2)
         zp = aux3(i1,1)*aux3(i1,2)

c        # transmitted part of down-going wave:
         a1 = (-asdq(i,1) + asdq(i,mv)*zz) /
     &         (zm + zz)

c        # transmitted part of up-going wave:
         a2 = (asdq(i,1) + asdq(i,mv)*zz) /
     &         (zz + zp)
c
c        # The down-going flux difference bmasdq is the product  -c * wave
c
         bmasdq(i,1) = cm * a1*zm
         bmasdq(i,mu) = 0.d0
         bmasdq(i,mv) = -cm * a1
c
c        # The up-going flux difference bpasdq is the product  c * wave
c
         bpasdq(i,1) = cp * a2*zp
         bpasdq(i,mu) = 0.d0
         bpasdq(i,mv) = cp * a2

      end do

      return
      end
