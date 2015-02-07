c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the Euler equations
c     #  with a tracer variable.
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were 
c     # computed in rpn2eu and stored in the common block comroe.
c
      dimension     ql(1-mbc:maxm+mbc, meqn)
      dimension     qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
c
      common /param/  gamma,gamma1
      dimension waveb(5,4),sb(4)
      parameter (maxm2 = 802)  !# assumes at most 600x600 grid with mbc=2
      common /comroe/ u2v2(-1:maxm2),
     &       u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2),
     &       g1a2(-1:maxm2),euv(-1:maxm2) 
c
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpt'
         stop
         endif
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
         do 20 i = 2-mbc, mx+mbc
            a3 = g1a2(i) * (euv(i)*asdq(i,1) 
     &             + u(i)*asdq(i,mu) + v(i)*asdq(i,mv) - asdq(i,4))
            a2 = asdq(i,mu) - u(i)*asdq(i,1)
            a4 = (asdq(i,mv) + (a(i)-v(i))*asdq(i,1) - a(i)*a3)
     &              / (2.d0*a(i))
            a1 = asdq(i,1) - a3 - a4
c
            waveb(1,1) = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*(v(i)-a(i))
            waveb(4,1) = a1*(enth(i) - v(i)*a(i))
            waveb(5,1) = 0.d0
            sb(1) = v(i) - a(i)
c
            waveb(1,2) = a3
            waveb(mu,2) = a3*u(i) + a2
            waveb(mv,2) = a3*v(i)
            waveb(4,2) = a3*0.5d0*u2v2(i) + a2*u(i)
            waveb(5,2) = 0.d0
            sb(2) = v(i)
c
            waveb(1,3) = a4
            waveb(mu,3) = a4*u(i)
            waveb(mv,3) = a4*(v(i)+a(i))
            waveb(4,3) = a4*(enth(i)+v(i)*a(i))
            waveb(5,3) = 0.d0
            sb(3) = v(i) + a(i)
c
            waveb(1,4) = 0.d0
            waveb(mu,4) = 0.d0
            waveb(mv,4) = 0.d0
            waveb(4,4) = 0.d0
            waveb(5,4) = asdq(i,5)
            sb(4) = v(i)
c
c           # compute the flux differences bmasdq and bpasdq
c
            do 10 m=1,meqn
               bmasdq(i,m) = 0.d0
               bpasdq(i,m) = 0.d0
               do 10 mw=1,4
                  bmasdq(i,m) = bmasdq(i,m) 
     &                         + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                  bpasdq(i,m) = bpasdq(i,m)
     &                         + dmax1(sb(mw), 0.d0) * waveb(m,mw)
   10             continue
c                 
   20       continue
c
      return
      end
