      subroutine clawpack46_rpt2_euler4(ixy,maxm,meqn,mwaves,mbc,
     &      mx,ql,qr,aux1,aux2,aux3,ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the Euler equations.
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were
c     # computed in rpn2eu and stored in the common block comroe.
c
      integer ixy, maxm, meqn, mwaves, mbc, mx, ilr
      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision aux1(1-mbc:maxm+mbc, *)
      double precision aux2(1-mbc:maxm+mbc, *)
      double precision aux3(1-mbc:maxm+mbc, *)

c
      double precision gamma, gamma1
      common /cparam/  gamma,gamma1


      double precision waveb(4,3),sb(3)

      integer maxm2
      parameter (maxm2 = 602)  !# assumes at most 200x200 grid with mbc=2

      double precision u2v2, u,v,enth,a,g1a2, euv
      common /comroe/ u2v2(-1:maxm2),
     &       u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2),
     &       g1a2(-1:maxm2),euv(-1:maxm2)


      integer mu, mv, i, m, mw
      double precision a1, a2, a3, a4
c
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpt2'
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
         do i = 2-mbc, mx+mbc
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
            sb(1) = v(i) - a(i)
c
            waveb(1,2) = a3
            waveb(mu,2) = a3*u(i) + a2
            waveb(mv,2) = a3*v(i)
            waveb(4,2) = a3*0.5d0*u2v2(i) + a2*u(i)
            sb(2) = v(i)
c
            waveb(1,3) = a4
            waveb(mu,3) = a4*u(i)
            waveb(mv,3) = a4*(v(i)+a(i))
            waveb(4,3) = a4*(enth(i)+v(i)*a(i))
            sb(3) = v(i) + a(i)
c
c           # compute the flux differences bmasdq and bpasdq
c
            do m=1,meqn
               bmasdq(i,m) = 0.d0
               bpasdq(i,m) = 0.d0
               do mw=1,3
                  bmasdq(i,m) = bmasdq(i,m)
     &                         + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                  bpasdq(i,m) = bpasdq(i,m)
     &                         + dmax1(sb(mw), 0.d0) * waveb(m,mw)
               end do
            end do
        end do
c
      return
      end
