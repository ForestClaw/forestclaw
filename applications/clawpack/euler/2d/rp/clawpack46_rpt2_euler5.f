      subroutine clawpack46_rpt2_euler5(ixy,maxm,meqn,mwaves,mbc,mx,
     &      ql,qr,aux1,aux2,aux3,ilr,asdq,bmasdq,bpasdq,maux)

      implicit none

      integer ixy, ilr, maxm, meqn,mwaves,mbc,mx, maux

      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision  asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision  aux1(1-mbc:maxm+mbc, maux)
      double precision  aux2(1-mbc:maxm+mbc, maux)
      double precision  aux3(1-mbc:maxm+mbc, maux)

      double precision delta(4), gamma, gamma1
      double precision rhsqrtl, rhsqrtr, pl,pr,rhsq2
      double precision a1,a2,a3,a4
      double precision waveb(5,4),sb(4)

      integer mu,mv,i, mw, m

      integer maxm2
      parameter (maxm2 = 520)

      double precision u2v2(-1:maxm2),
     &      u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),
     &      a(-1:maxm2),
     &      g1a2(-1:maxm2),euv(-1:maxm2)

      common /cparam/  gamma,gamma1


      if (-1 .gt. 1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'rpt (rpt2eu5.f) : need to increase maxm2'
         stop
      endif
c
      if (ixy .eq. 1) then
         mu = 2
         mv = 3
      else
         mu = 3
         mv = 2
      endif

      do i = 2-mbc, mx+mbc
         rhsqrtl = dsqrt(qr(i-1,1))
         rhsqrtr = dsqrt(ql(i,1))
         pl = gamma1*(qr(i-1,4) - 0.5d0*(qr(i-1,2)**2 +
     &         qr(i-1,3)**2)/qr(i-1,1))
         pr = gamma1*(ql(i,4) - 0.5d0*(ql(i,2)**2 +
     &         ql(i,3)**2)/ql(i,1))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (qr(i-1,mu)/rhsqrtl + ql(i,mu)/rhsqrtr) / rhsq2
         v(i) = (qr(i-1,mv)/rhsqrtl + ql(i,mv)/rhsqrtr) / rhsq2
         enth(i) = (((qr(i-1,4)+pl)/rhsqrtl
     &         + (ql(i,4)+pr)/rhsqrtr)) / rhsq2
         u2v2(i) = u(i)**2 + v(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2(i))
         a(i) = dsqrt(a2)
         g1a2(i) = gamma1 / a2
         euv(i) = enth(i) - u2v2(i)
      enddo


      do i = 2-mbc, mx+mbc
         a3 = g1a2(i) * (euv(i)*asdq(i,1)
     &         + u(i)*asdq(i,mu) + v(i)*asdq(i,mv) - asdq(i,4))
         a2 = asdq(i,mu) - u(i)*asdq(i,1)
         a4 = (asdq(i,mv) + (a(i)-v(i))*asdq(i,1) - a(i)*a3)
     &         / (2.d0*a(i))
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

         waveb(1,4) = 0.d0
         waveb(mu,4) = 0.d0
         waveb(mv,4) = 0.d0
         waveb(4,4) = 0.d0
         waveb(5,4) = asdq(i,5)
         sb(4) = v(i)

c
c        # compute the flux differences bmasdq and bpasdq
c
         do m=1,meqn
            bmasdq(i,m) = 0.d0
            bpasdq(i,m) = 0.d0
            do mw=1,4
               if (sb(mw) .lt. 0.d0) then
                  bmasdq(i,m) = bmasdq(i,m) + sb(mw) * waveb(m,mw)
               else
                  bpasdq(i,m) = bpasdq(i,m) + sb(mw) * waveb(m,mw)
               endif
            end do
         end do
        end do
c
c
         return
         end
