      subroutine clawpack46_rpn2_euler5(ixy,maxm,meqn,
     &      mwaves,mbc,mx,ql,qr,auxl,auxr,
     &      wave,s,amdq,apdq,maux)
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx, maux

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc, maux)
      double precision  auxr(1-mbc:maxm+mbc, maux)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)


c     # Assume at most a 64x64 grid
      integer maxm2
      parameter (maxm2 = 520)

      double precision u2v2(-1:maxm2),
     &       u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2),
     &       g1a2(-1:maxm2),euv(-1:maxm2)


      double precision gamma, gamma1
      common /cparam/  gamma,gamma1

      double precision rhsqrtl, rhsqrtr, pl,pr,rhsq2
      double precision s0
      double precision rho1,rhou1,rhov1,en1,p1,c1,s1
      double precision rhoi,pi,ci,s3,rho2,rhou2,rhov2,en2
      double precision p2,c2,s2,df,cim1,pim1,rhoim1
      double precision sfract
      double precision delta(4)

      double precision a1,a2,a3,a4

      integer mu,mv,i, mw, m

      logical efix
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
      if (-1 .gt. 1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'rpn2 (rpn2eu5.f) : need to increase maxm2'
         stop
      endif
c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv to the orthogonal
c     # momentum:
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
     &        qr(i-1,3)**2)/qr(i-1,1))
         pr = gamma1*(ql(i,4) - 0.5d0*(ql(i,2)**2 +
     &        ql(i,3)**2)/ql(i,1))
         rhsq2 = rhsqrtl + rhsqrtr
         u(i) = (qr(i-1,mu)/rhsqrtl + ql(i,mu)/rhsqrtr) / rhsq2
         v(i) = (qr(i-1,mv)/rhsqrtl + ql(i,mv)/rhsqrtr) / rhsq2
         enth(i) = (((qr(i-1,4)+pl)/rhsqrtl
     &             + (ql(i,4)+pr)/rhsqrtr)) / rhsq2
         u2v2(i) = u(i)**2 + v(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2(i))
         a(i) = dsqrt(a2)
         g1a2(i) = gamma1 / a2
         euv(i) = enth(i) - u2v2(i)
      enddo


c     # now split the jump in q at each interface into waves

c     # find a1 thru a4, the coefficients of the 4 eigenvectors:
      do i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,mu) - qr(i-1,mu)
         delta(3) = ql(i,mv) - qr(i-1,mv)
         delta(4) = ql(i,4) - qr(i-1,4)
         a3 = g1a2(i) * (euv(i)*delta(1)
     &      + u(i)*delta(2) + v(i)*delta(3) - delta(4))
         a2 = delta(3) - v(i)*delta(1)
         a4 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a3) / (2.d0*a(i))
         a1 = delta(1) - a3 - a4
c
c        # Compute the waves.
c
c        # acoustic:
         wave(i,1,1) = a1
         wave(i,mu,1) = a1*(u(i)-a(i))
         wave(i,mv,1) = a1*v(i)
         wave(i,4,1) = a1*(enth(i) - u(i)*a(i))
         wave(i,5,1) = 0.d0
         s(i,1) = u(i)-a(i)
c
c        # shear:
         wave(i,1,2) = 0.d0
         wave(i,mu,2) = 0.d0
         wave(i,mv,2) = a2
         wave(i,4,2) = a2*v(i)
         wave(i,5,2) = 0.d0
         s(i,2) = u(i)
c
c        # entropy:
         wave(i,1,3) = a3
         wave(i,mu,3) = a3*u(i)
         wave(i,mv,3) = a3*v(i)
         wave(i,4,3) = a3*0.5d0*u2v2(i)
         wave(i,5,3) = 0.d0
         s(i,3) = u(i)
c
c        # acoustic:
         wave(i,1,4) = a4
         wave(i,mu,4) = a4*(u(i)+a(i))
         wave(i,mv,4) = a4*v(i)
         wave(i,4,4) = a4*(enth(i)+u(i)*a(i))
         wave(i,5,4) = 0.d0
         s(i,4) = u(i)+a(i)
c
c        # Another wave added for tracer concentration:
c
c        # tracer:
         wave(i,1,5) = 0.d0
         wave(i,mu,5) = 0.d0
         wave(i,mv,5) = 0.d0
         wave(i,4,5) = 0.d0
         wave(i,5,5) = ql(i,5) - qr(i-1,5)
         s(i,5) = u(i)
      end do
c
c
c
c    # compute flux differences amdq and apdq.
c    ---------------------------------------
c
      if (efix) go to 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do m=1,meqn
         do i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw=1,mwaves
               if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
                 else
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
            end do
        end do
       end do
      go to 900
c
c-----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of s*wave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do i = 2-mbc, mx+mbc
c
c        # check 1-wave:
c        ---------------
c
         rhoim1 = qr(i-1,1)
         pim1 = gamma1*(qr(i-1,4) - 0.5d0*(qr(i-1,mu)**2
     &           + qr(i-1,mv)**2) / rhoim1)
         cim1 = dsqrt(gamma*pim1/rhoim1)
         s0 = qr(i-1,mu)/rhoim1 - cim1     !# u-c in left state (cell i-1)

c        # check for fully supersonic case:
         if (s0.ge.0.d0 .and. s(i,1).gt.0.d0)  then
c            # everything is right-going
             do m=1,meqn
                amdq(i,m) = 0.d0
             end do
             go to 200
         endif
c
         rho1 = qr(i-1,1) + wave(i,1,1)
         rhou1 = qr(i-1,mu) + wave(i,mu,1)
         rhov1 = qr(i-1,mv) + wave(i,mv,1)
         en1 = qr(i-1,4) + wave(i,4,1)
         p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2)/rho1)
         c1 = dsqrt(gamma*p1/rho1)
         s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0 .lt. 0.d0 .and. s1 .gt. 0.d0) then
c            # transonic rarefaction in the 1-wave
             sfract = s0 * (s1-s(i,1)) / (s1-s0)
           else if (s(i,1) .lt. 0.d0) then
c            # 1-wave is leftgoing
             sfract = s(i,1)
           else
c            # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
           endif
         do m=1,meqn
            amdq(i,m) = sfract*wave(i,m,1)
         end do
c
c        # check contact discontinuity:
c        ------------------------------
c
         if (s(i,2) .ge. 0.d0) go to 200  !# 2- 3- and 5-waves are rightgoing
         do m=1,meqn
            amdq(i,m) = amdq(i,m) + s(i,2)*wave(i,m,2)
            amdq(i,m) = amdq(i,m) + s(i,3)*wave(i,m,3)
            amdq(i,m) = amdq(i,m) + s(i,5)*wave(i,m,5)
         end do
c
c        # check 4-wave:
c        ---------------
c
         rhoi = ql(i,1)
         pi = gamma1*(ql(i,4) - 0.5d0*(ql(i,mu)**2
     &           + ql(i,mv)**2) / rhoi)
         ci = dsqrt(gamma*pi/rhoi)
         s3 = ql(i,mu)/rhoi + ci     !# u+c in right state  (cell i)
c
         rho2 = ql(i,1) - wave(i,1,4)
         rhou2 = ql(i,mu) - wave(i,mu,4)
         rhov2 = ql(i,mv) - wave(i,mv,4)
         en2 = ql(i,4) - wave(i,4,4)
         p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2)/rho2)
         c2 = dsqrt(gamma*p2/rho2)
         s2 = rhou2/rho2 + c2   !# u+c to left of 4-wave
         if (s2 .lt. 0.d0 .and. s3 .gt. 0.d0) then
c            # transonic rarefaction in the 4-wave
             sfract = s2 * (s3-s(i,4)) / (s3-s2)
           else if (s(i,4) .lt. 0.d0) then
c            # 4-wave is leftgoing
             sfract = s(i,4)
           else
c            # 4-wave is rightgoing
             go to 200
           endif
c
         do m=1,meqn
            amdq(i,m) = amdq(i,m) + sfract*wave(i,m,4)
         end do
  200    continue
        end do
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do m=1,meqn
         do i = 2-mbc, mx+mbc
            df = 0.d0
            do mw=1,mwaves
               df = df + s(i,mw)*wave(i,m,mw)
            end do
            apdq(i,m) = df - amdq(i,m)
         end do
       end do
c
  900 continue
      return
      end
