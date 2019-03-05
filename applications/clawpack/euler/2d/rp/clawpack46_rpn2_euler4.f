      subroutine clawpack46_rpn2_euler4(ixy,maxm,meqn,mwaves,
     &      mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension  apdq(1-mbc:maxm+mbc, meqn)
      dimension  amdq(1-mbc:maxm+mbc, meqn)
c
c     local arrays -- common block comroe is passed to rpt2eu
c     ------------
      parameter (maxm2 = 602)  !# assumes at most 200x200 grid with mbc=2
      dimension delta(4)
      logical efix
      common /cparam/  gamma,gamma1
      common /comroe/ u2v2(-1:maxm2),
     &       u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2),
     &       g1a2(-1:maxm2),euv(-1:maxm2)
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpn2'
         stop
      endif

      if (meqn .ne. 4) then
         write(6,*) 'clawpack46_rpn2_euler4.f : meqn must equal to 4'
         stop
      endif

c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv to the orthogonal
c     # momentum:
c
      if (ixy.eq.1) then
          mu = 2
          mv = 3
        else
          mu = 3
          mv = 2
        endif
c
c     # note that notation for u and v reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direciton and v in the orthogonal direcion, but with the above
c     # definitions of mu and mv the routine also works with ixy=2
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt2eu to do the transverse wave splitting.
c

      do 10 i = 2-mbc, mx+mbc
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
   10    continue
c
c
c     # now split the jump in q at each interface into waves
c
c     # find a1 thru a4, the coefficients of the 4 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
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
c        # Note that the 2-wave and 3-wave travel at the same speed and
c        # are lumped together in wave(.,.,2).  The 4-wave is then stored in
c        # wave(.,.,3).
c
c        # acoustic:
         wave(i,1,1) = a1
         wave(i,mu,1) = a1*(u(i)-a(i))
         wave(i,mv,1) = a1*v(i)
         wave(i,4,1) = a1*(enth(i) - u(i)*a(i))
         s(i,1) = u(i)-a(i)
c
c        # shear:
         wave(i,1,2) = 0.d0
         wave(i,mu,2) = 0.d0
         wave(i,mv,2) = a2
         wave(i,4,2) = a2*v(i)
         s(i,2) = u(i)
c
c        # entropy:
         wave(i,1,3) = a3
         wave(i,mu,3) = a3*u(i)
         wave(i,mv,3) = a3*v(i)
         wave(i,4,3) = a3*0.5d0*u2v2(i)
         s(i,3) = u(i)
c
c        # acoustic:
         wave(i,1,4) = a4
         wave(i,mu,4) = a4*(u(i)+a(i))
         wave(i,mv,4) = a4*v(i)
         wave(i,4,4) = a4*(enth(i)+u(i)*a(i))
         s(i,4) = u(i)+a(i)
c
   20    continue
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
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do 90 mw=1,mwaves
               if (s(i,mw) .lt. 0.d0) then
                   amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
                 else
                   apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                 endif
   90          continue
  100       continue
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
      do 200 i = 2-mbc, mx+mbc
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
             do 60 m=1,meqn
                amdq(i,m) = 0.d0
   60           continue
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
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
c            # transonic rarefaction in the 1-wave
             sfract = s0 * (s1-s(i,1)) / (s1-s0)
           else if (s(i,1) .lt. 0.d0) then
c            # 1-wave is leftgoing
             sfract = s(i,1)
           else
c            # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
           endif
         do 120 m=1,meqn
            amdq(i,m) = sfract*wave(i,m,1)
  120       continue
c
c        # check contact discontinuity:
c        ------------------------------
c
         if (s(i,2) .ge. 0.d0) go to 200  !# 2- and 3-waves are rightgoing
         do 140 m=1,meqn
            amdq(i,m) = amdq(i,m) + s(i,2)*wave(i,m,2)
            amdq(i,m) = amdq(i,m) + s(i,3)*wave(i,m,3)
  140       continue
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
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0) then
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
         do 160 m=1,meqn
            amdq(i,m) = amdq(i,m) + sfract*wave(i,m,4)
  160       continue
  200    continue
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
               df = df + s(i,mw)*wave(i,m,mw)
  210          continue
            apdq(i,m) = df - amdq(i,m)
  220       continue
c
  900 continue
      return
      end
