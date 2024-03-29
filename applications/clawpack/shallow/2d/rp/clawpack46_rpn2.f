      subroutine clawpack46_rpn2(ixy,maxm,meqn,mwaves,mbc,mx,
     &     ql,qr,auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c     
c     # Roe-solver for the 2D shallow water equations
c     # solve Riemann problems along one slice of data.
c     
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the decomposition of the flux difference
c     #   f(qr(i-1)) - f(ql(i))
c     # into leftgoing and rightgoing parts respectively.
c     # With the Roe solver we have
c     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
c     # where A is the Roe matrix.  An entropy fix can also be incorporated
c     # into the flux differences.
c     
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c     
c     
      implicit none

      integer ixy, maxm, meqn, mwaves, mbc, mx
c     
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc, *)
      double precision  auxr(1-mbc:maxm+mbc, *)

      double precision grav
      common /cparam/ grav      !# gravitational parameter

c     local arrays -- common block comroe is passed to rpt2sh
c     ------------

      integer maxm2
      parameter (maxm2 = 603)   !# assumes at most 600x600 grid with mbc=3

      double precision delta(3), hsqrtl, hsqrtr, hsq2
      double precision a1, a2, a3, s0, s1, s3, him1
      double precision h1, hu1, sfract, hi, s03, h3, hu3
      double precision df
      integer mu, mv, i, m, mw

      logical efix

      double precision, dimension(-2:maxm2) :: u,v,a,h
      common /comroe/ u,v,a,h
c     
      data efix /.true./       !# use entropy fix for transonic rarefactions
c     
      if (-2 .gt. 1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'Check dimensions of local arrays in rpn2'
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
c     # later used in routine rpt2sh to do the transverse wave splitting.
c     
      do i = 2-mbc, mx+mbc
         h(i) = (qr(i-1,1)+ql(i,1))/2.d0
         hsqrtl = dsqrt(qr(i-1,1))
         hsqrtr = dsqrt(ql(i,1))
         hsq2 = hsqrtl + hsqrtr
         u(i) = (qr(i-1,mu)/hsqrtl + ql(i,mu)/hsqrtr) / hsq2
         v(i) = (qr(i-1,mv)/hsqrtl + ql(i,mv)/hsqrtr) / hsq2
         a(i) = dsqrt(grav*h(i))
      enddo
c     
c     
c     # now split the jump in q at each interface into waves
c     
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,mu) - qr(i-1,mu)
         delta(3) = ql(i,mv) - qr(i-1,mv)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))/(2.d0*a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))/(2.d0*a(i))
c     
c     # Compute the waves.
c     
         wave(i,1,1)  = a1
         wave(i,mu,1) = a1*(u(i)-a(i))
         wave(i,mv,1) = a1*v(i)
         s(i,1) = u(i)-a(i)
c     
         wave(i,1,2)  = 0.0d0
         wave(i,mu,2) = 0.0d0
         wave(i,mv,2) = a2
         s(i,2) = u(i)
c     
         wave(i,1,3)  = a3
         wave(i,mu,3) = a3*(u(i)+a(i))
         wave(i,mv,3) = a3*v(i)
         s(i,3) = u(i)+a(i)
      enddo

c     
c     # compute flux differences amdq and apdq.
c     ---------------------------------------
c     
      if (efix) go to 110
c     
c     # no entropy fix
c     ----------------
c     
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c     
      do i = 2-mbc, mx+mbc
         do m = 1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw = 1,mwaves
               amdq(i,m) = amdq(i,m) + 
     &              dmin1(s(i,mw),0.0)*wave(i,m,mw)
               apdq(i,m) = apdq(i,m) + 
     &              dmax1(s(i,mw),0.0)*wave(i,m,mw)
            enddo
         enddo
      enddo
      
      go to 900
c     
c-----------------------------------------------------
c     
 110  continue
c     
c     # With entropy fix
c     ------------------
c     
c     # compute flux differences amdq and apdq.
c     # First compute amdq as sum of s*wave for left going waves.
c     # Incorporate entropy fix by adding a modified fraction of wave
c     # if s should change sign.
c     
      do i = 2-mbc,mx+mbc
c        # check 1-wave
         him1 = qr(i-1,1)
         s0 =  qr(i-1,mu)/him1 - dsqrt(grav*him1)
c        # check for fully supersonic case :
         if (s0 .gt. 0.0d0 .and. s(i,1) .gt. 0.0d0) then
            do m = 1,3
               amdq(i,m) = 0.0d0
            end do
            goto 200
         endif
c     
         h1 = qr(i-1,1) + wave(i,1,1)
         hu1= qr(i-1,mu) + wave(i,mu,1)
         s1 = hu1/h1 - dsqrt(grav*h1) !speed just to right of 1-wave
         if (s0 .lt. 0.0d0 .and. s1 .gt. 0.0d0) then
c           # transonic rarefaction in 1-wave
            sfract = s0*((s1-s(i,1))/(s1-s0))
         else if (s(i,1) .lt. 0.0d0) then
c           # 1-wave is leftgoing
            sfract = s(i,1)
         else
c           # 1-wave is rightgoing
            sfract = 0.0d0
         endif
         do m = 1,3
            amdq(i,m) = sfract*wave(i,m,1)
         end do

c        # check 2-wave
         if (s(i,2) .gt. 0.0d0) then
c           # 2 and 3 waves are right-going
            go to 200
         endif

         do m=1,3
            amdq(i,m) = amdq(i,m) + s(i,2)*wave(i,m,2)
         end do
c     
c     check 3-wave
c     
         hi = ql(i,1)
         s03 = ql(i,mu)/hi + dsqrt(grav*hi)
         h3 = ql(i,1)-wave(i,1,3)
         hu3 = ql(i,mu)-wave(i,mu,3)
         s3 = hu3/h3 + dsqrt(grav*h3)
         if (s3 .lt. 0.0d0 .and. s03 .gt. 0.0d0) then
c           # transonic rarefaction in 3-wave
            sfract = s3*((s03-s(i,3))/(s03-s3))
         else if (s(i,3).lt.0.0d0) then
c           # 3-wave is leftgoing
            sfract = s(i,3)
         else
c           # 3-wave is rightgoing
            goto 200
         endif
         do m = 1,3
            amdq(i,m) = amdq(i,m) + sfract*wave(i,m,3)
         end do
200      continue
      end do
c     
c     compute rightgoing flux differences :
c     
      do m = 1,3
         do i = 2-mbc,mx+mbc
            df = 0.0d0
            do mw = 1,mwaves
               df = df + s(i,mw)*wave(i,m,mw)
            end do
            apdq(i,m)=df - amdq(i,m)
         end do
      end do
c     
c     
 900  continue
      return
      end
