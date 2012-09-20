      double precision function psi(xd,yd,zd,t)
      implicit none

      double precision xd, yd, zd, t
      double precision l, th, lp, Tfinal, kappa, pi
      logical issphere, isdisk

      common /compi/ pi

      call get_wind_parms(kappa,tfinal)

      call map2polar(xd,yd,zd,l,th)

      lp = l - 2*pi*t/Tfinal

c     # Set kappa to zero to get solid body rotation
      if (issphere()) then
         psi = kappa*sin(lp)**2*cos(th)**2*cos(pi*t/Tfinal) -
     &         2*pi*sin(th)/Tfinal
      else if (isdisk()) then
         psi = -2*pi*(xd**2 + yd**2)/Tfinal
      endif

c     # Sign difference from Benchmark problem
      psi = -psi


      end


      subroutine get_vel(xp,yp,zp,vvec,t)
      implicit none

      double precision xp,yp,zp,vvec(3), t
      double precision rsphere, get_rsphere
      double precision kappa,Tfinal, th, l, lp, pi
      double precision Tl(3), Tth(3), u, v
      double precision cth, sth, cl, sl
      integer m

      common /compi/ pi

c      kappa = 2.d0
c      Tfinal = 5.d0

      call get_wind_parms(kappa,Tfinal)

      call map2polar(xp,yp,zp,l,th)

      lp = l - 2*pi*t/Tfinal

      cth = cos(th)
      sth = sin(th)
      cl = cos(l)
      sl = sin(l)

c     # Normalized basis vectors for spherical coordinates
      Tl(1) = -sl
      Tl(2) =  cl
      Tl(3) = 0

      Tth(1) = -sth*cl
      Tth(2) = -sth*sl
      Tth(3) = cth

c     # Velocity of spherical coordinates
      u = kappa*sin(lp)**2*sin(2*th)*cos(pi*t/Tfinal) +
     &      2*pi*cth/Tfinal
      v = kappa*sin(2*lp)*cth*cos(pi*t/Tfinal)

c     # Get Cartesian components of the velocity vector
      do m = 1,3
         vvec(m) = u*Tl(m) + v*Tth(m)
      enddo

      end


      subroutine map2polar(x,y,z,lambda,th)
      implicit none

      double precision x,y,z,th,lambda
      double precision r, pi

      common /compi/ pi

      r = sqrt(x**2 + y**2 + z**2)
      th = asin(z/r)

      lambda = atan2(y,x)
      if (lambda .lt. 0) then
         lambda = lambda + 2*pi
      endif

      end
