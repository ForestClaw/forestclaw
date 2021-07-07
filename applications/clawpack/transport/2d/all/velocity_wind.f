      double precision function psi(xd,yd,zd,t)
      implicit none

      double precision xd, yd, zd, t

      double precision kappa,tfinal
      common /windparms/ kappa,tfinal

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision l, th, lp

c      call get_wind_parms(kappa,tfinal)

      call map2polar(xd,yd,zd,l,th)

      lp = l - 2*pi*t/tfinal

      psi = kappa*sin(lp)**2*cos(th)**2*cos(pi*t/Tfinal) -
     &      2*pi*sin(th)/Tfinal

c     # Sign difference from Benchmark problem
      psi = -psi


      end

      subroutine get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3),t) -
     &      psi(xd2(1),xd2(2),xd2(3),t))/ds

      end


      subroutine map2polar(x,y,z,lambda,th)
      implicit none

      double precision x,y,z,th,lambda
      double precision r

      double precision pi, pi2
      common /compi/ pi, pi2

      r = sqrt(x**2 + y**2 + z**2)
      th = asin(z/r)

      lambda = atan2(y,x)
      if (lambda .lt. 0) then
         lambda = lambda + 2*pi
      endif

      end
