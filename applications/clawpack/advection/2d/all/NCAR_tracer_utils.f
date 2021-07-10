c     # ------------------------------------------------
c     # Utilities for three NCAR examples : 
c     #     -- gaussian
c     #     -- correlatedcb
c     #     -- slotted_disk
c     # ------------------------------------------------

      subroutine setprob()
      implicit none

      double precision th, lambda

      double precision kappa,tfinal
      common /windparms/ kappa,tfinal

      double precision wcloc(3,2)
      common /location_parms/ wcloc

      double precision a_gauss, hmax_gauss
      common /gaussian_parms/ a_gauss, hmax_gauss

      double precision r_cb, hmax_cb, b_cb, c_cb
      common /cosinebell_parms/ r_cb, hmax_cb, b_cb, c_cb

      double precision r_sd, hmax_sd, b_sd, c_sd
      common /slotteddisk_parms/ r_sd, hmax_sd, b_sd, c_sd

      double precision pi, pi2
      common /compi/ pi, pi2

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

c     # These govern the wind speed; don't change final time here
c     # but rather use configuration tfinal.  This is really just
c     # a parameter and determines when the flow condition returns
c     # to its initial position.
      kappa = 2.0
      tfinal = 5.0

c     # -------------------------------------------------
c     # Locations for cosine bell, Gaussian or slotted disks
c     # -------------------------------------------------

      th = 0
      lambda = pi/6.d0
      wcloc(1,1) = cos(th)*cos(lambda)
      wcloc(2,1) = cos(th)*sin(lambda)
      wcloc(3,1) = sin(th)

      th = 0
      lambda = -pi/6.d0
      wcloc(1,2) = cos(th)*cos(lambda)
      wcloc(2,2) = cos(th)*sin(lambda)
      wcloc(3,2) = sin(th)

c     # Gaussian parameters      
      hmax_gauss = 1
      a_gauss = 5

c     # Cosine bell parameters
      r_cb = 0.5
      hmax_cb = 1
      b_cb = 0.1d0
      c_cb = 0.9d0

c     # slotted disk parameters
      r_sd = 0.5
      hmax_sd = 1
      b_sd = 0.1d0
      c_sd = 0.9d0

      end

c     # ------------------------------------------------
c     # Velocity in terms of a streamfunction
c     # ------------------------------------------------

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


      subroutine get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3),t) -
     &      psi(xd2(1),xd2(2),xd2(3),t))/ds

      end

c     # ------------------------------------------------
c     # Great circle distance on sphere
c     # ------------------------------------------------
      double precision function get_gc_distance(x,y,z,w)
      implicit none

      double precision x,y,z,w(3)
      double precision p(3), pn, wn, ca, th, rsphere
      integer m

      rsphere = 1.0

      p(1) = x
      p(2) = y
      p(3) = z
      pn = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
      if (abs(pn - rsphere) .ge. 1e-3) then
         write(6,*) 'get_gc_distance : Point is not on the sphere; '
         write(6,*) x, y, z, abs(pn-rsphere)
         stop
      endif

      wn = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
      do m = 1,3
         w(m) = w(m)/wn
      enddo
      ca = (w(1)*x + w(2)*y + w(3)*z)/pn
      if (abs(ca) .gt. 1) then
         write(6,*) 'get_gc_distance : abs(ca) > 1; ', ca
         write(6,*) w(1), w(2), w(3), x,y,z
         stop
      endif
      th = acos(ca)
      get_gc_distance = th*rsphere

      end





