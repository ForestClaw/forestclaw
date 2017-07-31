      subroutine setprob()
      implicit none

      double precision kappa,tfinal
      integer n, m, vflag, ichoice

      double precision rp2, th, lambda, wc(3,10), rps, a, w(3)

      double precision rot_angle(2), scale
      double precision r, hmax,b,c
      integer meqn
      double precision g

      double precision pi
      common /compi/ pi

      common /sw/  g

c     # This should be a user parm...
      g = 1.d0

      pi = 4.d0*atan(1.d0)

c     # -------------------------------------------------
c     # Setup location of cosine or Gaussians in array wc.
c     # Array 'wc' will be set in common block 'prob.i'
c     # -------------------------------------------------
c     # Locations of cosine bell or Gaussian

      th = pi/2.d0
      lambda = 5.d0*pi/6.d0
      lambda = 0.0
      wc(1,1) = cos(th)*cos(lambda)
      wc(2,1) = cos(th)*sin(lambda)
      wc(3,1) = sin(th)

      th = 0
      lambda = 7.d0*pi/6.d0
      wc(1,2) = cos(th)*cos(lambda)
      wc(2,2) = cos(th)*sin(lambda)
      wc(3,2) = sin(th)

      n = 1

      call set_initial_gaussian_locations(wc,n)

      a = 20
      hmax = 0.05d0
      call set_initial_gaussian_parms(a,hmax)

      end
