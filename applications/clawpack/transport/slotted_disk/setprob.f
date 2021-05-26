      subroutine slotted_disk_setprob()
      implicit none

      double precision kappa,tfinal
      integer n, m

      double precision rp2, th, lambda, wc(3,10), rps, a, w(3)

      double precision r, hmax,b,c

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

      open(10,file='setprob.data')
      read(10,*) kappa
      read(10,*) tfinal
      close(10)

      call set_wind_parms(kappa,tfinal)

c     # -------------------------------------------------
c     # Setup location of cosine or Gaussians in array wc.
c     # Array 'wc' will be set in common block 'prob.i'
c     # -------------------------------------------------
c     # Locations of cosine bell or Gaussian

      th = 0
      lambda = pi/6.d0
      wc(1,1) = cos(th)*cos(lambda)
      wc(2,1) = cos(th)*sin(lambda)
      wc(3,1) = sin(th)

      th = 0
      lambda = -pi/6.d0
      wc(1,2) = cos(th)*cos(lambda)
      wc(2,2) = cos(th)*sin(lambda)
      wc(3,2) = sin(th)

c     # Number of cosine bells
      n = 2
      call set_initial_gaussian_locations(wc,n)

      r = 0.5d0
      hmax = 1.d0
      b = 0.1d0
      c = 1.d0
      call set_initial_sdisk_parms(r,hmax,b,c)

      end
