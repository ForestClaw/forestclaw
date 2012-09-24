      subroutine setprob()
      implicit none

      integer ichoice
      double precision kappa,tfinal
      integer n, m

      double precision rp2, th, lambda, wc(3,10), rps, a, w(3)

      double precision r, hmax,b,c, rot_angle(2), scale
      integer meqn, maptype

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # read in data to see whether we should run with fixed rotation or not.
      open(7,file='setprob.data')
      read(7,*) ichoice
      read(7,*) kappa
      read(7,*) Tfinal
      close(7)


      call set_init_choice(ichoice)

c     # -------------------------------------------------
c     # Mapping routines
c     # -------------------------------------------------

c     # Set mapping scaling
      scale = 1
      rot_angle(1) = pi/3.0
      rot_angle(2) = -pi/2.d0

      call setup_mappedgrid(rot_angle,scale)


c     # Tfinal is needed here because it is used in determining the velocity
c     # field.  Again, an options file would really help here because then I
c     # wouldn't have to pass it in.
      call set_wind_parms(kappa,tfinal)

c     # -------------------------------------------------
c     # Setup location of cosine or Gaussians in array wc.
c     # Array 'wc' will be set in common block 'prob.i'
c     # -------------------------------------------------
c     # Locations of cosine bell or Gaussian

      th = pi/12.d0
c      th = 0
      lambda = 5.d0*pi/6.d0
      wc(1,1) = cos(th)*cos(lambda)
      wc(2,1) = cos(th)*sin(lambda)
      wc(3,1) = sin(th)

      th = -pi/12.d0
c      th = 0
      lambda = 7.d0*pi/6.d0
      wc(1,2) = cos(th)*cos(lambda)
      wc(2,2) = cos(th)*sin(lambda)
      wc(3,2) = sin(th)

      n = 2

      call set_initial_gaussian_locations(wc,n)

      if (ichoice .eq. 1) then
         a = 5
         hmax = 1d0
         call set_initial_gaussian_parms(a,hmax)
      elseif (ichoice .eq. 2 .or. ichoice .eq. 3) then
         r = 0.5d0
         hmax = 1.d0
         b = 0.1d0
         c = 0.9d0
         call set_initial_cosbell_parms(r,hmax,b,c)
      elseif (ichoice .eq. 4) then
         r = 0.5d0
         hmax = 1.d0
         b = 0.1d0
         c = 1d0
         call set_initial_sdisk_parms(r,hmax,b,c)
      endif

      end
