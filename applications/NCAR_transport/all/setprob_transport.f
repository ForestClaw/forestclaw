      subroutine setprob()
c     # This is here because we dont' otherwise have a setprob in the clawpack
c     solver directory.  Need to fix this!
      end

      subroutine setprob_transport(vflag,ichoice)
      implicit none

      double precision kappa,tfinal
      integer n, m, vflag, ichoice

      double precision rp2, th, lambda, wc(3,10), rps, a, w(3)

      double precision rot_angle(2), scale
      double precision r, hmax,b,c
      integer meqn

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # TODO : Need a better way to read in data parameters in
c     # Fortran for parallel cases

c      open(7,file='setprob.data')
c      read(7,*) vflag
c      read(7,*) ichoice
c      close(7)

      call set_init_choice(ichoice)

c     # -------------------------------------------------
c     # Set velocity flag to cell-centered, or edge centered,
c     # and using stream-function.
c     # -------------------------------------------------
      call set_vflag(vflag)

      kappa = 2.d0
      tfinal = 5.d0
      call set_wind_parms(kappa,tfinal)

c     # -------------------------------------------------
c     # Set up sphere mapping
c     # -------------------------------------------------

c     # -------------------------------------------------
c     # Setup location of cosine or Gaussians in array wc.
c     # Array 'wc' will be set in common block 'prob.i'
c     # -------------------------------------------------
c     # Locations of cosine bell or Gaussian

      th = 0
      lambda = 5.d0*pi/6.d0
      wc(1,1) = cos(th)*cos(lambda)
      wc(2,1) = cos(th)*sin(lambda)
      wc(3,1) = sin(th)

      th = 0
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

      subroutine set_init_choice(ichoice)
      implicit none

      integer ichoice, ichoice_com

      common /cominit/ ichoice_com

      ichoice_com = ichoice


      end

      integer function get_init_choice()
      implicit none

      common /cominit/ ichoice_com
      integer ichoice_com

      get_init_choice = ichoice_com

      end
