      subroutine setprob()
      implicit none

      double precision th, lambda

      integer refine_criteria
      common /com_refine/ refine_criteria

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
