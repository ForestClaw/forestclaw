      subroutine setprob()
      implicit none

      double precision kappa,tfinal
      double precision rot_angle(2), scale

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # -------------------------------------------------
c     # Mapping routines
c     # -------------------------------------------------

c     # Set mapping scaling
      scale = 1
      rot_angle(1) = 1.34
      rot_angle(2) = 2.1
      rot_angle(1) = 0
      rot_angle(2) = 1.8
      call setup_mappedgrid(rot_angle,scale)
      end
