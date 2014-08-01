      subroutine setprob()
      implicit none

      double precision kappa,tfinal
      double precision rot_angle(2), scale

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # read in data to see whether we should run with fixed rotation or not.
c     # commented this out since it is overwritten below anyway.
c      open(7,file='setprob.data')
c      read(7,*) kappa
c      read(7,*) Tfinal
c      close(7)

      call set_maptype_cubedsphere()

c     # These are the values that work well with the hemisphere.  Other
c     # values don't work so well, since we don't have any inflow
c     # conditions set.

      kappa = 0
      tfinal = 5.0
      call set_wind_parms(kappa,tfinal);

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
