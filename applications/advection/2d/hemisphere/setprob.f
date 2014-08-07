      subroutine setprob()
      implicit none

      double precision kappa,tfinal
      double precision rot_angle(2), scale
      logical isflat
      integer get_map_value
      logical ishemisphere

      double precision pi
      common /compi/ pi

      pi = 4.d0*atan(1.d0)

c     # read in data to see whether we should run with fixed rotation or not.
      open(7,file='setprob.data')
      read(7,*) kappa
      read(7,*) Tfinal
      close(7)

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
      rot_angle(1) = 0
      rot_angle(2) = 0
      call setup_mappedgrid(rot_angle,scale)
      end
