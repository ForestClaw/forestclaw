      subroutine radialdam_setprob(grav_in,x0_in,y0_in,r0_in,
     &      hin_in,hout_in,example)
      implicit none

      double precision grav_in, x0_in, y0_in, r0_in, hin_in,hout_in
      integer example
      double precision grav,x0,y0,r0,hin,hout
      integer idisc

      common /cparam/grav    !# gravitational parameter
      common/cdisc/ x0,y0,r0
      common /comic/ hin,hout
      common /comex/ idisc

c     # These should be read in as options.
      idisc = example

      grav = grav_in

      x0 = x0_in
      y0 = y0_in
      r0 = r0_in
      hin = hin_in
      hout = hout_in

      return
      end
