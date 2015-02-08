      subroutine setprob
      implicit double precision (a-h,o-z)
      common /param/  g    !# gravitational parameter
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hin,hout

c     # These should be read in as options.
      idisc = 2

      g = 1.d0
      x0 = 0.d0
      y0 = 0.d0
      r0 = 0.5d0
      hin = 2.d0
      hout = 1.d0

      return
      end
