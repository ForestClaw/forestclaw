      subroutine radialdam_setprob(g,x0,y0,r0,
     &      hin,hout)
      implicit none

      double precision g,x0,y0,r0,hin,hout
      double precision g_com,x0_com,y0_com
      double precision r0_com,hin_com,hout_com
      double precision alf,beta
      integer idisc

      common /param/  g_com    !# gravitational parameter
      common/cdisc/ x0_com,y0_com,alf,beta,r0_com,idisc
      common /comic/ hin_com,hout_com

c     # These should be read in as options.
      idisc = 2

      g_com = g
      x0_com = x0
      y0_com = y0
      r0_com = r0
      hin_com = hin
      hout_com = hout

      return
      end
