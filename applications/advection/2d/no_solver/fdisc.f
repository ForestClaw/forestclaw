      double precision function fdisc(x,y)
      implicit none

      double precision x,y, r
      double precision xp, yp, rp, th, x0
      double precision y0, pi
      integer m, ichoice

      common /com_init/ ichoice

      pi = 4.d0*atan(1.d0)

      xp = x
      yp = y

      fdisc = -1
      rp = sqrt((xp)**2 + (yp)**2)

c     # Large inner circle
      if (rp .lt. 0.4d0) then
         fdisc = 1.d0
      endif

      if (ichoice .eq. 1) then
         return
      endif

c     # Smaller surrounding circles
      do m = 1,4
c        # Larger circles
         th = (m-1)*pi/2.d0
         x0 = 0.6*cos(th)
         y0 = 0.6*sin(th)
         rp = sqrt((xp - x0)**2 + (yp-y0)**2)
         if (rp .lt. 0.2) then
            fdisc = 1.d0
         endif

c        # Smaller circles
         th = pi/4.d0 + (m-1)*pi/2.d0
         x0 = 0.55*cos(th)
         y0 = 0.55*sin(th)
         rp = sqrt((xp - x0)**2 + (yp-y0)**2)
         if (rp .lt. 0.15) then
            fdisc = 1.d0
         endif
      enddo
      fdisc = -fdisc


      end
