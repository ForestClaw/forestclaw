c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision x0(2), y0(2)
      common /qinit_comm/ x0, y0      

      double precision z0, q0, r0, r
      double precision Hsmooth
      integer k

c     # Sphere centered at (0.5,0.5,0) on swirl
      r0 = 0.35
      z0 = 0
      q0 = 0
      do  k = 1,2
          r = sqrt((xp - x0(k))**2 + (yp-y0(k))**2 + (zp-z0)**2)
          q0 = q0 + Hsmooth(r + r0) - Hsmooth(r - r0)
      end do

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end



