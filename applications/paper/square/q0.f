c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision x0, y0, z0, q0, r0, r
      double precision xpp, ypp
      double precision Hsmooth



c     # Sphere centered at (0.5,0.5,0) on swirl
      if (initchoice .le. 1) then
          x0 = 0.5
          y0 = 0.5
          z0 = 0
          r0 = 0.2   !! radius of sphere

          xpp = modulo(xp,1.d0)
          ypp = modulo(yp,1.d0)

          r = sqrt((xpp - x0)**2 + (ypp-y0)**2)
          q0 = Hsmooth(r + r0) - Hsmooth(r - r0)
      elseif (initchoice .eq. 2) then
          q0 = 1.d0
      endif

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end



