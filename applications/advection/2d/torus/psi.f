      double precision function psi(x,y,z,t)
      implicit none

      double precision x, y, z, t
      double precision pi, r
      logical iscart

      common /compi/ pi


      if (iscart()) then
         psi = y - x
      else
         r = sqrt(x*x + y*y)
         psi = r**2
      endif

      end
