      double precision function filament_psi(x,y)
      implicit none

      double precision x,y,r

      double precision pi, pi2
      common /compi/ pi, pi2

c      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

c     # Rigid body rotation
c      psi = r**2

c     # Filament formation (negative for clockwise rotation)
      filament_psi = (4.d0/3.d0)*r**3

c      psi = x + y

      return
      end

      subroutine filament_get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2)) -
     &      psi(xd2(1),xd2(2)))/ds

      end
