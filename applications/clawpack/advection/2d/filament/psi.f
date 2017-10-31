      double precision function psi(x,y,z)
      implicit none

      double precision x,y,z,pi,r
      common /compi/ pi

c      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)

c     # Rigid body rotation
c      psi = r**2

c     # Filament formation (negative for clockwise rotation)
      psi = (4.d0/3.d0)*r**3

c      psi = x + y

      return
      end

      subroutine get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3)) -
     &      psi(xd2(1),xd2(2),xd2(3)))/ds

      end
