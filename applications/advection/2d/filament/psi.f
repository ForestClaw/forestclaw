

c     =================================================
      double precision function psi(x,y)
c     =================================================

c     # stream function

      implicit none

      double precision x,y,pi,r
      common /compi/ pi

c      psi = ((sin(pi*x))**2 * (sin(pi*y))**2) / pi

      r = sqrt((x-1.d0)**2 + (y-1.d0)**2)
      psi = (4.d0/3.d0)*r**3

c      psi = x - y

      return
      end

c     =================================================
c     This is called by compute_velocity_psi or
c     compute_velocity_psi_nomap
c     =================================================
      subroutine get_vel_psi(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3)) -
     &      psi(xd2(1),xd2(2),xd2(3)))/ds

      end
