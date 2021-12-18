      double precision function psi(xp,yp,zp,t)
      implicit none

      double precision xp,yp,zp,t
      double precision revs_per_s

      double precision pi, pi2
      common /compi/ pi, pi2

      revs_per_s = 0.5d0

      psi = pi2*revs_per_s*zp

      psi = -psi

      end

      subroutine get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi,t

      vn = (psi(xd1(1),xd1(2),xd1(3),t) -
     &      psi(xd2(1),xd2(2),xd2(3),t))/ds

      end
