      double precision function get_vel_psi(xd1,xd2,ds,vn)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, psi

      vn = (psi(xd1(1),xd1(2)) -
     &      psi(xd2(1),xd2(2)))/ds

      end

      double precision function psi(xc,yc)
      implicit none

      double precision xc,yc

      double precision xp, yp, zp
      double precision mountain_height, zm

      integer*8 cont, get_context
      integer blockno, fc2d_clawpack46_get_block

      cont = get_context()
      blockno = fc2d_clawpack46_get_block()

      call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
      zm = mountain_height(xp)

      psi = 1.d0/(2000-zm)*(yp - zm)

      end
