      subroutine mapc2m(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xc1,yc1,xp,yp,zp
      logical isdisk, iscart, ispillowsphere
      logical iscubedsphere

      xc1 = 2*xc - 1.d0
      yc1 = 2*yc - 1.d0

      if (iscart()) then
         call mapc2m_cart(xc1,yc1,xp,yp,zp)
      elseif (isdisk()) then
         call mapc2m_pillowdisk(xc1,yc1,xp,yp,zp)
      elseif (ispillowsphere()) then
         call mapc2m_pillowsphere(xc1,yc1,xp,yp,zp)
      elseif (iscubedsphere()) then
         call mapc2m_cubedsphere(xc,yc,xp,yp,zp)
      endif

      call scale_map(xp,yp,zp)
      call rotate_map(xp,yp,zp)

      end
