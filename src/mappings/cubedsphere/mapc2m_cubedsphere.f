      subroutine mapc2m_cubedsphere(blockno,xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno

      if (blockno .eq. 0) then
         call csphere_basic(xc,yc,yp,xp,zp)
         zp = -zp
      elseif (blockno .eq. 1) then
         call csphere_basic(xc,yc,zp,xp,yp)
      elseif (blockno .eq. 2) then
         call csphere_basic(xc,yc,zp,yp,xp)
         xp = -xp
      elseif (blockno .eq. 3) then
         call csphere_basic(xc,yc,xp,yp,zp)
      elseif (blockno .eq. 4) then
         call csphere_basic(xc,yc,xp,zp,yp)
         yp = -yp
      elseif (blockno .eq. 5) then
         call csphere_basic(xc,yc,yp,zp,xp)
      endif

      end

      subroutine csphere_basic(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp
      double precision R, tan_xi, tan_eta
      double precision pi

      common /compi/ pi

      R = 1.d0
      tan_xi = tan(0.5d0*pi*(xc-0.5d0))
      tan_eta = tan(0.5d0*pi*(yc-0.5d0))
      zp = R/sqrt(tan_xi**2 + tan_eta**2 + 1.d0)
      xp = zp*tan_xi
      yp = zp*tan_eta

      end
