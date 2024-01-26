C     # This returns a unit square map in [-1,1]x[-1,1]
      subroutine fclaw_map_2d_c2m_fivepatch(
     &   blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp,alpha
      double precision xc1, yc1

      if (blockno .eq. 2) then
         xp = (2*xc-1)*alpha
         yp = (2*yc-1)*alpha
      else
         if (blockno .eq. 0) then
            xc1 = xc
            yc1 = 1-yc
            call bilinear_help(alpha,xc1,yc1,xp,yp)
            yp = -yp
         elseif (blockno .eq. 1) then
            xc1 = yc
            yc1 = 1 -xc
            call bilinear_help(alpha,xc1,yc1,yp,xp)
            xp = -xp
         elseif (blockno .eq. 3) then
            xc1 = yc
            yc1 = xc
            call bilinear_help(alpha,xc1,yc1,yp,xp)
         elseif (blockno .eq. 4) then
            xc1 = xc
            yc1 = yc
            call bilinear_help(alpha,xc1,yc1,xp,yp)
         else
            write(6,'(A,I5)') 'mapc2m_fivepatch.f : blockno = ', blockno
            stop
         endif
      endif
      zp = 0
      end

      subroutine bilinear_help(alpha,xi,eta,xp,yp)
      implicit none

      double precision alpha, xi, eta, xp,yp

      double precision corner(4,2)
      double precision a(2),u1(2),v1(2),v2(2), pt(2)
      integer m

      corner(1,1) = -alpha
      corner(4,1) = alpha
      corner(3,1) = 1
      corner(2,1) = -1

      corner(1,2) = alpha
      corner(4,2) = alpha
      corner(3,2) = 1
      corner(2,2) = 1

      do m = 1,2
         a(m)  = corner(1,m)
         u1(m) = corner(4,m) - corner(1,m)
         v1(m) = corner(2,m) - corner(1,m)
         v2(m) = corner(3,m) - corner(4,m)
      enddo

      do m = 1,2
         pt(m) = a(m) + u1(m)*xi + v1(m)*eta + (v2(m)-v1(m))*xi*eta
      enddo

      xp = pt(1)
      yp = pt(2)

      end
