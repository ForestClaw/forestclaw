c     # -------------------------------------------------------
c     # adapted from the C-function 'fclaw2d_map_c2m_disk'
c     # in fclaw2d_map.c (C. Burstedde)
c     # -------------------------------------------------------

      subroutine mapc2m_squareddisk(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      double precision xc,yc,xp,yp,zp
      double precision alpha
      integer blockno
      double precision xc1, yc1

      double precision pi, pi2
      common /compi/ pi, pi2

      if (blockno .eq. 2) then
         xp = (2*xc - 1)*alpha/sqrt(2.d0);
         yp = (2*yc - 1)*alpha/sqrt(2.d0);
      else
         if (blockno .eq. 0) then
            xc1 = xc
            yc1 = 1.d0-yc
            call squareddisk_help(alpha,xc1,yc1,xp,yp)
            yp = -yp
         elseif (blockno .eq. 1) then
            xc1 = yc
            yc1 = 1.d0-xc
            call squareddisk_help(alpha,xc1,yc1,yp,xp)
            xp = -xp
         elseif (blockno .eq. 3) then
            xc1 = yc
            yc1 = xc
            call squareddisk_help(alpha,xc1,yc1,yp,xp)
         elseif (blockno .eq. 4) then
            xc1 = xc
            yc1 = yc
            call squareddisk_help(alpha,xc1,yc1,xp,yp)
         else
            write(6,*) 'mapc2m_squareddisk : Invalid block number'
            stop
         endif
      endif
      zp = 0

      end

      subroutine squareddisk_help(alpha,xi,eta,x,y)
      implicit none

      double precision xi, eta, x,y
      double precision alpha
      double precision R, tau, xi_prime

      double precision pi, pi2

      common /compi/ pi, pi2

c     # Original R - might be smoother
      R = alpha**(1-eta)

c     # This seems to lead to more uniform cells
c      R = (1-alpha)*eta + alpha

      tau = tan(0.5d0*pi*(xi - 0.5d0))
      xi_prime = 2.d0*(1.d0 - eta)*(xi-0.5d0)+eta*tau

c      y = R/sqrt(1.d0 + eta*tau**2 + (1.d0-eta))

      y = R/sqrt(2.d0 + eta*(tau**2 - 1))
      x = y*xi_prime

      end
