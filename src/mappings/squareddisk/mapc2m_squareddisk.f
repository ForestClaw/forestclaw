c     # -------------------------------------------------------
c     # adapted from the C-function 'fclaw2d_map_c2m_disk'
c     # in fclaw2d_map.c (C. Burstedde)
c     # -------------------------------------------------------

      subroutine mapc2m_squareddisk(xc,yc,xp,yp,zp,alpha)
      implicit none

      double precision xc_in,yc_in,xp,yp,zp
      double precision alpha
c     double precision half_length
c     double precision R2sqrbyR1, R1byR2, R1, R2
      double precision xc1, yc1, xc, yc
      integer blockno, get_block
      logical l1,l2, u1, u2, iscorner
      logical c0, c1, c2, c3
      double precision pi

      common /compi/ pi

      pi = 4.d0*atan(1.0)

      blockno = get_block()

      if (blockno .eq. 2) then
c         half_length = user_double(2);
         xp = (2*xc - 1)*alpha/sqrt(2.d0);
         yp = (2*yc - 1)*alpha/sqrt(2.d0);
      else
c         R2sqrbyR1 = user_double(0)
c         R1byR2 = user_double(1)
c         R1 = user_double(0)
c         R2 = user_double(1)

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
      double precision R2sqrbyR1, R1byR2, R1, R2, alpha
      double precision R, tan_xi, xi_prime

      double precision pi, a, eta1

      common /compi/ pi

c     # Assume for basic mapping, that R1 = 1; R2 = alpha
      eta1 = sin(pi*eta/2)
      eta1 = eta
      R = alpha**(1-eta1)
      tan_xi = tan(0.5d0*pi*(xi - 0.5d0))
      xi_prime = 2.d0*(1.d0 - eta1)*(xi-0.5d0)+eta1*tan_xi

      y = R/sqrt(1.d0 + eta*tan_xi**2 + (1.d0-eta1))
      x = y*xi_prime


      end
