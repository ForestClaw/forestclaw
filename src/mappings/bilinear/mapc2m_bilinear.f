      subroutine mapc2m_bilinear(blockno,xc,yc,xp,yp,zp,center)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno
      double precision center(2)

      double precision quad(0:1,0:1,2)
      double precision a00(2), a01(2),a10(2), a11(2)
      double precision ll(2), lr(2), ul(2), ur(2)
      double precision shiftx(0:3), shifty(0:3)
      integer m, i, j

      data ll/0,0/, ul/0, 1/, lr/1, 0/, ur/1,1/
      data shiftx/-1,0,-1,0/, shifty/-1,-1,0,0/

c     # Maps [0,1]x[0,1] to a bilinear map with physical 
c     # corners at [quad(0,0,:), quad(0,1,:), quad(1,0,:), quad(1,1,:)]
      

      do m = 1,2
         quad(0,0,m) = ll(m)
         quad(1,0,m) = lr(m)
         quad(0,1,m) = ul(m)
         quad(1,1,m) = ur(m)
      enddo

      do i = 0,1
         do j = 0,1
            quad(i,j,1) = quad(i,j,1) + shiftx(blockno)
            quad(i,j,2) = quad(i,j,2) + shifty(blockno)
         enddo
      enddo

      do m = 1,2
         if (blockno .eq. 0) then         
            quad(1,1,m) = center(m)
         elseif (blockno .eq. 1) then
            quad(0,1,m) = center(m)
         elseif (blockno .eq. 2) then
            quad(1,0,m) = center(m)
         elseif (blockno .eq. 3) then
            quad(0,0,m) = center(m)
         endif
      enddo

      do m = 1,2
         a00(m) = quad(0,0,m)
         a01(m) = quad(1,0,m) - quad(0,0,m)
         a10(m) = quad(0,1,m) - quad(0,0,m)
         a11(m) = quad(1,1,m) - quad(1,0,m) -
     &         quad(0,1,m) + quad(0,0,m)
      enddo

c     #  T(xi,eta) = a00 + a01*xi + a10*eta + a11*xi*eta

      xp = a00(1) + a01(1)*xc + a10(1)*yc + a11(1)*xc*yc
      yp = a00(2) + a01(2)*xc + a10(2)*yc + a11(2)*xc*yc
      zp = 0


      end
