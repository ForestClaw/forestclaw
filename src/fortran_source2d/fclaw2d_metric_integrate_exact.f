      subroutine fclaw2d_fort_integrate_exact(mx,my,mbc,
     &      dx,dy,xlower, ylower,blockno,area,f,favg,
     &      compute_avg,ghost_only)
      implicit none

      external f
      double precision f
      integer mx,my,deg,mbc,blockno
      double precision xlower, ylower,dx,dy
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision favg(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      integer compute_avg, ghost_only

      deg = 4
      call integrate_exact(mx,my,mbc,xlower,ylower,dx,dy,
     &      blockno,deg,area,f,favg,compute_avg,ghost_only)

      end

      double precision function area_element_exact(xc,yc,blockno)
      implicit none

      double precision xc,yc,txi(3),teta(3)
      integer blockno
      double precision norm_cross
      double precision x(3),y(3),z(3),h
      double precision xc1,yc1

      integer*8 map_context_ptr, get_context

      map_context_ptr = get_context()

      h = 1.0e-8

      xc1 = xc
      yc1 = yc
      call fclaw2d_map_c2m(map_context_ptr,
     &         blockno,xc1,yc1,x(1),y(1),z(1))

      xc1 = xc + h
      yc1 = yc
      call fclaw2d_map_c2m(map_context_ptr,
     &         blockno,xc1,yc1,x(2),y(2),z(2))

      xc1 = xc + h
      yc1 = yc
      call fclaw2d_map_c2m(map_context_ptr,
     &         blockno,xc1,yc1,x(3),y(3),z(3))

      txi(1) = (x(2) - x(1))/h
      txi(2) = (y(2) - y(1))/h
      txi(3) = (z(2) - z(1))/h

      teta(1) = (x(3) - x(1))/h
      teta(2) = (y(3) - y(1))/h
      teta(3) = (z(3) - z(1))/h

      area_element_exact = norm_cross(txi,teta)

      end

c     # This is used for computing the area, which is the integral of
c     # '1' times the area element.
      double precision function one(xc,yc)
      implicit none

      double precision xc,yc

      one = 1.d0

      end


c     # ------------------------------------------------------
c     # area using high order quadrature rule
c     # ------------------------------------------------------

      subroutine integrate_exact(mx,my,mbc,xlower,ylower,dx,dy,
     &      blockno, deg, area,f,favg,compute_avg,ghost_only)
      implicit none

      external one
      external f
      double precision f, one, area_element_exact
      integer mx,my, deg,mbc, blockno
      double precision xlower,ylower
      integer compute_avg, ghost_only
      logical is_area_interior2

      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision favg(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii,jj, k
      double precision xlow, ylow, a,b,c,d,dx,dy, area1
      double precision quad_square

      do i = -mbc,mx+mbc+1
         do j = -mbc,my+mbc+1
            if (is_area_interior2(mx,my,i,j)) then
               cycle
            endif
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            a = xlow
            b = xlow + dx
            c = ylow
            d = ylow + dy
            if (compute_avg .eq. 1) then
               favg(i,j) = quad_square(f,a,b,c,d,deg,
     &               blockno,area1,compute_avg)
               area(i,j) = area1
            else
               area(i,j) = quad_square(one,a,b,c,d,deg,
     &               blockno,area1,compute_avg)
            endif
         enddo
      enddo

      end


      logical function is_area_interior2(mx,my,i,j)
      implicit none
      integer mx,my,i,j

      is_area_interior2 = i .ge. 0 .and. i .le. mx+1 .and.
     &      j .ge. 0 .and. j .le. my+1

      end

c     # ------------------------------------------------------
c     # Misc. functions needed for all quad calculations
c     # ------------------------------------------------------
      double precision function quad_square(f,a,b,c,d,deg,
     &      blockno,area,compute_avg)
      implicit none

      double precision f
      integer deg
      integer compute_avg
      integer blockno
      double precision a,b, c,d, area
      double precision x(8),wx(8), xc,yc,qn
      double precision y(8),wy(8), qn2,fc
      double precision area1,fint1,area2,fint2, h
      double precision area_element_exact, ae
      integer i,j

      call gltable2(x,wx,deg,a,b)
      call gltable2(y,wy,deg,c,d)

      fint2 = 0
      area2 = 0
      do i = 1,deg
         fint1 = 0
         area1 = 0
         xc = x(i)
         do j = 1,deg
            yc = y(j)
            h = wy(j)
            area1 = area1 + h
            ae = area_element_exact(xc,yc,blockno)
            fint1 = fint1 + f(xc,yc)*ae*h
         enddo
         area2 = area2 + wx(i)*area1
         fint2 = fint2 + wx(i)*fint1
      enddo

      quad_square = fint2
      if (compute_avg .eq. 1) then
         quad_square = fint2/area2
      endif
      area = area2

      end


      subroutine gltable2(x,w,d,a,b)
      implicit none

      double precision x(8), w(8), a,b
      integer d, k

      if (d .eq. 1) then
         x(1) = 0
         w(1) = 2.d0
      elseif (d .eq. 2) then
         x(1) = -1.d0/sqrt(3.d0)
         x(2) = -x(1)

         w(1) =  1.d0
         w(2) =  w(1)

      elseif (d .eq. 3) then
         x(1) = -sqrt(3.d0/5.d0)
         x(2) =  0.d0
         x(3) = -x(1)
         w(1) =  5.d0/9.d0;
         w(2) =  8.d0/9.d0;
         w(3) =  w(1)

      elseif (d. eq. 4) then
         x(1) = -0.861136311594053d0
         x(4) = -x(1)
         x(2) = -0.339981043584856d0
         x(3) = -x(2)
         w(1) =  0.347854845137454d0
         w(4) =  w(1)
         w(2) =  0.652145154862546d0
         w(3) =  w(2)

      elseif (d .eq. 5) then
         x(1) = -0.906179845938664d0
         x(5) = -x(1)
         x(2) = -0.538469310105683d0
         x(4) = -x(2)
         x(3) =  0
         w(1) =  0.236926885056189d0
         w(5) =  w(1)
         w(2) =  0.478628670499366d0
         w(4) =  w(2)
         w(3) =  0.568888888888889d0

      elseif (d. eq. 6) then
         x(1) = -0.932469514203152d0
         x(6) = -x(1)
         x(2) = -0.661209386466265d0
         x(5) = -x(2)
         x(3) = -0.238619186083197d0
         x(4) = -x(3)
         w(1) =  0.171324492379170d0
         w(6) =  w(1)
         w(2) =  0.360761573048139d0
         w(5) =  w(2)
         w(3) =  0.467913934572691d0
         w(4) =  w(3)

      elseif (d .eq. 7) then
         x(1) = -0.949107912342759d0
         x(7) = -x(1)
         x(2) = -0.741531185599394d0
         x(6) = -x(2)
         x(3) = -0.405845151377397d0
         x(5) = -x(3)
         x(4) =  0.d0
         w(1) =  0.129484966168870d0
         w(7) =  w(1)
         w(2) =  0.279705391489277d0
         w(6) =  w(2)
         w(3) =  0.381830050505119d0
         w(5) =  w(3)
         w(4) =  0.417959183673469d0

      elseif (d. eq. 8) then
         x(1) = -0.960289856497536d0
         x(8) = -x(1)
         x(2) = -0.796666477413627d0
         x(7) = -x(2)
         x(3) = -0.525532409916329d0
         x(6) = -x(3)
         x(4) = -0.183434642495650d0
         x(5) = -x(4)

         w(1) =  0.101228536290376d0
         w(8) =  w(1)
         w(2) =  0.222381034453374d0
         w(7) =  w(2)
         w(3) =  0.313706645877887d0
         w(6) =  w(3)
         w(4) =  0.362683783378362d0
         w(5) =  w(4)
      else
         write(6,*) 'Degree must be <= 8'
         stop
      endif

c     # Scale nodes and weights to interval [a,b]
      do k = 1,d
         x(k) = (b-a)*x(k)/2.d0 + (a+b)/2.d0
         w(k) = (b-a)*w(k)/2.d0
      enddo

      end
