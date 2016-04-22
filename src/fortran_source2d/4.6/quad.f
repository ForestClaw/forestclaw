c     # ------------------------------------------------------
c     # Misc. functions needed for all quad calculations
c     # ------------------------------------------------------
      double precision function quad_square(f,a,b,c,d,deg,isavg,area)
      implicit none

      external f
      integer deg
      logical isavg
      double precision a,b, c,d,f, area
      double precision x(8),wx(8), xc,yc,qn
      double precision y(8),wy(8), qn2,fc
      double precision area1,fint1,area2,fint2, h
      double precision area_element_exact, ae
      integer i,j
      double precision a_com,b_com,c_com,d_com

c     # These are needed by the routine approximating the surfaces by
c     # trilinear functions
      common /comquad/ a_com,b_com,c_com,d_com

      a_com = a
      b_com = b
      c_com = c
      d_com = d

      call gltable(x,wx,deg,a,b)
      call gltable(y,wy,deg,c,d)

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
            ae = area_element_exact(xc,yc)
            fint1 = fint1 + f(xc,yc)*ae*h
         enddo
         area2 = area2 + wx(i)*area1
         fint2 = fint2 + wx(i)*fint1
      enddo

      quad_square = fint2
      if (isavg) then
         quad_square = quad_square/area2
      endif
      area = area2

      end

c      double precision function area_element(xc,yc)
c      end


c     # Detect mesh cells along the diagonals of the computational
c     # grid. These cells have to be divided into triangles and
c     # integrated over each triangle.
      logical function isdiag(mx,my,i,j)
      implicit none

      integer mx,my,i,j
      logical issphere

      if (.not. issphere()) then
         isdiag = .false.
         if (i - j .eq. 0) then
            isdiag = .true.
         elseif (i + j .eq. my + 1) then
            isdiag = .true.
         endif
      else
         isdiag = .false.
         if (i - j .eq. 0) then
            isdiag = .true.
         elseif (i + j .eq. my + 1) then
            isdiag = .true.
         elseif (i - j .eq. mx/2) then
            isdiag = .true.
         elseif (i + j .eq. mx/2 + my + 1) then
            isdiag = .true.
         endif
      endif

      end

c     # This is used for computing the area, which is the integral of
c     # '1' times the area element.
      double precision function one(xc,yc)
      implicit none

      double precision xc,yc

      one = 1.d0

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
