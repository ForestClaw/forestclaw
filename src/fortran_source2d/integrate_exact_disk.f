      subroutine integrate_exact(mx,my,xlower, ylower, dx,dy,
     &      f,deg,qp,qd,isavg)

      external f
      integer mx,my,deg
      double precision xlower, ylower,dx,dy
      double precision qp(0:mx+1,0:my+1)
      double precision qd(mx+1,my+1)
      logical isavg

      call integrate_exact_p(mx,my,xlower,ylower,dx,dy,f,deg,qp,isavg)
      call integrate_exact_d(mx,my,xlower,ylower,dx,dy,f,deg,qd,isavg)

      end


c     # ------------------------------------------------------
c     # area using high order quadrature rule
c     # ------------------------------------------------------

      subroutine integrate_exact_p(mx,my,xlower,ylower,dx,dy,
     &      f,deg,ap,isavg)
      implicit none

      external f,area_element_exact
      integer mx,my, deg
      double precision xlower,ylower
      logical isavg

      double precision ap(0:mx+1,0:my+1)

      integer i,j,ii,jj, k
      double precision xlow, ylow, a,b,c,d,dx,dy, area
      double precision quad_square,quad_triangle
      logical isdiag_p

c     # Primary cells - easy, except for diagonal cells
      do i = 1,mx
         do j = 1,my
            xlow = xlower + (i-1)*dx
            ylow = ylower + (j-1)*dy
            a = xlow
            b = xlow + dx
            c = ylow
            d = ylow + dy
            if (.not. isdiag_p(mx,my,i,j)) then
c              # map is smooth over this cell
               ap(i,j) = quad_square(f,area_element_exact,a,b,c,d,
     &               deg,isavg,area)
            else
c              # map has singularities along some diagonal
               ap(i,j) = quad_triangle(f,area_element_exact,a,b,c,d,
     &               deg,isavg)
            endif
         enddo
      enddo

      end

      subroutine integrate_exact_d(mx,my,xlower,ylower,dx,dy,
     &      f,deg,ad,isavg)
      implicit none

      external f,area_element_exact
      integer mx,my, deg
      logical isavg
      double precision xlower,ylower

      double precision ad(mx+1,my+1)

      integer i,j
      double precision xlow, ylow, a,b,c,d,dx,dy,dx2,dy2, area
      double precision quad_square,quad_triangle
      logical isdiag_d

      double precision xupper,yupper

      dx2 = dx/2.d0
      dy2 = dy/2.d0

      xupper = xlower + mx*dx
      yupper = ylower + my*dy

c     # Dual cells
      do i = 1,mx+1
         do j = 1,my+1
            xlow = xlower + (i-1.5)*dx
            ylow = ylower + (j-1.5)*dy
            call fix_edges(xlow,ylow,xlower,ylower,xupper,yupper)
            a = xlow
            c = ylow
            if (i .eq. 1 .or. i .eq. mx+1) then
               b = xlow + dx2
            else
               b = xlow + dx
            endif
            if (j .eq. 1 .or. j .eq. my+1) then
               d = ylow + dy2
            else
               d = ylow + dy
            endif
            if (.not. isdiag_d(mx,my,i,j)) then
c              # map is smooth over this cell
               ad(i,j) = quad_square(f,area_element_exact,a,b,c,d,
     &               deg,isavg,area)
            else
c              # map has singularities along some diagonal
               ad(i,j) = quad_triangle(f,area_element_exact,a,b,c,d,
     &               deg,isavg)
            endif
         enddo
      enddo



      end

      double precision function area_element_exact(xc,yc)
      implicit none

      double precision xc,yc,txi(3),teta(3)
      double precision norm_cross

      call basis_vectors(xc,yc,txi,teta)

      area_element_exact = norm_cross(txi,teta)

      end


      double precision function quad_triangle(f,ae,a,b,c,d,deg,isavg)
      implicit none

      external f,ae
      double precision a,b,c,d
      integer deg
      logical isavg
      double precision x(8), y(8), wx(8),wy(8)
      integer i,j,k
      double precision a1,b1,c1,d1
      double precision xc,yc,xm,ym,xv(3),yv(3),sv(2),s
      double precision fint1,fint2,area1,area2, h

      double precision f,ae

      xm = (a + b)/2.d0
      ym = (c + d)/2.d0

c     # Bottom and top triangles

      yv(1) = c
      yv(2)= ym
      yv(3) = d
      sv(1) = 1
      sv(2) = -1
      fint2 = 0
      area2 = 0
      do k = 1,2
         call gltable(y,wy,deg,yv(k),yv(k+1))
         s = sv(k)
         do j = 1,deg
            yc = y(j)
c           # (x-xm = (y-ym) at right edge
c           # (x-xm = -(y-ym) at left edge
            a1 =  (yc-ym) + xm
            b1 = -(yc-ym) + xm
            call gltable(x,wx,deg,a1,b1)
            area1 = 0.d0
            fint1 = 0.d0
            do i = 1,deg
               xc = x(i)
               h = s*wx(i)*ae(xc,yc)
               area1 = area1 + h
               fint1 = fint1 + f(xc,yc)*h
            enddo
            area2 = area2 + wy(j)*area1
            fint2 = fint2 + wy(j)*fint1
         enddo
      enddo

      xv(1) = a
      xv(2)= xm
      xv(3) = b
      sv(1) = 1
      sv(2) = -1
      do k = 1,2
         call gltable(x,wx,deg,xv(k),xv(k+1))
         s = sv(k)
         do i = 1,deg
            xc = x(i)
c           # (y-ym) = (x-xm) at right edge
c           # (y-ym) = -(x-xm) at left edge
            c1 =  (xc-xm) + ym
            d1 = -(xc-xm) + ym
            call gltable(y,wy,deg,c1,d1)
            area1 = 0.d0
            fint1 = 0.d0
            do j = 1,deg
               yc = y(j)
               h = s*wy(j)*ae(xc,yc)
               area1 = area1 + h
               fint1 = fint1 + f(xc,yc)*h
            enddo
            fint2 = fint2 + wx(i)*fint1
            area2 = area2 + wx(i)*area1
         enddo
      enddo

      quad_triangle = fint2
      if (isavg) then
         quad_triangle = quad_triangle/area2
      endif

      end


      double precision function quad_general(f,ae,quad,deg,isavg,area)
      implicit none

      external f,ae
      double precision xv(5),yv(5), area
      integer deg
      logical isavg
      double precision a00(2),a01(2),a10(2),a11(2)
      double precision x(8),wx(8),y(8),wy(8)
      double precision xc,yc,txi(3), teta(3)
      double precision norm_cross, quad(0:1,0:1,2)
      double precision f,ae
      double precision h,nc,fint1,fint2,area1,area2
      integer i,j, k

c     # This is NOT a trilinear approximation on the surface, but
c     # rather just a way to integrate over a general quadrilateral
c     # in computational space, so it is in fact still exact for the
c     # surface
      do k = 1,2
         a00(k) = quad(0,0,k)
         a01(k) = quad(1,0,k) - quad(0,0,k)
         a10(k) = quad(0,1,k) - quad(0,0,k)
         a11(k) = quad(1,1,k) - quad(1,0,k) - quad(0,1,k) + quad(0,0,k)
      enddo

      call gltable(x,wx,deg,0.d0,1.d0)
      call gltable(y,wy,deg,0.d0,1.d0)

      area2 = 0
      fint2 = 0
      do i = 1,deg
         fint1 = 0
         area1 = 0
         do j = 1,deg
            xc = a00(1) + a01(1)*x(i) + a10(1)*y(j) + a11(1)*x(i)*y(j)
            yc = a00(2) + a01(2)*x(i) + a10(2)*y(j) + a11(2)*x(i)*y(j)
            do k = 1,2
               txi(k)  = a01(k) + a11(k)*y(j)
               teta(k) = a10(k) + a11(k)*x(i)
            enddo
            txi(3) = 0
            teta(3) = 0
            nc = norm_cross(txi,teta)
            h = wy(j)*nc*ae(xc,yc)
            area1 = area1 + h
            fint1 = fint1 + f(xc,yc)*h
         enddo
         fint2 = fint2 + wx(i)*fint1
         area2 = area2 + wx(i)*area1
      enddo

      quad_general = fint2
      if (isavg) then
         quad_general = quad_general/area2
      endif
      area = area2

      end
