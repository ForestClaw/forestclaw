      subroutine setup_mesh(mx,my,mbc,xlower,ylower,dx,dy,
     &      xp,yp,zp,xd,yd,zd)
      implicit none

      integer mx,my, mbc
      double precision dx,dy,xlower,ylower

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      integer i,j
      double precision dxf,dyf, xc,yc,xd1, yd1,zd1

c     # We need both cell centered and node locations to
c     # compute the normals at cell edges.
      dxf = dx/2.d0
      dyf = dy/2.d0
      do i = -2*mbc-1,2*(mx+mbc+1)+1
         do j = -2*mbc-1,2*(my+mbc+1)+1
            if (abs(mod(i,2)) .ne. abs(mod(j,2))) then
               cycle
            endif
            xc = xlower + (i-1)*dxf
            yc = ylower + (j-1)*dyf

            call mapc2m(xc,yc ,xd1,yd1,zd1)

            if (abs(mod(i,2)) .eq. 1) then
c              # Physical location of cell vertices
               xd((i-1)/2 + 1, (j-1)/2 + 1) = xd1
               yd((i-1)/2 + 1, (j-1)/2 + 1) = yd1
               zd((i-1)/2 + 1, (j-1)/2 + 1) = zd1
            else
c              # Physical locations of cell centers
               xp(i/2,j/2) = xd1
               yp(i/2,j/2) = yd1
               zp(i/2,j/2) = zd1
            endif
         enddo
      enddo
      end

      subroutine compute_area(mx,my,mbc,dx,dy,xlower, ylower,
     &      area, level,maxlevel,refratio)
      implicit none

      integer mx,my,mbc,level, refratio,maxlevel
      double precision dx,dy, xlower, ylower


      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii,jj, ir, rfactor, icell, jcell
      double precision sum_area, xcorner, ycorner
      double precision xp1,yp1,zp1
      double precision quad(0:1,0:1,3)
      double precision get_area_approx

      double precision dxf, dyf, xdf, ydf, total_area, a
      double precision xef, yef, xe,ye
      integer k, m

      logical debug

      rfactor = 1
      do ir = level,maxlevel-1
         rfactor = rfactor*refratio
      enddo
      dxf = dx/rfactor
      dyf = dy/rfactor

c     # Primary cells.  Note that we don't do anything special
c     # in the diagonal cells - so the areas there are less accurate
c     # than in the rest of the mesh.

c     # Primary cells.  Note that we don't do anything special
c     # in the diagonal cells - so the areas there are less accurate
c     # than in the rest of the mesh.
      do j = -mbc,my+mbc
         do i = -mbc,mx+mbc
            xe = xlower + (i-1)*dx
            ye = ylower + (j-1)*dy
            sum_area = 0.d0
            do ii = 1,rfactor
               do jj = 1,rfactor
                  xef = xe + (ii - 1)*dxf
                  yef = ye + (jj - 1)*dyf
                  do icell = 0,1
                     do jcell = 0,1
                        xcorner = xef + icell*dxf
                        ycorner = yef + jcell*dyf
                        call mapc2m(xcorner,ycorner,xp1,yp1,zp1)
                        quad(icell,jcell,1) = xp1
                        quad(icell,jcell,2) = yp1
                        quad(icell,jcell,3) = zp1
                     enddo
                  enddo
                  sum_area = sum_area + get_area_approx(quad)
               enddo
            enddo
            area(i,j) = sum_area
         enddo
      enddo

c      open(10,file='area.out');
c      do i = 1,mx
c         do j = 1,my
c            write(10,'(2I5,E16.8)') i,j,area(i,j)
c         enddo
c         write(10,*) ' '
c      enddo
c      close(10)
c
c      open(10,file='points.out');
c      do i = 1,mx+1
c         do j = 1,my+1
c            write(10,'(2I5,3E16.8)') i,j,xd(i,j),yd(i,j),zd(i,j)
c         enddo
c      enddo
c      close(10)


      end


c     # This is the area element based on a bilinear approximation to
c     # the surface defined by the corners of the mesh cell.
      double precision function get_area_approx(quad)
      implicit none

      double precision quad(0:1,0:1,3)

      double precision a00(3), a01(3),a10(3), a11(3)
      double precision Txi(3), Teta(3)
      double precision norm_cross, xi,eta

      integer m

c     # Coefficients for bilinear approximation to surface
      do m = 1,3
         a00(m) = quad(0,0,m)
         a01(m) = quad(1,0,m) - quad(0,0,m)
         a10(m) = quad(0,1,m) - quad(0,0,m)
         a11(m) = quad(1,1,m) - quad(1,0,m) -
     &         quad(0,1,m) + quad(0,0,m)
      enddo

      eta = 0.5d0
      xi  = 0.5d0

c     # Mesh square is approximated by
c     #       T(xi,eta) = a00 + a01*xi + a10*eta + a11*xi*eta
c     #
c     # Area is the length of the cross product of T_xi and T_eta.
c     # computed at the center of the mesh cell.
      do m = 1,3
         txi(m)  = (a01(m) + a11(m)*eta)
         teta(m) = (a10(m) + a11(m)*xi)
      enddo

      get_area_approx = norm_cross(txi,teta)

      end

      double precision function norm_cross(u,v)
      implicit none

      double precision u(3), v(3), w(3)

      w(1) =  u(2)*v(3) - u(3)*v(2)
      w(2) = -u(1)*v(3) + u(3)*v(1)
      w(3) =  u(1)*v(2) - u(2)*v(1)

      norm_cross = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))

      end
