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
c      do i = 1-mbc,mx+mbc
c         do j = 1-mbc,my+mbc
c            write(10,'(2I5,E16.8)') i,j,area(i,j)
c         enddo
c         write(10,*) ' '
c      enddo
c      close(10)
c      stop
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



      SUBROUTINE compute_normals(mx,my,mbc,xp,yp,zp,xd,yd,zd,
     &      xnormals,ynormals)
      IMPLICIT NONE

      INTEGER mx,my,mbc

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)

      INTEGER i,j,m, ibc, jbc
      DOUBLE PRECISION taud(3),taup(3),nv(3), sp, xv(3)

c     # Compute normals at all interior edges.


c     # Get x-face normals
      DO j =  -mbc,my+mbc+1
         DO i = 1-mbc,mx+mbc+1

            taud(1) = xp(i,j) - xp(i-1,j)
            taud(2) = yp(i,j) - yp(i-1,j)
            taud(3) = zp(i,j) - zp(i-1,j)

            taup(1) = xd(i,j+1) - xd(i,j)
            taup(2) = yd(i,j+1) - yd(i,j)
            taup(3) = zd(i,j+1) - zd(i,j)

            CALL get_normal(taup,taud,nv,sp)

            DO m = 1,3
               xnormals(i,j,m) = nv(m)
            ENDDO
         ENDDO
      ENDDO


      DO j = 1-mbc,my+mbc+1
         DO i =  -mbc,mx+mbc+1
c           # Now do y-faces
            taud(1) = xp(i,j) - xp(i,j-1)
            taud(2) = yp(i,j) - yp(i,j-1)
            taud(3) = zp(i,j) - zp(i,j-1)

            taup(1) = xd(i+1,j) - xd(i,j)
            taup(2) = yd(i+1,j) - yd(i,j)
            taup(3) = zd(i+1,j) - zd(i,j)

            CALL get_normal(taup,taud,nv,sp)

c           # nv has unit length
            DO m = 1,3
               ynormals(i,j,m) = nv(m)
            ENDDO
         ENDDO
      ENDDO


      END SUBROUTINE compute_normals


      SUBROUTINE compute_tangents(mx,my,mbc,xd,yd,zd,
     &      xtangents,ytangents,edge_lengths)
      IMPLICIT NONE

      INTEGER mx,my,mbc

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      double precision xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision edge_lengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)

      INTEGER i,j,m
      DOUBLE PRECISION taup(3),tlen

c     # Compute normals at all interior edges.

c     # Get x-face normals
      DO j = -mbc,my+mbc+1
         DO i = -mbc,mx+mbc+2

            taup(1) = xd(i,j+1) - xd(i,j)
            taup(2) = yd(i,j+1) - yd(i,j)
            taup(3) = zd(i,j+1) - zd(i,j)
            tlen = sqrt(taup(1)**2 + taup(2)**2 + taup(3)**2)

            DO m = 1,3
               xtangents(i,j,m) = taup(m)
            ENDDO
            edge_lengths(i,j,1) = tlen
         ENDDO
      ENDDO

      DO j = -mbc,my+mbc+2
         DO i = -mbc,mx+mbc+1
c           # Now do y-faces
            taup(1) = xd(i+1,j) - xd(i,j)
            taup(2) = yd(i+1,j) - yd(i,j)
            taup(3) = zd(i+1,j) - zd(i,j)
            tlen = sqrt(taup(1)**2 + taup(2)**2 + taup(3)**2)

c           # nv has unit length
            DO m = 1,3
               ytangents(i,j,m) = taup(m)
            ENDDO
            edge_lengths(i,j,2) = tlen
         ENDDO
      ENDDO

      END SUBROUTINE compute_tangents


c     # Compute an approximate unit normal to cell edge
c     #
c     # Inputs :
c     # taud     Vector between adjoining cell centers ("dual" cell edge)
c     # taup     Vector between edge nodes ("primal" cell edge)
c     #
c     # Ouputs :
c     # nv       Unit vector vector normal to taup, and tangent to the surface
c     # sp       Length of taup
c     #
c     # Idea is to construct normal to edge vector using
c     #
c     #     t^i = a^{i1} t_1 + a^{i2} t_2
c     #
c     # where t_1, t_2 are coordinate basis vectors T_xi, T_eta
c     # and t^i dot t_j = 1 (i == j), 0 (i /= j)
c     #
c     # At left edge, we have
c     #       taud ~ T_xi = t_1  and   taup ~ T_eta = t_2
c     #
c     # At bottom edge, we have
c     #       taup ~ T_xi = t_1  and   taud ~ T_eta = t_2
c     #
c     # Metric a_{ij} = t_i dot t_j
c     # Metric a^{ij} is the inverse of a_{ij}
c     #
c     # Both versions use this idea.  That is, both compute
c     # coefficients c1, c2 so that the unit normal at the edge
c     # is expressed as
c     #
c     #        n = c1*taup + c2*taud
c     #
c     # In version 1, trig. identities are used to express these coefficients
c     # in terms of the angle between taup and taud.  In version 2, entries of
c     # the conjugate tensor are used directly.   Both compute the same vector, but
c     # version 2 is about 20%-30% faster.
c     #
c     # Formulas work for both left and bottom edges.
c     #
c     # Neither version depends on any knowledge of the surface normals.

      subroutine get_normal(taup,taud,nv,sp)
      implicit none

      double precision taup(3),taud(3),nv(3), sp
      integer version

      data version /2/

      if (version .eq .1) then
         call get_normal_ver1(taup,taud,nv)
      elseif (version .eq. 2) then
         call get_normal_ver2(taup,taud,nv,sp)
      endif

      end

      subroutine get_normal_ver1(taup,taud,nv)
      implicit none

      double precision taup(3),taud(3),nv(3)
      double precision sp,sd,dt,ct,st,c1,c2
      double precision tp(3), td(3), nlen

      integer m

      sp = 0
      sd = 0
      do m = 1,3
         sp = sp + taup(m)*taup(m)
         sd = sd + taud(m)*taud(m)
      enddo
      sp = sqrt(sp)
      sd = sqrt(sd)


c     # st = sin(theta)
c     # ct = cos(theta)
c     # theta = angle between taup and taud
      ct = 0
      do m = 1,3
         tp(m) = taup(m)/sp
         td(m) = taud(m)/sd
         ct = ct + tp(m)*td(m)
      enddo
      st = sqrt(1.d0-ct*ct)

      c1 = 1.d0/st
      c2 = -ct/st
      do m = 1,3
         nv(m) = c1*td(m) + c2*tp(m)
         nlen = nlen + nv(m)*nv(m)
      enddo
c      if (nlen .eq. 0) then
c         write(6,*) 'nlen == 0'
c      endif


      end

      subroutine get_normal_ver2(taup,taud,nv,sp)
      implicit none

      double precision taup(3),taud(3),nv(3), sp

      integer m
      double precision nlen,aii,ai2,c1,c2

c     # we leave out any scaling by the determinant of the metric,
c     # since we just normalize the vector anyway.
      aii = 0
      ai2 = 0
      do m = 1,3
         aii = aii + taup(m)*taup(m)
         ai2 = ai2 - taud(m)*taup(m)
      enddo

      nlen = 0
      do m = 1,3
         nv(m) = aii*taud(m) + ai2*taup(m)
         nlen = nlen + nv(m)*nv(m)
      enddo
      nlen = sqrt(nlen)

c      c1 = aii/nlen
c      c2 = ai2/nlen

      do m = 1,3
c         nv(m) = c1*taud(m) + c2*taup(m)
         nv(m) = nv(m)/nlen
      enddo
      sp = sqrt(aii)


      end

      SUBROUTINE compute_surf_normals(mx,my,mbc,
     &      xnormals,ynormals,edge_lengths,curvature,surfnormals)
      IMPLICIT NONE

      INTEGER mx,my,mbc

      double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision edge_lengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
      double precision  surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
      double precision    curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      INTEGER i,j,m
      DOUBLE PRECISION nv(3), sp, nlen

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            do m = 1,3
               nv(m) = 0
            enddo
            nlen = 0
            do m = 1,3
               nv(m) =
     &               edge_lengths(i+1,j,1)*xnormals(i+1,j,m) -
     &               edge_lengths(i,j,1)*xnormals(i,j,m) +
     &               edge_lengths(i,j+1,2)*ynormals(i,j+1,m) -
     &               edge_lengths(i,j,2)*ynormals(i,j,m)
               nlen = nlen + nv(m)*nv(m)
            enddo
            nlen = sqrt(nlen)
            do m = 1,3
               surfnormals(i,j,m) = nv(m)/nlen
            enddo
            curvature(i,j) = nlen
         enddo
      enddo



      end
