      subroutine fclaw2d_fort_compute_mesh(mx,my,mbc,
     &      xlower,ylower, dx,dy,blockno,xp,yp,zp,xd,yd,zd)
      implicit none

      integer mx,my, mbc, blockno
      double precision dx,dy,xlower,ylower

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

      integer i,j
      double precision dxf,dyf, xc,yc,xd1, yd1,zd1

      integer*8 map_context_ptr, get_context

      map_context_ptr = get_context()

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


c            call mapc2m(xc,yc ,xd1,yd1,zd1)

            call fclaw2d_map_c2m(map_context_ptr,
     &            blockno,xc,yc,xd1,yd1,zd1)

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

      subroutine fclaw2d_fort_compute_area(mx,my,mbc,dx,dy,
     &      xlower, ylower, blockno,area, quadsize, quadstore,
     &      ghost_only)
      implicit none

      integer mx,my,mbc,blockno, quadsize
      double precision quadstore(0:quadsize,0:quadsize,3)
      double precision dx,dy, xlower, ylower
      integer ghost_only
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer*8 cont, get_context
      logical isaffine

      cont = get_context()
      if (isaffine()) then
c        # We don't need to compute areas all the way to the
c        # finest level.
         call fclaw2d_fort_compute_area_affine(mx,my,mbc,dx,dy,
     &         xlower, ylower, blockno,area, ghost_only)
      else
         call fclaw2d_fort_compute_area_general(mx,my,mbc,dx,dy,
     &         xlower, ylower, blockno,area,quadsize,quadstore,
     &         ghost_only)
      endif

      end

      subroutine fclaw2d_fort_compute_area_general(mx,my,mbc,
     &      dx,dy, xlower, ylower, blockno,area,
     &      quadsize, quadstore,ghost_only)
      implicit none

      integer mx,my,mbc,blockno, quadsize
      double precision dx,dy, xlower, ylower
      double precision quadstore(0:quadsize,0:quadsize,3)
      integer ghost_only


      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii,jj, rfactor, icell, jcell
      double precision sum_area
      double precision xp1,yp1,zp1
      double precision quad(0:1,0:1,3)
      double precision get_area_approx

      double precision dxf, dyf
      double precision xef, yef, xe,ye
      logical is_area_interior

      integer*8 cont, get_context

      cont = get_context()

      rfactor = quadsize
      dxf = dx/rfactor
      dyf = dy/rfactor

c     # Primary cells.  Note that we don't do anything special
c     # in the diagonal cells - so the areas there are less accurate
c     # than in the rest of the mesh.
      do j = -mbc,my+mbc+1
         do i = -mbc,mx+mbc+1
            if (is_area_interior(mx,my,i,j) .and.
     &            ghost_only .eq. 1) then
               cycle
            endif
            xe = xlower + (i-1)*dx
            ye = ylower + (j-1)*dy

            do ii = 0,rfactor
               do jj = 0,rfactor
                  xef = xe + ii*dxf
                  yef = ye + jj*dyf

                  call fclaw2d_map_c2m(cont,
     &                  blockno,xef,yef,xp1,yp1,zp1)

                  quadstore(ii,jj,1) = xp1
                  quadstore(ii,jj,2) = yp1
                  quadstore(ii,jj,3) = zp1
               enddo
            enddo

            sum_area = 0.d0
            do ii = 0,rfactor-1
               do jj = 0,rfactor-1
                  do icell = 0,1
                     do jcell = 0,1
                        quad(icell,jcell,1) =
     &                        quadstore(ii+icell,jj+jcell,1)
                        quad(icell,jcell,2) =
     &                        quadstore(ii+icell,jj+jcell,2)
                        quad(icell,jcell,3) =
     &                        quadstore(ii+icell,jj+jcell,3)
                     enddo
                  enddo
                  sum_area = sum_area + get_area_approx(quad)
               enddo
            enddo
            area(i,j) = sum_area
         enddo
      enddo

      write(6,'(A,E24.16)') 'area = ', area(1,1)

      end


c     # If the mapping is affine, (e.g. Ax + b) then we don't need to sum
c     # the finer level areas.
      subroutine fclaw2d_fort_compute_area_affine(mx,my,mbc,dx,dy,
     &      xlower, ylower, blockno,area,ghost_only)
      implicit none

      integer mx,my,mbc,blockno
      double precision dx,dy, xlower, ylower
      integer ghost_only


      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, icell, jcell
      double precision xcorner, ycorner
      double precision xe,ye, xp1,yp1,zp1
      double precision quad(0:1,0:1,3)
      double precision get_area_approx
      logical is_area_interior

      integer*8 map_context_ptr, get_context

      map_context_ptr = get_context()

      do j = -mbc,my+mbc+1
         do i = -mbc,mx+mbc+1
            if (is_area_interior(mx,my,i,j) .and.
     &            ghost_only .eq. 1) then
               cycle
            endif
            xe = xlower + (i-1)*dx
            ye = ylower + (j-1)*dy
            do icell = 0,1
               do jcell = 0,1
                  xcorner = xe + icell*dx
                  ycorner = ye + jcell*dy
                  call fclaw2d_map_c2m(map_context_ptr,
     &                  blockno,xcorner,ycorner,xp1,yp1,zp1)
                  quad(icell,jcell,1) = xp1
                  quad(icell,jcell,2) = yp1
                  quad(icell,jcell,3) = zp1
               enddo
            enddo
            area(i,j) = get_area_approx(quad)
         enddo
      enddo

      end

      logical function is_area_interior(mx,my,i,j)
      implicit none
      integer mx,my,i,j

      is_area_interior = i .ge. 0 .and. i .le. mx+1 .and.
     &      j .ge. 0 .and. j .le. my+1

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

      subroutine fclaw2d_fort_compute_normals(mx,my,mbc,
     &      xp,yp,zp,xd,yd,zd,xnormals,ynormals)
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

      INTEGER i,j,m
      DOUBLE PRECISION taud(3),taup(3),nv(3), sp

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

            call get_normal(taup,taud,nv,sp)

c           # nv has unit length
            do m = 1,3
               ynormals(i,j,m) = nv(m)
            enddo
         enddo
      enddo

      end subroutine


      subroutine fclaw2d_fort_compute_tangents(mx,my,mbc,
     &      xd,yd,zd, xtangents,ytangents,edge_lengths)
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
               xtangents(i,j,m) = taup(m)/tlen
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
            do m = 1,3
               ytangents(i,j,m) = taup(m)/tlen
            enddo
            edge_lengths(i,j,2) = tlen
         enddo
      enddo

      end subroutine


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
      double precision sp,sd,ct,st,c1,c2
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
      double precision nlen,aii,ai2

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

      subroutine fclaw2d_fort_compute_surf_normals(mx,my,mbc,
     &      xnormals,ynormals,edge_lengths,curvature,
     &      surfnormals,area)
      IMPLICIT NONE

      INTEGER mx,my,mbc

      double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision edge_lengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
      double precision  surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)
      double precision    curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision         area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,m
      double precision nv(3), nlen, av(3), bv(3), kappa

      logical isflat

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            if (isflat()) then
c              # The surface normal can be computed from the cross product of
c              # an xnormal and a ynormal.  The curvature is zero in this
c              # case
               do m = 1,3
c                 # It shouldn't matter which ones we pick.
                  av(m) = xnormals(i,j,m)
                  bv(m) = ynormals(i,j,m)
               enddo
c              # construct cross product
               nv(1) =   av(2)*bv(3) - av(3)*bv(2)
               nv(2) = -(av(1)*bv(3) - av(3)*bv(1))
               nv(3) =   av(1)*bv(2) - av(2)*bv(1)
               nlen = sqrt(nv(1)*nv(1) + nv(2)*nv(2) + nv(3)*nv(3))
               kappa = 0
            else
               nlen = 0
               do m = 1,3
                  nv(m) = 0
               enddo
               do m = 1,3
                  nv(m) =
     &                  edge_lengths(i+1,j,1)*xnormals(i+1,j,m) -
     &                  edge_lengths(i,j,1)*xnormals(i,j,m) +
     &                  edge_lengths(i,j+1,2)*ynormals(i,j+1,m) -
     &                  edge_lengths(i,j,2)*ynormals(i,j,m)
                  nlen = nlen + nv(m)*nv(m)
               enddo
               nlen = sqrt(nlen)
               kappa = 0.5d0*nlen/area(i,j)
            endif
            do m = 1,3
               surfnormals(i,j,m) = nv(m)/nlen
            enddo
            curvature(i,j) = kappa
         enddo
      enddo


      end

c> \ingroup  Averaging
c> Average area of fine grid siblings to parent coarse grid.
      subroutine fclaw2d_fort_average_area(mx,my,mbc,
     &      areacoarse, areafine, igrid)
      implicit none

      integer mx,my,mbc,p4est_refineFactor, refratio, igrid

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj
      double precision sum

c     # This should be refratio*refratio.
      integer i1,j1, r2, m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      double precision kf

      p4est_refineFactor = 2
      refratio = 2

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*mx/p4est_refineFactor

      r2 = refratio*refratio
      do j = 0,my/p4est_refineFactor+1
         do i = 0,mx/p4est_refineFactor +1
            i1 = i+ic_add
            j1 = j+jc_add
            m = 0
            do jj = 1,refratio
               do ii = 1,refratio
                  i2(m) = (i-1)*refratio + ii
                  j2(m) = (j-1)*refratio + jj
                  m = m + 1
               enddo
            enddo
            sum = 0
            do m = 0,r2-1
               kf = areafine(i2(m),j2(m))
               sum = sum + kf
            enddo
            areacoarse(i1,j1) = sum
         enddo
      enddo

c     # Compute area in the ghost cells

      end
