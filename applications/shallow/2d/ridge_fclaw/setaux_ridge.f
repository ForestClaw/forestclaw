      subroutine setaux_ridge(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &      maux,aux,xnormals,ynormals,xtangents,ytangents,
     &      surfnormals,level,maxlevel,refratio)
      implicit none

      integer maxmx, maxmy, mbc, mx,my, maux
      integer level, maxlevel, refratio
      double precision xlower, ylower, dx,dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      double precision    xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision    ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision   xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision   ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

c     The aux array has the following elements:
c     1  kappa = ratio of cell area to dxc*dyc
c     2  enx = x-component of normal vector to left edge in tangent plane
c     3  eny = y-component of normal vector to left edge in tangent plane
c     4  enz = z-component of normal vector to left edge in tangent plane
c     5  etx = x-component of tangent vector to left edge in tangent plane
c     6  ety = y-component of tangent vector to left edge in tangent plane
c     7  etz = z-component of tangent vector to left edge in tangent plane
c     8  enx = x-component of normal vector to bottom edge in tangent plane
c     9  eny = y-component of normal vector to bottom edge in tangent plane
c     10  enz = z-component of normal vector to bottom edge in tangent plane
c     11  etx = x-component of tangent vector to bottom edge in tangent plane
c     12  ety = y-component of tangent vector to bottom edge in tangent plane
c     13  etz = z-component of tangent vector to bottom edge in tangent plane
c     14  erx = x-component of unit vector in radial direction at cell ctr
c     15  ery = y-component of unit vector in radial direction at cell ctr
c     16  erz = z-component of unit vector in radial direction at cell ctr
c     17  xc at cell center
c     18  yc at cell center
c     19  bathymetry - averaged over all possible finer cells

      call  assign_comp_centers(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)

      call assign_capacity(maxmx, maxmy, mx,my,mbc,dx,dy,xlower, ylower,
     &      level,maxlevel,refratio,aux,maux)

      call assign_normals(maxmx,maxmy,mx,my,mbc,
     &      xnormals,ynormals,xtangents,ytangents,
     &      surfnormals,aux,maux);

      end



      subroutine assign_comp_centers(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer maxmx, maxmy,mbc,mx,my,maux
      double precision xlower,ylower,dx,dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i,j
      double precision xc,yc

c     17  xc at cell center
c     18  yc at cell center
      do j = 1-mbc,maxmy+mbc
         do i = 1-mbc,maxmx + mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            aux(i,j,17) = xc
            aux(i,j,18) = yc
         enddo
      enddo

      end

      subroutine assign_normals(maxmx, maxmy, mx,my,mbc,
     &      xnormals,ynormals, xtangents,ytangents,
     &      surfnormals,aux,maux)
      implicit none

      integer mx,my,mbc,maxmx,maxmy,maux

      double precision xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision   xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision   ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i,j,m

c     # Compute normals at all interior edges.

c     # Assign x-face normals and tangents
c     2  enx = x-component of normal vector to left edge in tangent plane
c     3  eny = y-component of normal vector to left edge in tangent plane
c     4  enz = z-component of normal vector to left edge in tangent plane
c     5  etx = x-component of tangent vector to left edge in tangent plane
c     6  ety = y-component of tangent vector to left edge in tangent plane
c     7  etz = z-component of tangent vector to left edge in tangent plane
      do j = 1-mbc,my+mbc
         do i = 2-mbc,mx+mbc
            do m = 1,3
               aux(i,j,1+m) = xnormals(i,j,m)
               aux(i,j,4+m) = xtangents(i,j,m)
            enddo
         enddo
      enddo

c     # Assign y-face normals and tangents
c     8  enx = x-component of normal vector to bottom edge in tangent plane
c     9  eny = y-component of normal vector to bottom edge in tangent plane
c     10  enz = z-component of normal vector to bottom edge in tangent plane
c     11  etx = x-component of tangent vector to bottom edge in tangent plane
c     12  ety = y-component of tangent vector to bottom edge in tangent plane
c     13  etz = z-component of tangent vector to bottom edge in tangent plane
      do j = 2-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            do m = 1,3
               aux(i,j,7+m) = ynormals(i,j,m)
               aux(i,j,10+m) = ytangents(i,j,m)
            enddo
         enddo
      enddo

c     # Compute a surface normal. Use the sum of normals.
c     14  erx = x-component of unit vector in radial direction at cell ctr
c     15  ery = y-component of unit vector in radial direction at cell ctr
c     16  erz = z-component of unit vector in radial direction at cell ctr
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            do m = 1,3
               aux(i,j,13+m) = surfnormals(i,j,m)
            enddo
         enddo
      enddo

      end


      subroutine assign_capacity(maxmx, maxmy, mx,my,mbc,dx,dy,
     &      xlower, ylower, level,maxlevel,refratio,aux,maux)
      implicit none

      integer maxmx, maxmy, mx,my,mbc,level, refratio,maxlevel, maux
      double precision dx,dy, xlower, ylower
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii,jj,icell, jcell
      double precision sum_area, xcorner, ycorner
      double precision xp1,yp1,zp1
      double precision quad(0:1,0:1,3)
      double precision get_area_approx
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      double precision dxf, dyf, xdf, ydf, total_area, a
      double precision xef, yef, xe,ye, area1, areab
      integer k, m, ir, jr

      integer rmax, rfactor
      parameter(rmax = 64)
      double precision area_ref(rmax,rmax), avg_bath, areafact

      rfactor = 1
      do ir = level,maxlevel-1
         rfactor = rfactor*refratio
      enddo
      dxf = dx/rfactor
      dyf = dy/rfactor

c     # Set up capacity function and average bathymetry
c
c     1  kappa = ratio of cell area to dxc*dyc
c     19  bathymetry - averaged over all possible finer cells

      areafact = 1.d0/(dx*dy)

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
                  area1 = get_area_approx(quad)
                  area_ref(ir,jr) = area1
                  sum_area = sum_area + area1
               enddo
            enddo
            areab = avg_bath(rmax,rfactor,area_ref,
     &            xe,ye,dx,dy)
            aux(i,j,1) = sum_area*areafact
            aux(i,j,19) = areab/sum_area
         enddo
      enddo

      end



      double precision function avg_bath(rmax,rfactor,
     &      area_ref,xe,ye,dx,dy)
      implicit none

      integer rmax, rfactor
      double precision area_ref(rmax,rmax)
      double precision xe,ye,dx,dy

      integer ii,jj
      double precision xc,yc,bathysum, dxf,dyf
      double precision bmount

      dxf = dx/rfactor
      dyf = dy/rfactor

      bathysum = 0
      do ii = 1,rfactor
         xc = xe + (ii-0.5)*dxf
         do jj = 1,rfactor
            yc = ye + (jj-0.5)*dyf
            bathysum = bathysum + area_ref(ii,jj)*bmount(xc,yc)
         enddo
      enddo

      avg_bath = bathysum

      end
