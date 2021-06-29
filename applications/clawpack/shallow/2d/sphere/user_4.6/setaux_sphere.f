      subroutine setaux_sphere(mx,my,mbc,xlower,ylower,
     &      dx,dy,maux,aux,xnormals,ynormals,
     &      xtangents,ytangents,surfnormals)
      implicit none

      integer mbc, mx,my, maux
      integer level, maxlevel, refratio
      double precision xlower, ylower, dx,dy
      double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j,m
      double precision xc,yc

      include 'fclaw2d_metric_terms.i'

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

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            do m = 1,3
               aux(i,j,1+m) = xnormals(i,j,m)
               aux(i,j,4+m) = xtangents(i,j,m)
            enddo
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            do m = 1,3
               aux(i,j,7+m) = ynormals(i,j,m)
               aux(i,j,10+m) = ytangents(i,j,m)
            enddo
         enddo
      enddo

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            do m = 1,3
               aux(i,j,13+m) = surfnormals(i,j,m)
            enddo
         enddo
      enddo

      end
