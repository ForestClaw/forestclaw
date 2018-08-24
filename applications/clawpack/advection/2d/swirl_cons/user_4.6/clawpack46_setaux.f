      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit none

      integer mbc, mx, my, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer i, j
      double precision xc,yc,u,v

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            call velocity_field(xc,yc,u,v)
            aux(i,j,1) = u
            aux(i,j,2) = v
 
         enddo
      enddo

      return
      end


      subroutine velocity_field(xc,yc,u,v)
      implicit none

      double precision xc,yc,u,v

      double precision s, pi
      integer example

      common /compi/ pi
      common /comex/ example


      s = sqrt(2.d0)
      if (example .eq. 0) then
c        # Conservative for all solvers (rp=1,2,3,4)              
c         u = cos(2*pi*xc) + 2
c         v = cos(2*pi*yc) + 2
         u = 1
         v = 1
      elseif (example .eq. 1) then
c        # Conservative for all solvers (rp=1,2,3,4)               
         u = s*(cos(pi*xc)**2 + 0.5d0)
         v = s*(sin(pi*yc)**2 + 0.5d0)
      elseif (example .eq. 2) then
         u = s*(cos(pi*xc)**2 - 0.5d0)
         v = s*(sin(pi*yc)**2 - 0.5d0)
      else if (example .eq. 3) then
         u = 1.d0
         v = 0.d0
      else
         write(6,'(A,A)') 'clawpack46_setaux : ',
     &              'No valid example provided'
         stop
      endif



      end
