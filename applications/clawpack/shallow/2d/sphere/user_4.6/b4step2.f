      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer maxmx, maxmy, mbc, mx, my, meqn, maux
      double precision xlower, ylower, dx,dy,time,dt

      double precision  q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      double precision xlow,ylow
      double precision erx, ery, erz, qn

      integer i,j, rf

c     # We may not even need this if we project out the correct surface
c     # normal from the velocities.
      return

      do i=1,mx
         do j=1,my
            erx = aux(i,j,14)
            ery = aux(i,j,15)
            erz = aux(i,j,16)
            qn = erx*q(i,j,2) + ery*q(i,j,3) + erz*q(i,j,4)

            q(i,j,2) = q(i,j,2) - qn*erx
            q(i,j,3) = q(i,j,3) - qn*ery
            q(i,j,4) = q(i,j,4) - qn*erz
         enddo
      enddo

      return
      end
