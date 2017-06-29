      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, rc,xloc(0:4),yloc(0:4)
      integer blockno
      integer i

      data xloc /0, 1, 1, 0, 0.5d0/
      data yloc /0, 0, 1, 1, 0.5d0/

c     # Map each brick to a [0,1]x[0,1] domain and duplicate
c     # initial conditions.
      do i = 0,4
         rc = sqrt((xc-xloc(i))**2 + (yc-yloc(i))**2)
         fdisc = rc-0.3d0
         if (fdisc .lt. 0) then
            return
         endif
      enddo

      end
